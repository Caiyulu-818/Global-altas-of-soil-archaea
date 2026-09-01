#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Direct-to-GeoTIFF ridge-Mahalanobis OOD layers for Guild 1-5.

This optimized version avoids the slow CSV -> transposed memmap -> CSV workflow.
It reads global environmental CSVs row-by-row (18000 latitude rows x 36000
longitude columns), calculates ridge-Mahalanobis distance, and writes GeoTIFF
rows directly.

It does NOT rerun RF prediction. It reads the existing prediction CSV from
g*_results/pca_multilevel1, uses a temporary raster-oriented memmap, and writes
prediction_ridge_reliable_*.tif directly.

Outputs:
  g*_results/ridge_mahalanobis_ood_tif_direct1/
    mahalanobis_distance(.tif / _abundance.tif)
    environmental_extrapolation_risk_ratio(.tif / _abundance.tif)
    environmental_extrapolation_risk_0_1(.tif / _abundance.tif)
    environmental_extrapolation_risk_class(.tif / _abundance.tif)
    environmental_extrapolation_mask_<strict|moderate|lenient>(.tif / _abundance.tif)
    prediction_ridge_reliable_<strict|moderate|lenient>.tif
    prediction_abundance_ridge_reliable_<strict|moderate|lenient>.tif

Risk class:
  0 = low risk,    MD <= 99% training threshold
  1 = medium risk, 99% < MD <= 99.5% training threshold
  2 = high risk,   MD > 99.5% training threshold
  255 = NoData

Run one guild on a node, e.g.:
  nohup python make_g4_ridge_mahalanobis_ood_direct_tif.py > g4_ridge_direct_tif.log 2>&1 &
"""

import csv
import os
import time
import numpy as np
import pandas as pd
import tempfile

try:
    from osgeo import gdal
except Exception:
    import gdal

gdal.UseExceptions()


# -----------------------------
# 0. User settings
# -----------------------------
BASE_DIR = "/public/home/ylwang/Gisdata/TC_New/lucy"
RESULT_ROOT = os.path.join(BASE_DIR, "uncertainty")
ENV_CSV_ROOT = "/public/home/ZhangFZ/Gisdata/resample/csv"
TIF_ROOT = "/public/home/ZhangFZ/Gisdata/resample/tif"
REFERENCE_TIF = os.path.join(TIF_ROOT, "New_TotalNumber_Bootstrap_CoefVar.tif")
LAND_MASK_TIF = REFERENCE_TIF

MODELS_TO_RUN = ["g1", "g1_abundance", "g2", "g3", "g3_abundance", "g4", "g5"]

# CPU / GDAL settings.
# This mainly accelerates GeoTIFF compression/writing. CSV parsing remains
# mostly single-process I/O, so the biggest speed-up is still running different
# guild scripts on different nodes.
CPU_THREADS = 8
GDAL_CACHEMAX_MB = 2048

gdal.SetConfigOption("GDAL_NUM_THREADS", str(CPU_THREADS))
gdal.SetConfigOption("GDAL_CACHEMAX", str(GDAL_CACHEMAX_MB))

RIDGE = 1e-3
PROFILE_PCTS = {
    "strict": 97.5,
    "moderate": 99.0,
    "lenient": 99.5,
}

LOW_PROFILE = "moderate"
HIGH_PROFILE = "lenient"

FLOAT_NODATA = -9999.0
BYTE_NODATA = 255
PROGRESS_EVERY = 200

# If disk is tight, set WRITE_MAHALANOBIS_DISTANCE = False and keep only risk layers.
WRITE_MAHALANOBIS_DISTANCE = True

# If True, also write prediction_ridge_reliable_*.tif by applying the ridge
# environmental OOD masks to predictions.
WRITE_MASKED_PREDICTION_TIFS = True

# Prediction input is always read from pca_multilevel1/*.csv. prediction.tif is
# deliberately ignored to keep the workflow consistent on the server.
PREDICTION_INPUT_DIRNAME = "pca_multilevel1"

# "finite": land/valid pixels are non-NaN and not explicit NoData.
# "positive": use for binary 1=land, 0=ocean masks.
LAND_MASK_MODE = "finite"


TRAINING_FILES = {
    "g1": os.path.join(BASE_DIR, "env.meta.red1_scaled.csv"),
    "g2": os.path.join(BASE_DIR, "env.meta.red2_scaled.csv"),
    "g3": os.path.join(BASE_DIR, "env.meta.red3_scaled.csv"),
    "g4": os.path.join(BASE_DIR, "env.meta.red4_scaled.csv"),
    "g5": os.path.join(BASE_DIR, "env.meta.red5_scaled.csv"),
}

GLOBAL_FILES = {
    "Elevation": os.path.join(ENV_CSV_ROOT, "New_dem.csv"),
    "evenness": os.path.join(ENV_CSV_ROOT, "New_evenness_01_05_1km_uint16.csv"),
    "northness": os.path.join(ENV_CSV_ROOT, "New_northness.csv"),
    "ref_1": os.path.join(ENV_CSV_ROOT, "New_Nadir_Reflectance_1.csv"),
    "ECE": os.path.join(ENV_CSV_ROOT, "New_Mean_of_annual_predicted_ECe_1980_2018.csv"),
    "biomass": os.path.join(ENV_CSV_ROOT, "New_soilbiomass_1km.csv"),
    "ph": os.path.join(ENV_CSV_ROOT, "New_Average_phh2o2.csv"),
    "soc": os.path.join(ENV_CSV_ROOT, "New_Average_soc2.csv"),
    "ndvimean": os.path.join(ENV_CSV_ROOT, "New_NDVImean.csv"),
}

TIF_CATALOG = {
    "Elevation": os.path.join(TIF_ROOT, "New_dem.tif"),
    "evenness": os.path.join(TIF_ROOT, "New_evenness_01_05_1km_uint16.tif"),
    "northness": os.path.join(TIF_ROOT, "New_northness.tif"),
    "ref_1": os.path.join(TIF_ROOT, "New_Nadir_Reflectance_1.tif"),
    "ECE": os.path.join(TIF_ROOT, "New_Mean_of_annual_predicted_ECe_1980_2018.tif"),
    "biomass": os.path.join(TIF_ROOT, "New_soilbiomass_1km.tif"),
    "ph": os.path.join(TIF_ROOT, "New_Average_phh2o2.tif"),
    "soc": os.path.join(TIF_ROOT, "New_Average_soc2.tif"),
    "ndvimean": os.path.join(TIF_ROOT, "New_NDVImean.tif"),
}

CONFIGS = {
    "g1": {
        "base_guild": "g1",
        "output_subdir": "g1_results",
        "tag": "",
        "features": ["ref_1", "ndvimean", "soc"],
        "prediction_csv": "prediction.csv",
    },
    "g1_abundance": {
        "base_guild": "g1",
        "output_subdir": "g1_results",
        "tag": "_abundance",
        # Final G1 abundance model after removing the incomplete ECE layer.
        "features": ["ref_1", "ndvimean", "Elevation"],
        "prediction_csv": "prediction_abundance.csv",
    },
    "g2": {
        "base_guild": "g2",
        "output_subdir": "g2_results",
        "tag": "",
        "features": ["ph", "ref_1", "ndvimean"],
        "prediction_csv": "prediction.csv",
    },
    "g3": {
        "base_guild": "g3",
        "output_subdir": "g3_results",
        "tag": "",
        "features": ["ref_1", "ndvimean"],
        "prediction_csv": "prediction.csv",
    },
    "g3_abundance": {
        "base_guild": "g3",
        "output_subdir": "g3_results",
        "tag": "_abundance",
        "features": ["ndvimean", "ref_1"],
        "prediction_csv": "prediction_abundance.csv",
    },
    "g4": {
        "base_guild": "g4",
        "output_subdir": "g4_results",
        "tag": "",
        "features": ["ph", "soc", "ref_1", "biomass"],
        "prediction_csv": "prediction.csv",
    },
    "g5": {
        "base_guild": "g5",
        "output_subdir": "g5_results",
        "tag": "",
        "features": ["ref_1", "ndvimean", "biomass", "ph"],
        "prediction_csv": "prediction.csv",
    },
}


# -----------------------------
# 1. Helpers
# -----------------------------
def tagged_name(prefix, tag):
    return f"{prefix}{tag}.tif"


def tagged_profile_name(prefix, tag, profile):
    return f"{prefix}{tag}_{profile}.tif"


def masked_prediction_tif_name(tag, profile):
    if tag:
        return f"prediction{tag}_ridge_reliable_{profile}.tif"
    return f"prediction_ridge_reliable_{profile}.tif"


def prediction_csv_path(cfg):
    return os.path.join(
        RESULT_ROOT,
        cfg["output_subdir"],
        PREDICTION_INPUT_DIRNAME,
        cfg["prediction_csv"],
    )


def open_ds(path):
    ds = gdal.Open(path)
    if ds is None:
        raise FileNotFoundError(path)
    return ds


def parse_csv_row(line, expected_cols):
    line = line.strip()
    if not line:
        return None
    line = line.replace(b"NA", b"nan")
    arr = np.fromstring(line.decode("ascii"), sep=",", dtype=np.float32)
    if arr.size != expected_cols:
        raise ValueError(f"Expected {expected_cols} columns, got {arr.size}")
    return arr


def count_nonempty_rows(path):
    n = 0
    with open(path, "rb") as f:
        for line in f:
            if line.strip():
                n += 1
    return n


def build_prediction_memmap_from_csv(pred_csv, xsize, ysize, tmp_dir, model_name):
    """Read prediction CSV orientation 36000x18000 into raster-row memmap 18000x36000."""
    if not os.path.exists(pred_csv):
        print(f"  [prediction CSV skip] missing: {pred_csv}", flush=True)
        return None, None

    n_rows = count_nonempty_rows(pred_csv)
    if n_rows != xsize:
        print(f"  [prediction CSV skip] row count {n_rows} != {xsize}: {pred_csv}", flush=True)
        return None, None

    os.makedirs(tmp_dir, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(prefix=f"{model_name}_prediction_raster_", suffix=".dat", dir=tmp_dir)
    os.close(fd)
    mm = np.memmap(tmp_path, dtype=np.float32, mode="w+", shape=(ysize, xsize))

    print(f"  prediction CSV -> temporary raster memmap: {pred_csv}", flush=True)
    start = time.time()
    lon_idx = 0
    with open(pred_csv, "rb") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            vals = parse_csv_row(line, ysize)
            mm[:, lon_idx] = vals
            lon_idx += 1
            if lon_idx % 1000 == 0:
                elapsed_h = (time.time() - start) / 3600.0
                done = lon_idx / xsize
                eta_h = elapsed_h / max(done, 1e-9) - elapsed_h
                print(
                    f"    prediction memmap columns {lon_idx}/{xsize} "
                    f"| elapsed h {elapsed_h:.2f} | ETA h {eta_h:.2f}",
                    flush=True,
                )

    if lon_idx != xsize:
        del mm
        try:
            os.remove(tmp_path)
        except OSError:
            pass
        raise RuntimeError(f"{model_name}: prediction CSV columns {lon_idx} != {xsize}")

    mm.flush()
    return mm, tmp_path


def load_training_params(model_name, cfg):
    training_file = TRAINING_FILES[cfg["base_guild"]]
    features = cfg["features"]
    if not os.path.exists(training_file):
        raise FileNotFoundError(training_file)

    train = pd.read_csv(training_file)
    missing = [x for x in features if x not in train.columns]
    if missing:
        raise RuntimeError(f"Missing columns in {training_file}: {missing}")

    X = train[features].apply(pd.to_numeric, errors="coerce")
    med = X.median(axis=0, skipna=True)
    X = X.fillna(med)
    mu = X.mean(axis=0)
    sd = X.std(axis=0, ddof=1).replace(0, 1.0)
    Z = ((X - mu) / sd).to_numpy(dtype=np.float64)

    cov = np.cov(Z, rowvar=False)
    cov = np.atleast_2d(cov) + RIDGE * np.eye(len(features))
    inv_cov = np.linalg.inv(cov)
    d_train = np.sqrt(np.einsum("ij,jk,ik->i", Z, inv_cov, Z))
    thresholds = {name: float(np.percentile(d_train, pct)) for name, pct in PROFILE_PCTS.items()}

    return {
        "model": model_name,
        "training_file": training_file,
        "features": features,
        "mu": mu.to_numpy(dtype=np.float64),
        "sd": sd.to_numpy(dtype=np.float64),
        "inv_cov": inv_cov,
        "thresholds": thresholds,
    }


def write_threshold_csv(path, params, tag):
    rows = []
    for profile, pct in PROFILE_PCTS.items():
        rows.append({
            "model": params["model"],
            "tag": tag,
            "threshold_type": "training_empirical",
            "profile": profile,
            "percentile": pct,
            "threshold": params["thresholds"][profile],
            "ridge": RIDGE,
            "space": "standardized_raw_predictor_space",
            "features": ";".join(params["features"]),
        })
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def open_output_tif(path, ref_ds, dtype, nodata):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    driver = gdal.GetDriverByName("GTiff")
    out = driver.Create(
        path,
        ref_ds.RasterXSize,
        ref_ds.RasterYSize,
        1,
        dtype,
        options=[
            "TILED=YES",
            "BLOCKXSIZE=256",
            "BLOCKYSIZE=256",
            "COMPRESS=LZW",
            f"NUM_THREADS={CPU_THREADS}",
            "BIGTIFF=YES",
        ],
    )
    out.SetGeoTransform(ref_ds.GetGeoTransform())
    out.SetProjection(ref_ds.GetProjection())
    band = out.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    return out, band


def valid_row_from_band(band, y, xsize, top_left=None):
    arr = band.ReadAsArray(0, y, xsize, 1).reshape(-1)
    nodata = band.GetNoDataValue()
    valid = np.isfinite(arr)
    if nodata is not None and np.isfinite(nodata):
        valid &= arr != nodata
    if top_left is not None and np.isfinite(top_left):
        valid &= arr != top_left
    return valid


def open_feature_tif_info(features):
    out = {}
    for feature in features:
        path = TIF_CATALOG[feature]
        if not os.path.exists(path):
            raise FileNotFoundError(f"Missing feature TIF for mask: {feature}: {path}")
        ds = open_ds(path)
        band = ds.GetRasterBand(1)
        top_left = band.ReadAsArray(0, 0, 1, 1)[0, 0]
        out[feature] = (ds, band, top_left)
    return out


def close_feature_tif_info(info):
    for feature, (ds, band, top_left) in info.items():
        band = None
        ds = None


def make_land_valid_row(land_band, y, xsize):
    arr = land_band.ReadAsArray(0, y, xsize, 1).reshape(-1)
    nodata = land_band.GetNoDataValue()
    if LAND_MASK_MODE == "finite":
        valid = np.isfinite(arr)
        if nodata is not None and np.isfinite(nodata):
            valid &= arr != nodata
        return valid
    if LAND_MASK_MODE == "positive":
        valid = np.isfinite(arr) & (arr > 0)
        if nodata is not None and np.isfinite(nodata):
            valid &= arr != nodata
        return valid
    raise ValueError("LAND_MASK_MODE must be 'finite' or 'positive'")


def ridge_md_for_row(feature_arrays, params):
    X = np.vstack(feature_arrays).T.astype(np.float64)
    valid = np.all(np.isfinite(X), axis=1)
    md = np.full(X.shape[0], np.nan, dtype=np.float32)
    if np.any(valid):
        Z = (X[valid, :] - params["mu"]) / params["sd"]
        d = np.sqrt(np.einsum("ij,jk,ik->i", Z, params["inv_cov"], Z))
        md[valid] = d.astype(np.float32)
    return md, valid


# -----------------------------
# 2. Main processing
# -----------------------------
def process_model(model_name, ref_ds, land_band):
    cfg = CONFIGS[model_name]
    tag = cfg["tag"]
    features = cfg["features"]
    output_dir = os.path.join(RESULT_ROOT, cfg["output_subdir"], "ridge_mahalanobis_ood_tif_direct2")
    os.makedirs(output_dir, exist_ok=True)

    print(f"\n=== {model_name} direct ridge-Mahalanobis OOD TIF ===", flush=True)
    print("Training file:", TRAINING_FILES[cfg["base_guild"]], flush=True)
    print("Features:", ", ".join(features), flush=True)
    print("Output dir:", output_dir, flush=True)

    params = load_training_params(model_name, cfg)
    threshold_csv = os.path.join(output_dir, f"{model_name}_ridge_mahalanobis_thresholds.csv")
    write_threshold_csv(threshold_csv, params, tag)
    for profile, threshold in params["thresholds"].items():
        print(f"  {profile}: p{PROFILE_PCTS[profile]} threshold = {threshold:.6g}", flush=True)

    xsize = ref_ds.RasterXSize
    ysize = ref_ds.RasterYSize

    global_handles = []
    feature_tifs = {}
    outputs = {}
    pred_mm = None
    pred_mm_path = None
    try:
        for feature in features:
            csv_path = GLOBAL_FILES[feature]
            if not os.path.exists(csv_path):
                raise FileNotFoundError(f"Missing global CSV for {feature}: {csv_path}")
            print(f"  reading global CSV: {feature} <- {csv_path}", flush=True)
            global_handles.append(open(csv_path, "rb"))

        feature_tifs = open_feature_tif_info(features)

        use_prediction_input = False
        pred_csv = prediction_csv_path(cfg)
        if WRITE_MASKED_PREDICTION_TIFS:
            pred_mm, pred_mm_path = build_prediction_memmap_from_csv(
                pred_csv=pred_csv,
                xsize=xsize,
                ysize=ysize,
                tmp_dir=os.path.join(output_dir, "_tmp_prediction_memmap"),
                model_name=model_name,
            )
            if pred_mm is not None:
                use_prediction_input = True
                print(f"  using prediction CSV for masked outputs: {pred_csv}", flush=True)

        if WRITE_MAHALANOBIS_DISTANCE:
            outputs["md"] = open_output_tif(
                os.path.join(output_dir, tagged_name("mahalanobis_distance", tag)),
                ref_ds, gdal.GDT_Float32, FLOAT_NODATA
            )
        outputs["ratio"] = open_output_tif(
            os.path.join(output_dir, tagged_name("environmental_extrapolation_risk_ratio", tag)),
            ref_ds, gdal.GDT_Float32, FLOAT_NODATA
        )
        outputs["risk01"] = open_output_tif(
            os.path.join(output_dir, tagged_name("environmental_extrapolation_risk_0_1", tag)),
            ref_ds, gdal.GDT_Float32, FLOAT_NODATA
        )
        outputs["class"] = open_output_tif(
            os.path.join(output_dir, tagged_name("environmental_extrapolation_risk_class", tag)),
            ref_ds, gdal.GDT_Byte, BYTE_NODATA
        )
        for profile in PROFILE_PCTS:
            outputs[f"mask_{profile}"] = open_output_tif(
                os.path.join(output_dir, tagged_profile_name("environmental_extrapolation_mask", tag, profile)),
                ref_ds, gdal.GDT_Byte, BYTE_NODATA
            )
            if use_prediction_input:
                outputs[f"pred_{profile}"] = open_output_tif(
                    os.path.join(output_dir, masked_prediction_tif_name(tag, profile)),
                    ref_ds, gdal.GDT_Float32, FLOAT_NODATA
                )

        class_counts = {0: 0, 1: 0, 2: 0}
        keep_counts = {profile: 0 for profile in PROFILE_PCTS}
        valid_total = 0
        start = time.time()

        for y in range(ysize):
            feature_arrays = []
            for handle, feature in zip(global_handles, features):
                line = handle.readline()
                if not line:
                    raise RuntimeError(f"{model_name}: {feature} global CSV ended at row {y}")
                feature_arrays.append(parse_csv_row(line, xsize))

            md, valid = ridge_md_for_row(feature_arrays, params)

            # Land + model-feature valid mask from TIFs.
            row_valid = valid & make_land_valid_row(land_band, y, xsize)
            for feature in features:
                ds, band, top_left = feature_tifs[feature]
                row_valid &= valid_row_from_band(band, y, xsize, top_left=top_left)

            t_low = params["thresholds"][LOW_PROFILE]
            t_high = params["thresholds"][HIGH_PROFILE]

            ratio = np.full(xsize, FLOAT_NODATA, dtype=np.float32)
            risk01 = np.full(xsize, FLOAT_NODATA, dtype=np.float32)
            klass = np.full(xsize, BYTE_NODATA, dtype=np.uint8)

            ratio[row_valid] = md[row_valid] / t_low
            risk01[row_valid] = np.minimum(md[row_valid] / t_high, 1.0)
            klass[row_valid & (md <= t_low)] = 0
            klass[row_valid & (md > t_low) & (md <= t_high)] = 1
            klass[row_valid & (md > t_high)] = 2

            valid_total += int(np.sum(row_valid))
            class_counts[0] += int(np.sum(klass == 0))
            class_counts[1] += int(np.sum(klass == 1))
            class_counts[2] += int(np.sum(klass == 2))

            if WRITE_MAHALANOBIS_DISTANCE:
                md_out = np.full(xsize, FLOAT_NODATA, dtype=np.float32)
                md_out[row_valid] = md[row_valid]
                outputs["md"][1].WriteArray(md_out.reshape(1, -1), 0, y)
            outputs["ratio"][1].WriteArray(ratio.reshape(1, -1), 0, y)
            outputs["risk01"][1].WriteArray(risk01.reshape(1, -1), 0, y)
            outputs["class"][1].WriteArray(klass.reshape(1, -1), 0, y)

            pred_row = None
            pred_valid = None
            if use_prediction_input:
                pred_row = np.asarray(pred_mm[y, :], dtype=np.float32)
                pred_valid = np.isfinite(pred_row)

            for profile, threshold in params["thresholds"].items():
                mask = np.full(xsize, BYTE_NODATA, dtype=np.uint8)
                keep = row_valid & (md <= threshold)
                mask[row_valid] = np.where(keep[row_valid], 0, 1).astype(np.uint8)
                keep_counts[profile] += int(np.sum(keep))
                outputs[f"mask_{profile}"][1].WriteArray(mask.reshape(1, -1), 0, y)
                if use_prediction_input:
                    pred_out = np.full(xsize, FLOAT_NODATA, dtype=np.float32)
                    pred_keep = keep & pred_valid
                    pred_out[pred_keep] = pred_row[pred_keep]
                    outputs[f"pred_{profile}"][1].WriteArray(pred_out.reshape(1, -1), 0, y)

            if (y + 1) % PROGRESS_EVERY == 0 or (y + 1) == ysize:
                elapsed_h = (time.time() - start) / 3600.0
                done = (y + 1) / ysize
                eta_h = elapsed_h / max(done, 1e-9) - elapsed_h
                print(
                    f"[{model_name}] raster rows {y + 1}/{ysize} "
                    f"| elapsed h {elapsed_h:.2f} | ETA h {eta_h:.2f}",
                    flush=True,
                )

        for key, (ds, band) in outputs.items():
            band.FlushCache()
            ds.FlushCache()
            outputs[key] = (None, None)

        summary_rows = []
        for cls, count in class_counts.items():
            summary_rows.append({
                "model": model_name,
                "metric": f"risk_class_{cls}_fraction",
                "value": count / valid_total if valid_total else np.nan,
            })
        for profile, count in keep_counts.items():
            summary_rows.append({
                "model": model_name,
                "metric": f"{profile}_in_distribution_fraction",
                "value": count / valid_total if valid_total else np.nan,
            })
        pd.DataFrame(summary_rows).to_csv(
            os.path.join(output_dir, f"{model_name}_ridge_ood_pixel_summary.csv"),
            index=False
        )

        print(f"[{model_name}] done raster rows={ysize}", flush=True)
        for row in summary_rows:
            print(f"  {row['metric']} = {row['value']:.4f}", flush=True)

    finally:
        for handle in global_handles:
            try:
                handle.close()
            except Exception:
                pass
        close_feature_tif_info(feature_tifs)
        if pred_mm is not None:
            del pred_mm
        if pred_mm_path is not None:
            try:
                os.remove(pred_mm_path)
            except OSError:
                pass
        for key, value in list(outputs.items()):
            ds, band = value
            band = None
            ds = None


def main():
    ref_ds = open_ds(REFERENCE_TIF)
    land_ds = open_ds(LAND_MASK_TIF)
    land_band = land_ds.GetRasterBand(1)

    print("Reference grid:", REFERENCE_TIF, flush=True)
    print("Raster size:", ref_ds.RasterXSize, ref_ds.RasterYSize, flush=True)
    for model_name in MODELS_TO_RUN:
        process_model(model_name, ref_ds, land_band)

    land_band = None
    land_ds = None
    ref_ds = None
    print("\nAll direct ridge-Mahalanobis OOD TIF layers finished.", flush=True)


if __name__ == "__main__":
    main()
