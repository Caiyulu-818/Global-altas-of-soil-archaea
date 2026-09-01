#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
TD/GFM scenario ridge-Mahalanobis OOD layers + ridge-masked prediction TIFs.

This is a post-processing script. It does NOT rerun TD/GFM RF prediction.
It reads:
  1) training metadata for TD or GFM,
  2) static global environmental CSV layers,
  3) the scenario-specific cropland CSV layer,
  4) the already completed scenario prediction.csv.

It writes direct GeoTIFF outputs without creating huge intermediate CSVs.

Output directory:
  /public/home/ylwang/Gisdata/TC_New/lucy/cropland/scenario_extrapolation/
    <TD|GFM>/<scenario>/ridge_mahalanobis_ood_tif_direct/

Main outputs:
  mahalanobis_distance.tif
  environmental_extrapolation_risk_ratio.tif
  environmental_extrapolation_risk_0_1.tif
  environmental_extrapolation_risk_class.tif
  environmental_extrapolation_mask_<strict|moderate|lenient>.tif
  prediction_ridge_reliable_<strict|moderate|lenient>.tif

Risk class:
  0 = low risk,    MD <= 99% training threshold
  1 = medium risk, 99% < MD <= 99.5% training threshold
  2 = high risk,   MD > 99.5% training threshold
  255 = NoData

Run examples:
  python td_gfm_ridge_mahalanobis_ood_direct_tif.py TD cropland2020_SSP245
  python td_gfm_ridge_mahalanobis_ood_direct_tif.py GFM cropland2100_SSP370

The wrapper scripts set DATASET_OVERRIDE and SCENARIO_OVERRIDE automatically.
"""

import csv
import os
import sys
import time
import tempfile
import numpy as np
import pandas as pd

try:
    from osgeo import gdal
except Exception:
    import gdal

gdal.UseExceptions()


# -----------------------------
# 0. User settings
# -----------------------------
BASE_DIR = "/public/home/ylwang/Gisdata/TC_New/lucy"
TRAINING_DIR = os.path.join(BASE_DIR, "cropland")
ENV_CSV_ROOT = "/public/home/ZhangFZ/Gisdata/resample/csv"
TIF_ROOT = "/public/home/ZhangFZ/Gisdata/resample/tif"
OUT_ROOT = os.path.join(TRAINING_DIR, "scenario_extrapolation")

REFERENCE_TIF = os.path.join(TIF_ROOT, "New_TotalNumber_Bootstrap_CoefVar.tif")
LAND_MASK_TIF = REFERENCE_TIF

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
ROW_BLOCK_SIZE = 256
PRINT_EVERY_BLOCK = True
LAND_MASK_MODE = "finite"

WRITE_MAHALANOBIS_DISTANCE = True
WRITE_MASKED_PREDICTIONS = True

# Fast prediction.csv fallback.
# Fortran-order memmap makes each CSV longitude row a contiguous write.
PREDICTION_MEMMAP_ORDER = "F"
CHECK_PREDICTION_ROW_COUNT = False

# Prefer node-local scratch for the temporary 18000 x 36000 prediction memmap.
LOCAL_TMP_ROOT = os.environ.get("SLURM_TMPDIR") or os.environ.get("TMPDIR") or None


DATASET_CONFIGS = {
    "TD": {
        "training_file": os.path.join(TRAINING_DIR, "TD.env.meta.csv"),
        "target": "Shannon.TD",
        "features_static": [
            "New_NDVImean",
            "New_wc2.1_30s_bio_12",
            "aridity",
            "Soil.pH.x.10.in.H2O",
            "New_wc2.1_30s_bio_1",
        ],
        "training_cropland_col": "cropland2020_SSP245",
    },
    "GFM": {
        "training_file": os.path.join(TRAINING_DIR, "GFM.env.meta.with_lonlat.csv"),
        "target": "GFM",
        "features_static": [
            "aridity",
            "New_NDVImean",
            "New_northness",
        ],
        "training_cropland_col": "cropland2020_SSP245",
    },
}

GLOBAL_FILE_CATALOG = {
    "aridity": os.path.join(ENV_CSV_ROOT, "New_ai_et0.csv"),
    "Soil.pH.x.10.in.H2O": os.path.join(ENV_CSV_ROOT, "New_Average_phh2o2.csv"),
    "New_northness": os.path.join(ENV_CSV_ROOT, "New_northness.csv"),
    "New_NDVImean": os.path.join(ENV_CSV_ROOT, "New_NDVImean.csv"),
    "New_wc2.1_30s_bio_1": os.path.join(ENV_CSV_ROOT, "New_wc2.1_30s_bio_1.csv"),
    "New_wc2.1_30s_bio_12": os.path.join(ENV_CSV_ROOT, "New_wc2.1_30s_bio_12.csv"),
}

CROPLAND_SCENARIOS = {
    "cropland2100_SSP126": os.path.join(
        TRAINING_DIR,
        "New_globalCropland_2100CE_SSP126_LandOnly.csv",
        "New_globalCropland_2100CE_SSP126_LandOnly.csv",
    ),
    "cropland2100_SSP370": os.path.join(
        TRAINING_DIR,
        "New_globalCropland_2100CE_SSP370_LandOnly.csv",
        "New_globalCropland_2100CE_SSP370_LandOnly.csv",
    ),
    "cropland2020_SSP245": os.path.join(
        TRAINING_DIR,
        "New_globalCropland_2020CE_SSP245_LandOnly.csv",
        "New_globalCropland_2020CE_SSP245_LandOnly.csv",
    ),
    "cropland2100_SSP245": os.path.join(
        TRAINING_DIR,
        "New_globalCropland_2100CE_SSP245_LandOnly.csv",
        "New_globalCropland_2100CE_SSP245_LandOnly.csv",
    ),
}


# Optional wrapper overrides.
DATASET_OVERRIDE = globals().get("DATASET_OVERRIDE", None)
SCENARIO_OVERRIDE = globals().get("SCENARIO_OVERRIDE", None)


# -----------------------------
# 1. Helpers
# -----------------------------
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


def choose_tmp_dir(tmp_dir, label):
    if LOCAL_TMP_ROOT:
        return os.path.join(LOCAL_TMP_ROOT, f"ridge_{label}_prediction_memmap")
    return tmp_dir


def open_output_tif(path, ref_ds, dtype, nodata):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    driver = gdal.GetDriverByName("GTiff")
    out = driver.Create(
        path,
        ref_ds.RasterXSize,
        ref_ds.RasterYSize,
        1,
        dtype,
        options=["TILED=YES", "COMPRESS=LZW", "BIGTIFF=YES"],
    )
    out.SetGeoTransform(ref_ds.GetGeoTransform())
    out.SetProjection(ref_ds.GetProjection())
    band = out.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    return out, band


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


def build_prediction_memmap_from_csv(pred_csv, xsize, ysize, tmp_dir, label):
    if not os.path.exists(pred_csv):
        raise FileNotFoundError(f"Missing prediction CSV: {pred_csv}")

    if CHECK_PREDICTION_ROW_COUNT:
        n_rows = count_nonempty_rows(pred_csv)
        if n_rows != xsize:
            raise RuntimeError(f"{label}: prediction.csv rows={n_rows}, expected={xsize}: {pred_csv}")

    tmp_dir = choose_tmp_dir(tmp_dir, label)
    os.makedirs(tmp_dir, exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(prefix=f"{label}_prediction_raster_", suffix=".dat", dir=tmp_dir)
    os.close(fd)
    mm = np.memmap(tmp_path, dtype=np.float32, mode="w+", shape=(ysize, xsize), order=PREDICTION_MEMMAP_ORDER)

    print(f"Prediction CSV -> temporary raster memmap: {pred_csv}", flush=True)
    print(f"  memmap path: {tmp_path}", flush=True)
    print(f"  memmap order: {PREDICTION_MEMMAP_ORDER}", flush=True)
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
                    f"  prediction columns {lon_idx}/{xsize} | elapsed h {elapsed_h:.2f} | ETA h {eta_h:.2f}",
                    flush=True,
                )

    if lon_idx != xsize:
        del mm
        try:
            os.remove(tmp_path)
        except OSError:
            pass
        raise RuntimeError(f"{label}: prediction CSV columns={lon_idx}, expected={xsize}")

    mm.flush()
    return mm, tmp_path


def load_training_params(dataset):
    cfg = DATASET_CONFIGS[dataset]
    features = cfg["features_static"] + ["cropland_scenario"]
    train = pd.read_csv(cfg["training_file"])
    missing = [x for x in cfg["features_static"] + [cfg["training_cropland_col"], cfg["target"]] if x not in train.columns]
    if missing:
        raise RuntimeError(f"Missing training columns in {cfg['training_file']}: {missing}")

    train = train.copy()
    train["cropland_scenario"] = pd.to_numeric(train[cfg["training_cropland_col"]], errors="coerce")
    X = train[features].apply(pd.to_numeric, errors="coerce")
    y = pd.to_numeric(train[cfg["target"]], errors="coerce")
    ok = X.notna().all(axis=1) & np.isfinite(y)
    X = X.loc[ok]

    if len(X) < len(features) + 2:
        raise RuntimeError(f"{dataset}: too few complete training rows for ridge-Mahalanobis.")

    mu = X.mean(axis=0)
    sd = X.std(axis=0, ddof=1).replace(0, 1.0)
    Z = ((X - mu) / sd).to_numpy(dtype=np.float64)

    cov = np.cov(Z, rowvar=False)
    cov = np.atleast_2d(cov) + RIDGE * np.eye(len(features))
    inv_cov = np.linalg.inv(cov)
    d_train = np.sqrt(np.einsum("ij,jk,ik->i", Z, inv_cov, Z))
    thresholds = {name: float(np.percentile(d_train, pct)) for name, pct in PROFILE_PCTS.items()}

    return {
        "features": features,
        "features_static": cfg["features_static"],
        "mu": mu.to_numpy(dtype=np.float64),
        "sd": sd.to_numpy(dtype=np.float64),
        "inv_cov": inv_cov,
        "thresholds": thresholds,
        "n_training_rows": int(len(X)),
    }


def write_threshold_csv(path, dataset, scenario, params):
    rows = []
    for profile, pct in PROFILE_PCTS.items():
        rows.append({
            "dataset": dataset,
            "scenario": scenario,
            "threshold_type": "training_empirical",
            "profile": profile,
            "percentile": pct,
            "threshold": params["thresholds"][profile],
            "ridge": RIDGE,
            "space": "standardized_raw_predictor_space",
            "features": ";".join(params["features"]),
            "n_training_rows": params["n_training_rows"],
        })
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def ridge_md_for_row(feature_arrays, params):
    X = np.vstack(feature_arrays).T.astype(np.float64)
    valid = np.all(np.isfinite(X), axis=1)
    md = np.full(X.shape[0], np.nan, dtype=np.float32)
    if np.any(valid):
        Z = (X[valid, :] - params["mu"]) / params["sd"]
        d = np.sqrt(np.einsum("ij,jk,ik->i", Z, params["inv_cov"], Z))
        md[valid] = d.astype(np.float32)
    return md, valid


def run_one(dataset, scenario):
    dataset = dataset.upper()
    if dataset not in DATASET_CONFIGS:
        raise RuntimeError(f"dataset must be TD or GFM, got {dataset}")
    if scenario not in CROPLAND_SCENARIOS:
        raise RuntimeError(f"Unknown scenario: {scenario}")

    cfg = DATASET_CONFIGS[dataset]
    scenario_dir = os.path.join(OUT_ROOT, dataset, scenario)
    prediction_csv = os.path.join(scenario_dir, "prediction.csv")
    output_dir = os.path.join(scenario_dir, "ridge_mahalanobis_ood_tif_direct")
    os.makedirs(output_dir, exist_ok=True)

    ref_ds = open_ds(REFERENCE_TIF)
    land_ds = open_ds(LAND_MASK_TIF)
    land_band = land_ds.GetRasterBand(1)
    xsize = ref_ds.RasterXSize
    ysize = ref_ds.RasterYSize

    label = f"{dataset}_{scenario}"
    print("\n==============================", flush=True)
    print("Ridge-Mahalanobis OOD:", label, flush=True)
    print("==============================", flush=True)
    print("Prediction CSV:", prediction_csv, flush=True)
    print("Output dir:", output_dir, flush=True)

    params = load_training_params(dataset)
    threshold_csv = os.path.join(output_dir, f"{label}_ridge_mahalanobis_thresholds.csv")
    write_threshold_csv(threshold_csv, dataset, scenario, params)
    for profile, threshold in params["thresholds"].items():
        print(f"  {profile}: p{PROFILE_PCTS[profile]} threshold = {threshold:.6g}", flush=True)

    pred_mm = None
    pred_mm_path = None
    global_handles = []
    outputs = {}

    try:
        pred_mm, pred_mm_path = build_prediction_memmap_from_csv(
            pred_csv=prediction_csv,
            xsize=xsize,
            ysize=ysize,
            tmp_dir=os.path.join(output_dir, "_tmp_prediction_memmap"),
            label=label,
        )

        for feature in cfg["features_static"]:
            path = GLOBAL_FILE_CATALOG[feature]
            if not os.path.exists(path):
                raise FileNotFoundError(f"Missing global CSV for {feature}: {path}")
            print(f"Reading static CSV: {feature} <- {path}", flush=True)
            global_handles.append((feature, open(path, "rb")))
        crop_path = CROPLAND_SCENARIOS[scenario]
        if not os.path.exists(crop_path):
            raise FileNotFoundError(f"Missing cropland CSV: {crop_path}")
        print(f"Reading cropland CSV: cropland_scenario <- {crop_path}", flush=True)
        global_handles.append(("cropland_scenario", open(crop_path, "rb")))

        if WRITE_MAHALANOBIS_DISTANCE:
            outputs["md"] = open_output_tif(
                os.path.join(output_dir, "mahalanobis_distance.tif"),
                ref_ds, gdal.GDT_Float32, FLOAT_NODATA
            )
        outputs["ratio"] = open_output_tif(
            os.path.join(output_dir, "environmental_extrapolation_risk_ratio.tif"),
            ref_ds, gdal.GDT_Float32, FLOAT_NODATA
        )
        outputs["risk01"] = open_output_tif(
            os.path.join(output_dir, "environmental_extrapolation_risk_0_1.tif"),
            ref_ds, gdal.GDT_Float32, FLOAT_NODATA
        )
        outputs["class"] = open_output_tif(
            os.path.join(output_dir, "environmental_extrapolation_risk_class.tif"),
            ref_ds, gdal.GDT_Byte, BYTE_NODATA
        )
        for profile in PROFILE_PCTS:
            outputs[f"mask_{profile}"] = open_output_tif(
                os.path.join(output_dir, f"environmental_extrapolation_mask_{profile}.tif"),
                ref_ds, gdal.GDT_Byte, BYTE_NODATA
            )
            if WRITE_MASKED_PREDICTIONS:
                outputs[f"pred_{profile}"] = open_output_tif(
                    os.path.join(output_dir, f"prediction_ridge_reliable_{profile}.tif"),
                    ref_ds, gdal.GDT_Float32, FLOAT_NODATA
                )

        t_low = params["thresholds"][LOW_PROFILE]
        t_high = params["thresholds"][HIGH_PROFILE]
        class_counts = {0: 0, 1: 0, 2: 0}
        keep_counts = {profile: 0 for profile in PROFILE_PCTS}
        valid_total = 0
        start = time.time()

        for y0 in range(0, ysize, ROW_BLOCK_SIZE):
            nrows = min(ROW_BLOCK_SIZE, ysize - y0)
            y1 = y0 + nrows

            md_block = np.full((nrows, xsize), np.nan, dtype=np.float32)
            valid_block = np.zeros((nrows, xsize), dtype=bool)

            for r in range(nrows):
                y = y0 + r
                row_by_feature = {}
                for feature, handle in global_handles:
                    line = handle.readline()
                    if not line:
                        raise RuntimeError(f"{label}: {feature} CSV ended at raster row {y}")
                    row_by_feature[feature] = parse_csv_row(line, xsize)

                feature_arrays = [row_by_feature[f] for f in params["features"]]
                md, valid = ridge_md_for_row(feature_arrays, params)
                row_valid = valid & make_land_valid_row(land_band, y, xsize)
                md_block[r, :] = md
                valid_block[r, :] = row_valid

            ratio = np.full((nrows, xsize), FLOAT_NODATA, dtype=np.float32)
            risk01 = np.full((nrows, xsize), FLOAT_NODATA, dtype=np.float32)
            klass = np.full((nrows, xsize), BYTE_NODATA, dtype=np.uint8)
            ratio[valid_block] = md_block[valid_block] / t_low
            risk01[valid_block] = np.minimum(md_block[valid_block] / t_high, 1.0)
            klass[valid_block & (md_block <= t_low)] = 0
            klass[valid_block & (md_block > t_low) & (md_block <= t_high)] = 1
            klass[valid_block & (md_block > t_high)] = 2

            valid_total += int(np.sum(valid_block))
            class_counts[0] += int(np.sum(klass == 0))
            class_counts[1] += int(np.sum(klass == 1))
            class_counts[2] += int(np.sum(klass == 2))

            if WRITE_MAHALANOBIS_DISTANCE:
                md_out = np.full((nrows, xsize), FLOAT_NODATA, dtype=np.float32)
                md_out[valid_block] = md_block[valid_block]
                outputs["md"][1].WriteArray(md_out, 0, y0)
            outputs["ratio"][1].WriteArray(ratio, 0, y0)
            outputs["risk01"][1].WriteArray(risk01, 0, y0)
            outputs["class"][1].WriteArray(klass, 0, y0)

            pred_block = np.array(pred_mm[y0:y1, :], dtype=np.float32, order="C", copy=True)
            pred_valid = np.isfinite(pred_block) & (pred_block != FLOAT_NODATA)

            for profile, threshold in params["thresholds"].items():
                keep = valid_block & (md_block <= threshold)
                mask = np.full((nrows, xsize), BYTE_NODATA, dtype=np.uint8)
                mask[valid_block] = np.where(keep[valid_block], 0, 1).astype(np.uint8)
                keep_counts[profile] += int(np.sum(keep))
                outputs[f"mask_{profile}"][1].WriteArray(mask, 0, y0)

                if WRITE_MASKED_PREDICTIONS:
                    pred_out = np.full((nrows, xsize), FLOAT_NODATA, dtype=np.float32)
                    pred_keep = keep & pred_valid
                    pred_out[pred_keep] = pred_block[pred_keep]
                    outputs[f"pred_{profile}"][1].WriteArray(pred_out, 0, y0)

            if PRINT_EVERY_BLOCK or y1 % PROGRESS_EVERY == 0 or y1 == ysize:
                elapsed_h = (time.time() - start) / 3600.0
                done = y1 / ysize
                eta_h = elapsed_h / max(done, 1e-9) - elapsed_h
                print(
                    f"{label} | raster rows {y1}/{ysize} "
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
                "dataset": dataset,
                "scenario": scenario,
                "metric": f"risk_class_{cls}_fraction",
                "value": count / valid_total if valid_total else np.nan,
            })
        for profile, count in keep_counts.items():
            summary_rows.append({
                "dataset": dataset,
                "scenario": scenario,
                "metric": f"{profile}_in_distribution_fraction",
                "value": count / valid_total if valid_total else np.nan,
            })
        pd.DataFrame(summary_rows).to_csv(
            os.path.join(output_dir, f"{label}_ridge_ood_pixel_summary.csv"),
            index=False
        )

        print(f"{label} done raster rows={ysize}", flush=True)
        for row in summary_rows:
            print(f"  {row['metric']} = {row['value']:.4f}", flush=True)

    finally:
        for feature, handle in global_handles:
            try:
                handle.close()
            except Exception:
                pass
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
        land_band = None
        land_ds = None
        ref_ds = None


def main():
    if DATASET_OVERRIDE is not None and SCENARIO_OVERRIDE is not None:
        dataset = DATASET_OVERRIDE
        scenario = SCENARIO_OVERRIDE
    else:
        args = sys.argv[1:]
        if len(args) < 2:
            raise SystemExit("Usage: python td_gfm_ridge_mahalanobis_ood_direct_tif.py <TD|GFM> <scenario>")
        dataset, scenario = args[0], args[1]
    run_one(dataset, scenario)


if __name__ == "__main__":
    main()
