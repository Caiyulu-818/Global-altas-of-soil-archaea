#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Add extra TD/GFM ridge-Mahalanobis 90% OOD profiles from existing outputs.

This script does NOT rerun TD/GFM RF prediction and does NOT recompute global
Mahalanobis distance from environmental CSV layers.

It uses existing files from:

  /public/home/ylwang/Gisdata/TC_New/lucy/cropland/scenario_extrapolation/
    <TD|GFM>/<scenario>/
      prediction.csv
      ridge_mahalanobis_ood_tif_direct/mahalanobis_distance.tif

Outputs are added to the same ridge output directory:

  environmental_extrapolation_mask_p90.tif
  prediction_ridge_reliable_p90.tif

or:

  environmental_extrapolation_mask_chi90.tif
  prediction_ridge_reliable_chi90.tif

Mask convention:
  0   = in-distribution / retained
  1   = environmental outlier / masked
  255 = NoData
"""

import csv
import os
import sys
import tempfile
import time

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
OUT_ROOT = os.path.join(TRAINING_DIR, "scenario_extrapolation")

# Existing ridge output directory used as input for mahalanobis_distance.tif.
RIDGE_INPUT_DIRNAME = "ridge_mahalanobis_ood_tif_direct"

# New 90% profile outputs are written here, leaving previous ridge outputs intact.
RIDGE_OUTPUT_DIRNAME = "ridge_mahalanobis_ood_tif_direct_90_profiles"

RIDGE = 1e-3
FLOAT_NODATA = -9999.0
BYTE_NODATA = 255
PROGRESS_EVERY = 200
ROW_BLOCK_SIZE = 256
PRINT_EVERY_BLOCK = True

# The prediction CSV is written as longitude columns. A Fortran-order memmap
# makes each CSV row a contiguous column write; row blocks are copied to C-order
# before writing TIF rows, so both phases stay fast.
PREDICTION_MEMMAP_ORDER = "F"
CHECK_PREDICTION_ROW_COUNT = False
LOCAL_TMP_ROOT = os.environ.get("SLURM_TMPDIR") or os.environ.get("TMPDIR") or None

# Wrapper scripts set these. Direct command-line usage is also supported.
DATASET_OVERRIDE = globals().get("DATASET_OVERRIDE", None)
SCENARIO_OVERRIDE = globals().get("SCENARIO_OVERRIDE", None)
PROFILE_MODE = globals().get("PROFILE_MODE", "empirical_p90")


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

CROPLAND_SCENARIOS = {
    "cropland2100_SSP126",
    "cropland2100_SSP245",
    "cropland2100_SSP370",
    "cropland2020_SSP245",
}


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
    mm = np.memmap(
        tmp_path,
        dtype=np.float32,
        mode="w+",
        shape=(ysize, xsize),
        order=PREDICTION_MEMMAP_ORDER,
    )

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


def md_valid_row(md, nodata):
    valid = np.isfinite(md)
    if nodata is not None and np.isfinite(nodata):
        valid &= md != nodata
    valid &= md != FLOAT_NODATA
    return valid


def chi_square_quantile_0p90(df):
    table = {
        1: 2.705543454095404,
        2: 4.605170185988092,
        3: 6.251388631170325,
        4: 7.779440339734858,
        5: 9.236356899781123,
        6: 10.644640675668422,
        7: 12.017036623780532,
        8: 13.36156613651173,
        9: 14.683656573259837,
        10: 15.987179172105265,
    }
    if df not in table:
        raise RuntimeError(f"No built-in qchisq(0.90, df) value for df={df}.")
    return table[df]


def training_distances(dataset):
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
    return d_train, features, int(len(X))


def threshold_for(dataset):
    d_train, features, n_training_rows = training_distances(dataset)
    df = len(features)

    if PROFILE_MODE == "empirical_p90":
        return {
            "profile": "p90",
            "threshold_type": "training_empirical",
            "threshold": float(np.percentile(d_train, 90.0)),
            "percentile": 90.0,
            "probability": np.nan,
            "df": df,
            "chi_square_threshold_d2": np.nan,
            "features": features,
            "n_training_rows": n_training_rows,
        }

    if PROFILE_MODE == "chisq90":
        chi2 = chi_square_quantile_0p90(df)
        return {
            "profile": "chi90",
            "threshold_type": "chi_square",
            "threshold": float(np.sqrt(chi2)),
            "percentile": np.nan,
            "probability": 0.90,
            "df": df,
            "chi_square_threshold_d2": chi2,
            "features": features,
            "n_training_rows": n_training_rows,
        }

    raise RuntimeError("PROFILE_MODE must be 'empirical_p90' or 'chisq90'")


def write_threshold_csv(path, dataset, scenario, threshold_info):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    row = {
        "dataset": dataset,
        "scenario": scenario,
        "threshold_type": threshold_info["threshold_type"],
        "profile": threshold_info["profile"],
        "percentile": threshold_info["percentile"],
        "probability": threshold_info["probability"],
        "df": threshold_info["df"],
        "chi_square_threshold_d2": threshold_info["chi_square_threshold_d2"],
        "threshold": threshold_info["threshold"],
        "ridge": RIDGE,
        "space": "standardized_raw_predictor_space_sqrt_distance",
        "features": ";".join(threshold_info["features"]),
        "n_training_rows": threshold_info["n_training_rows"],
    }
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)


def run_one(dataset, scenario):
    dataset = dataset.upper()
    if dataset not in DATASET_CONFIGS:
        raise RuntimeError(f"dataset must be TD or GFM, got {dataset}")
    if scenario not in CROPLAND_SCENARIOS:
        raise RuntimeError(f"Unknown scenario: {scenario}")

    scenario_dir = os.path.join(OUT_ROOT, dataset, scenario)
    prediction_csv = os.path.join(scenario_dir, "prediction.csv")
    ridge_input_dir = os.path.join(scenario_dir, RIDGE_INPUT_DIRNAME)
    ridge_output_dir = os.path.join(scenario_dir, RIDGE_OUTPUT_DIRNAME)
    md_tif = os.path.join(ridge_input_dir, "mahalanobis_distance.tif")

    label = f"{dataset}_{scenario}"
    threshold_info = threshold_for(dataset)
    profile = threshold_info["profile"]

    print("\n==============================", flush=True)
    print(f"Add TD/GFM ridge profile: {label} | {profile}", flush=True)
    print("==============================", flush=True)
    print("Prediction CSV:", prediction_csv, flush=True)
    print("Mahalanobis TIF:", md_tif, flush=True)
    print("Ridge input dir:", ridge_input_dir, flush=True)
    print("Ridge output dir:", ridge_output_dir, flush=True)
    print(f"Threshold ({threshold_info['threshold_type']}): {threshold_info['threshold']:.6g}", flush=True)

    if not os.path.exists(md_tif):
        raise FileNotFoundError(f"Missing existing Mahalanobis TIF: {md_tif}")

    md_ds = open_ds(md_tif)
    md_band = md_ds.GetRasterBand(1)
    md_nodata = md_band.GetNoDataValue()
    xsize = md_ds.RasterXSize
    ysize = md_ds.RasterYSize

    write_threshold_csv(
        os.path.join(ridge_output_dir, f"{label}_ridge_mahalanobis_threshold_{profile}.csv"),
        dataset,
        scenario,
        threshold_info,
    )

    pred_mm = None
    pred_mm_path = None
    mask_ds = None
    mask_band = None
    pred_ds = None
    pred_band = None

    try:
        pred_mm, pred_mm_path = build_prediction_memmap_from_csv(
            pred_csv=prediction_csv,
            xsize=xsize,
            ysize=ysize,
            tmp_dir=os.path.join(ridge_output_dir, "_tmp_prediction_memmap"),
            label=f"{label}_{profile}",
        )

        mask_ds, mask_band = open_output_tif(
            os.path.join(ridge_output_dir, f"environmental_extrapolation_mask_{profile}.tif"),
            md_ds,
            gdal.GDT_Byte,
            BYTE_NODATA,
        )
        pred_ds, pred_band = open_output_tif(
            os.path.join(ridge_output_dir, f"prediction_ridge_reliable_{profile}.tif"),
            md_ds,
            gdal.GDT_Float32,
            FLOAT_NODATA,
        )

        kept = 0
        valid_total = 0
        threshold = threshold_info["threshold"]
        start = time.time()

        for y0 in range(0, ysize, ROW_BLOCK_SIZE):
            nrows = min(ROW_BLOCK_SIZE, ysize - y0)
            y1 = y0 + nrows

            md = md_band.ReadAsArray(0, y0, xsize, nrows).astype(np.float32)
            if md.ndim == 1:
                md = md.reshape(1, -1)
            valid = np.isfinite(md)
            if md_nodata is not None and np.isfinite(md_nodata):
                valid &= md != md_nodata
            valid &= md != FLOAT_NODATA
            keep = valid & (md <= threshold)

            mask = np.full((nrows, xsize), BYTE_NODATA, dtype=np.uint8)
            mask[valid] = np.where(keep[valid], 0, 1).astype(np.uint8)
            mask_band.WriteArray(mask, 0, y0)

            pred_block = np.array(pred_mm[y0:y1, :], dtype=np.float32, order="C", copy=True)
            pred_valid = np.isfinite(pred_block) & (pred_block != FLOAT_NODATA)
            pred_out = np.full((nrows, xsize), FLOAT_NODATA, dtype=np.float32)
            pred_keep = keep & pred_valid
            pred_out[pred_keep] = pred_block[pred_keep]
            pred_band.WriteArray(pred_out, 0, y0)

            kept += int(np.sum(keep))
            valid_total += int(np.sum(valid))

            if PRINT_EVERY_BLOCK or y1 % PROGRESS_EVERY == 0 or y1 == ysize:
                elapsed_h = (time.time() - start) / 3600.0
                done = y1 / ysize
                eta_h = elapsed_h / max(done, 1e-9) - elapsed_h
                print(
                    f"{label} {profile} | rows {y1}/{ysize} "
                    f"| elapsed h {elapsed_h:.2f} | ETA h {eta_h:.2f}",
                    flush=True,
                )

        mask_band.FlushCache()
        mask_ds.FlushCache()
        pred_band.FlushCache()
        pred_ds.FlushCache()

        summary = pd.DataFrame([{
            "dataset": dataset,
            "scenario": scenario,
            "profile": profile,
            "threshold_type": threshold_info["threshold_type"],
            "threshold": threshold,
            "valid_pixels": valid_total,
            "kept_pixels": kept,
            "in_distribution_fraction": kept / valid_total if valid_total else np.nan,
        }])
        summary.to_csv(
            os.path.join(ridge_output_dir, f"{label}_ridge_{profile}_pixel_summary.csv"),
            index=False,
        )
        print(
            f"{label} {profile} done | in-distribution fraction = "
            f"{kept / valid_total if valid_total else np.nan:.4f}",
            flush=True,
        )

    finally:
        md_band = None
        md_ds = None
        mask_band = None
        mask_ds = None
        pred_band = None
        pred_ds = None
        if pred_mm is not None:
            del pred_mm
        if pred_mm_path is not None:
            try:
                os.remove(pred_mm_path)
            except OSError:
                pass


def main():
    if DATASET_OVERRIDE is not None and SCENARIO_OVERRIDE is not None:
        run_one(DATASET_OVERRIDE, SCENARIO_OVERRIDE)
        return

    args = sys.argv[1:]
    if len(args) != 2:
        raise SystemExit("Usage: python td_gfm_add_ridge_90_profile_core.py <TD|GFM> <scenario>")
    run_one(args[0], args[1])


if __name__ == "__main__":
    main()
