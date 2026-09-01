#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fast G2 chi90 risk ratio and mask from an existing Mahalanobis GeoTIFF.

This does not refit random forest or read environmental/prediction CSV files.
It assumes mahalanobis_distance.tif stores D = sqrt(D^2) and was calculated
for the G2 predictors ph, ref_1 and ndvimean.
"""

import csv
import os
import time

import numpy as np

try:
    from osgeo import gdal
except Exception:
    import gdal


gdal.UseExceptions()
gdal.SetConfigOption("GDAL_CACHEMAX", "512")

RESULT_DIR = "/public/home/ylwang/Gisdata/TC_New/lucy/uncertainty/g2_results"
MD_TIF = os.path.join(RESULT_DIR, "mahalanobis_distance.tif")
OUTPUT_DIR = os.path.join(RESULT_DIR, "chi90_from_existing_mahalanobis_20260730")

FLOAT_NODATA = -9999.0
BYTE_NODATA = 255
ROW_BLOCK = 512

# G2 has p = 3 final predictors: ph, ref_1 and ndvimean.
CHI_PROBABILITY = 0.90
DF = 3
CHI_SQUARE_D2 = 6.251388631170325  # qchisq(0.90, df = 3)
CHI90_THRESHOLD = float(np.sqrt(CHI_SQUARE_D2))


def create_output(path, reference, dtype, nodata):
    dataset = gdal.GetDriverByName("GTiff").Create(
        path,
        reference.RasterXSize,
        reference.RasterYSize,
        1,
        dtype,
        options=["TILED=YES", "COMPRESS=LZW", "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS"],
    )
    dataset.SetGeoTransform(reference.GetGeoTransform())
    dataset.SetProjection(reference.GetProjection())
    band = dataset.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    return dataset, band


def main():
    if not os.path.exists(MD_TIF):
        raise FileNotFoundError("Missing input Mahalanobis TIFF: {}".format(MD_TIF))
    if os.path.exists(OUTPUT_DIR) and os.listdir(OUTPUT_DIR):
        raise FileExistsError(
            "Output directory is not empty; stopping to protect its contents: {}".format(
                OUTPUT_DIR
            )
        )
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    source = gdal.Open(MD_TIF)
    if source is None:
        raise RuntimeError("Unable to open: {}".format(MD_TIF))
    source_band = source.GetRasterBand(1)
    source_nodata = source_band.GetNoDataValue()
    xsize, ysize = source.RasterXSize, source.RasterYSize

    ratio_path = os.path.join(OUTPUT_DIR, "environmental_extrapolation_risk_ratio.tif")
    mask_path = os.path.join(OUTPUT_DIR, "environmental_extrapolation_mask_chi90.tif")
    ratio_dataset, ratio_band = create_output(ratio_path, source, gdal.GDT_Float32, FLOAT_NODATA)
    mask_dataset, mask_band = create_output(mask_path, source, gdal.GDT_Byte, BYTE_NODATA)
    ratio_dataset.SetMetadataItem(
        "RISK_DEFINITION", "ridge_Mahalanobis_distance / sqrt(qchisq(0.90, df=3))"
    )
    ratio_dataset.SetMetadataItem("CHI90_DISTANCE_THRESHOLD", "{:.10f}".format(CHI90_THRESHOLD))

    valid_pixels = 0
    in_distribution_pixels = 0
    minimum_ratio = np.inf
    maximum_ratio = -np.inf
    started = time.time()
    try:
        print("=== G2 chi90 from existing Mahalanobis distance ===", flush=True)
        print("Input:", MD_TIF, flush=True)
        print("Threshold sqrt(qchisq(0.90, df=3)): {:.10f}".format(CHI90_THRESHOLD), flush=True)
        print("Output directory:", OUTPUT_DIR, flush=True)

        for y0 in range(0, ysize, ROW_BLOCK):
            height = min(ROW_BLOCK, ysize - y0)
            distance = source_band.ReadAsArray(0, y0, xsize, height).astype(np.float32)
            valid = np.isfinite(distance) & (distance != FLOAT_NODATA)
            if source_nodata is not None and np.isfinite(source_nodata):
                valid &= distance != source_nodata

            ratio = np.full((height, xsize), FLOAT_NODATA, dtype=np.float32)
            ratio[valid] = distance[valid] / CHI90_THRESHOLD
            mask = np.full((height, xsize), BYTE_NODATA, dtype=np.uint8)
            mask[valid] = np.where(ratio[valid] <= 1.0, 0, 1).astype(np.uint8)

            ratio_band.WriteArray(ratio, 0, y0)
            mask_band.WriteArray(mask, 0, y0)
            valid_pixels += int(valid.sum())
            in_distribution_pixels += int(np.sum(valid & (ratio <= 1.0)))
            if np.any(valid):
                minimum_ratio = min(minimum_ratio, float(np.min(ratio[valid])))
                maximum_ratio = max(maximum_ratio, float(np.max(ratio[valid])))

            elapsed = (time.time() - started) / 3600.0
            done = (y0 + height) / float(ysize)
            eta = elapsed / max(done, 1e-12) - elapsed
            print(
                "  rows {}/{} | elapsed h {:.2f} | ETA h {:.2f}".format(
                    y0 + height, ysize, elapsed, eta
                ),
                flush=True,
            )

        ratio_band.FlushCache()
        ratio_dataset.FlushCache()
        mask_band.FlushCache()
        mask_dataset.FlushCache()

        with open(os.path.join(OUTPUT_DIR, "g2_chi90_threshold_and_summary.csv"), "w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "model", "threshold_type", "probability", "df", "chi_square_threshold_d2",
                    "distance_threshold", "risk_definition", "valid_pixels", "in_distribution_pixels",
                    "in_distribution_fraction", "minimum_ratio", "maximum_ratio",
                ],
            )
            writer.writeheader()
            writer.writerow({
                "model": "g2",
                "threshold_type": "chi_square",
                "probability": CHI_PROBABILITY,
                "df": DF,
                "chi_square_threshold_d2": CHI_SQUARE_D2,
                "distance_threshold": CHI90_THRESHOLD,
                "risk_definition": "ridge_Mahalanobis_distance_over_sqrt_qchisq_0.90_df3",
                "valid_pixels": valid_pixels,
                "in_distribution_pixels": in_distribution_pixels,
                "in_distribution_fraction": in_distribution_pixels / valid_pixels if valid_pixels else np.nan,
                "minimum_ratio": minimum_ratio if valid_pixels else np.nan,
                "maximum_ratio": maximum_ratio if valid_pixels else np.nan,
            })
        print("Done.", flush=True)
        print("Risk ratio:", ratio_path, flush=True)
        print("chi90 mask:", mask_path, flush=True)
    finally:
        ratio_band = None
        ratio_dataset = None
        mask_band = None
        mask_dataset = None
        source_band = None
        source = None


if __name__ == "__main__":
    main()
