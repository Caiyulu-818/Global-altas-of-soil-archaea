#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Create the missing G1-abundance chi90 environmental-risk-ratio GeoTIFF.

This is a post-processing step only. It does not retrain the random forest,
does not read the prediction CSV and does not modify any existing prediction
or mask. It divides an existing ridge-Mahalanobis distance GeoTIFF by the
theoretical chi-square 0.90 threshold for three predictors.

The input distance layer must be based on the same G1 abundance model:
ref_1, ndvimean and Elevation. Distances are assumed to be square-root
Mahalanobis distances, D = sqrt(D^2), as produced by the current workflow.
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

BASE_DIR = "/public/home/ylwang/Gisdata/TC_New/lucy"
RESULT_DIR = os.path.join(BASE_DIR, "uncertainty", "g1_results")

# The first existing candidate is used. Set G1_ABUNDANCE_MD_TIF before running
# to explicitly use a different valid Mahalanobis-distance GeoTIFF.
MD_CANDIDATES = [
    os.environ.get("G1_ABUNDANCE_MD_TIF", ""),
    os.path.join(
        RESULT_DIR,
        "ridge_mahalanobis_ood_tif_direct0728",
        "mahalanobis_distance_abundance.tif",
    ),
    os.path.join(
        RESULT_DIR,
        "ridge_mahalanobis_ood_tif_direct1",
        "mahalanobis_distance_abundance.tif",
    ),
]

# The missing ratio is added beside the matching distance layer. No previous
# output is deleted or overwritten unless OVERWRITE is explicitly set to True.
OUTPUT_NAME = "environmental_extrapolation_risk_ratio_abundance.tif"
SUMMARY_NAME = "g1_abundance_chi90_environmental_risk_ratio_summary.csv"
OVERWRITE = False

FLOAT_NODATA = -9999.0
ROW_BLOCK = 256
CHI_SQUARE_0P90_DF3 = 6.251388631170325  # qchisq(0.90, df = 3)
CHI90_THRESHOLD = float(np.sqrt(CHI_SQUARE_0P90_DF3))


def find_md_tif():
    for candidate in MD_CANDIDATES:
        if candidate and os.path.exists(candidate):
            return candidate
    raise FileNotFoundError(
        "Could not find mahalanobis_distance_abundance.tif. Checked:\n{}\n"
        "Set G1_ABUNDANCE_MD_TIF to the correct distance-layer path, then rerun.".format(
            "\n".join(path for path in MD_CANDIDATES if path)
        )
    )


def create_output(path, reference):
    if os.path.exists(path):
        if not OVERWRITE:
            raise FileExistsError(
                "Output already exists and will not be overwritten: {}".format(path)
            )
        os.remove(path)
    if os.path.exists(path + ".aux.xml"):
        os.remove(path + ".aux.xml")

    dataset = gdal.GetDriverByName("GTiff").Create(
        path,
        reference.RasterXSize,
        reference.RasterYSize,
        1,
        gdal.GDT_Float32,
        options=["TILED=YES", "COMPRESS=LZW", "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS"],
    )
    dataset.SetGeoTransform(reference.GetGeoTransform())
    dataset.SetProjection(reference.GetProjection())
    dataset.SetMetadataItem(
        "RISK_DEFINITION",
        "ridge_Mahalanobis_distance / sqrt(qchisq(0.90, df=3))",
    )
    dataset.SetMetadataItem("CHI90_DISTANCE_THRESHOLD", "{:.10f}".format(CHI90_THRESHOLD))
    band = dataset.GetRasterBand(1)
    band.SetNoDataValue(FLOAT_NODATA)
    return dataset, band


def main():
    md_path = find_md_tif()
    output_dir = os.path.dirname(md_path)
    output_path = os.path.join(output_dir, OUTPUT_NAME)
    summary_path = os.path.join(output_dir, SUMMARY_NAME)

    source = gdal.Open(md_path)
    if source is None:
        raise RuntimeError("Unable to open: {}".format(md_path))
    source_band = source.GetRasterBand(1)
    source_nodata = source_band.GetNoDataValue()
    xsize, ysize = source.RasterXSize, source.RasterYSize
    output, output_band = create_output(output_path, source)

    valid_pixels = 0
    in_distribution_pixels = 0
    max_ratio = -np.inf
    min_ratio = np.inf
    started = time.time()
    try:
        print("=== G1 abundance chi90 environmental extrapolation risk ratio ===", flush=True)
        print("Input Mahalanobis distance:", md_path, flush=True)
        print("chi90 distance threshold: {:.10f}".format(CHI90_THRESHOLD), flush=True)
        print("Output risk ratio:", output_path, flush=True)

        for y0 in range(0, ysize, ROW_BLOCK):
            height = min(ROW_BLOCK, ysize - y0)
            distance = source_band.ReadAsArray(0, y0, xsize, height).astype(np.float32)
            valid = np.isfinite(distance) & (distance != FLOAT_NODATA)
            if source_nodata is not None and np.isfinite(source_nodata):
                valid &= distance != source_nodata

            ratio = np.full((height, xsize), FLOAT_NODATA, dtype=np.float32)
            ratio[valid] = distance[valid] / CHI90_THRESHOLD
            output_band.WriteArray(ratio, 0, y0)

            valid_pixels += int(valid.sum())
            in_distribution_pixels += int(np.sum(valid & (ratio <= 1.0)))
            if np.any(valid):
                min_ratio = min(min_ratio, float(np.min(ratio[valid])))
                max_ratio = max(max_ratio, float(np.max(ratio[valid])))

            done = (y0 + height) / float(ysize)
            elapsed = (time.time() - started) / 3600.0
            eta = elapsed / max(done, 1e-12) - elapsed
            print(
                "  rows {}/{} | elapsed h {:.2f} | ETA h {:.2f}".format(
                    y0 + height, ysize, elapsed, eta
                ),
                flush=True,
            )

        output_band.FlushCache()
        output.FlushCache()
        with open(summary_path, "w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "model", "risk_definition", "chi_square_probability", "df",
                    "chi90_distance_threshold", "valid_pixels", "in_distribution_pixels",
                    "in_distribution_fraction", "minimum_ratio", "maximum_ratio",
                ],
            )
            writer.writeheader()
            writer.writerow({
                "model": "g1_abundance",
                "risk_definition": "ridge_Mahalanobis_distance_over_chi90_threshold",
                "chi_square_probability": 0.90,
                "df": 3,
                "chi90_distance_threshold": CHI90_THRESHOLD,
                "valid_pixels": valid_pixels,
                "in_distribution_pixels": in_distribution_pixels,
                "in_distribution_fraction": (
                    in_distribution_pixels / valid_pixels if valid_pixels else np.nan
                ),
                "minimum_ratio": min_ratio if valid_pixels else np.nan,
                "maximum_ratio": max_ratio if valid_pixels else np.nan,
            })

        print("Done.", flush=True)
        print("Summary:", summary_path, flush=True)
    finally:
        output_band = None
        output = None
        source_band = None
        source = None


if __name__ == "__main__":
    main()
