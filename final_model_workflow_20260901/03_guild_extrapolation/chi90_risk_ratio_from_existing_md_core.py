#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Shared streaming helper for chi90 risk products from an MD GeoTIFF."""

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

FLOAT_NODATA = -9999.0
BYTE_NODATA = 255
ROW_BLOCK = 512


def _create_output(path, reference, data_type, nodata):
    dataset = gdal.GetDriverByName("GTiff").Create(
        path,
        reference.RasterXSize,
        reference.RasterYSize,
        1,
        data_type,
        options=["TILED=YES", "COMPRESS=LZW", "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS"],
    )
    dataset.SetGeoTransform(reference.GetGeoTransform())
    dataset.SetProjection(reference.GetProjection())
    band = dataset.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    return dataset, band


def run(config):
    """Write a chi90 risk ratio, a chi90 mask and a small summary CSV."""
    md_path = config["md_tif"]
    output_dir = config["output_dir"]
    degrees_of_freedom = int(config["df"])
    chi_square_d2 = float(config["chi_square_d2"])
    threshold = float(np.sqrt(chi_square_d2))
    abundance_tag = "_abundance" if config.get("abundance", False) else ""

    if not os.path.exists(md_path):
        raise FileNotFoundError("Missing input Mahalanobis TIFF: {}".format(md_path))
    if os.path.exists(output_dir) and os.listdir(output_dir):
        raise FileExistsError(
            "Output directory is not empty; stopping to protect its contents: {}".format(
                output_dir
            )
        )
    os.makedirs(output_dir, exist_ok=True)

    source = gdal.Open(md_path)
    if source is None:
        raise RuntimeError("Unable to open: {}".format(md_path))
    source_band = source.GetRasterBand(1)
    source_nodata = source_band.GetNoDataValue()
    xsize, ysize = source.RasterXSize, source.RasterYSize

    ratio_path = os.path.join(
        output_dir, "environmental_extrapolation_risk_ratio{}.tif".format(abundance_tag)
    )
    mask_path = os.path.join(
        output_dir, "environmental_extrapolation_mask{}_chi90.tif".format(abundance_tag)
    )
    ratio_dataset, ratio_band = _create_output(
        ratio_path, source, gdal.GDT_Float32, FLOAT_NODATA
    )
    mask_dataset, mask_band = _create_output(
        mask_path, source, gdal.GDT_Byte, BYTE_NODATA
    )
    ratio_dataset.SetMetadataItem(
        "RISK_DEFINITION",
        "ridge_Mahalanobis_distance / sqrt(qchisq(0.90, df={}))".format(degrees_of_freedom),
    )
    ratio_dataset.SetMetadataItem("CHI90_DISTANCE_THRESHOLD", "{:.10f}".format(threshold))

    valid_pixels = 0
    in_distribution_pixels = 0
    minimum_ratio = np.inf
    maximum_ratio = -np.inf
    started = time.time()
    try:
        print("=== {} chi90 from existing Mahalanobis distance ===".format(config["model"]), flush=True)
        print("Input:", md_path, flush=True)
        print("Threshold sqrt(qchisq(0.90, df={})): {:.10f}".format(degrees_of_freedom, threshold), flush=True)
        print("Output directory:", output_dir, flush=True)

        for y0 in range(0, ysize, ROW_BLOCK):
            height = min(ROW_BLOCK, ysize - y0)
            distance = source_band.ReadAsArray(0, y0, xsize, height).astype(np.float32)
            valid = np.isfinite(distance) & (distance != FLOAT_NODATA)
            if source_nodata is not None and np.isfinite(source_nodata):
                valid &= distance != source_nodata

            ratio = np.full((height, xsize), FLOAT_NODATA, dtype=np.float32)
            ratio[valid] = distance[valid] / threshold
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
        with open(
            os.path.join(output_dir, "{}_chi90_threshold_and_summary.csv".format(config["model"])),
            "w",
            newline="",
        ) as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "model", "threshold_type", "probability", "df", "chi_square_threshold_d2",
                    "distance_threshold", "features", "valid_pixels", "in_distribution_pixels",
                    "in_distribution_fraction", "minimum_ratio", "maximum_ratio",
                ],
            )
            writer.writeheader()
            writer.writerow({
                "model": config["model"],
                "threshold_type": "chi_square",
                "probability": 0.90,
                "df": degrees_of_freedom,
                "chi_square_threshold_d2": chi_square_d2,
                "distance_threshold": threshold,
                "features": ";".join(config["features"]),
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
