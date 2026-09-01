#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compile reviewer-facing supplementary tables for TD and GFM.

The source folders contain one set of ridge-Mahalanobis summary files for each
dataset and cropland scenario.  This script keeps model-validation information
separate from projection/OOD information so that sample sizes, predictors and
threshold definitions cannot be conflated.
"""

from __future__ import annotations

import csv
import glob
import math
import os
import re
from collections import defaultdict


BASE = "/Users/caiyulu/Documents/Codex/2026-06-15/files-mentioned-by-the-user-env"
OUT_DIR = os.path.join(BASE, "supplementary_extrapolation_tables_td_gfm_20260801")
SUMMARY_DIRS = {
    "TD": "/Volumes/DEEP/extra/TD/summary",
    "GFM": "/Volumes/DEEP/extra/GFM/summary",
}
VALIDATION_FILE = os.path.join(
    BASE, "supplementary_tables", "supplementary_model_validation_summary_guilds_TD_GFM.csv"
)


def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows, fieldnames):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def value(row, key, default=""):
    return row.get(key, default) if row else default


def clean(value_text):
    if value_text is None:
        return ""
    text = str(value_text).strip()
    if text.lower() in {"nan", "na", "null", "none"}:
        return ""
    return text


def as_float(text):
    text = clean(text)
    if text == "":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def fmt_number(value_number, digits=6):
    if value_number is None:
        return ""
    return f"{value_number:.{digits}f}".rstrip("0").rstrip(".")


def scenario_parts(scenario):
    year_match = re.search(r"(2020|2100)", scenario)
    ssp_match = re.search(r"SSP(126|245|370)", scenario)
    return (
        year_match.group(1) if year_match else "",
        f"SSP{ssp_match.group(1)}" if ssp_match else "",
    )


def scenario_sort_key(row):
    dataset_order = {"TD": 0, "GFM": 1}
    scenario = row.get("scenario", "")
    year, ssp = scenario_parts(scenario)
    return (
        dataset_order.get(row.get("dataset", ""), 99),
        0 if year == "2020" else 1,
        ssp,
    )


def parse_source_files():
    records = defaultdict(dict)
    for dataset, root in SUMMARY_DIRS.items():
        paths = sorted(glob.glob(os.path.join(root, f"{dataset}_*_ridge_*.csv")))
        for path in paths:
            basename = os.path.basename(path)
            match = re.match(rf"{dataset}_(.+?)_ridge_(.+)\.csv$", basename)
            if not match:
                continue
            scenario, suffix = match.groups()
            key = (dataset, scenario)
            records[key][suffix] = (path, read_csv(path))
    return records


def select_one(records_for_scenario, suffix):
    item = records_for_scenario.get(suffix)
    if not item:
        return {}
    rows = item[1]
    return rows[0] if rows else {}


def split_features(text):
    return [clean(x) for x in clean(text).split(";") if clean(x)]


def model_metadata(validation_rows):
    return {row.get("dataset", "").strip(): row for row in validation_rows}


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    validation_rows = read_csv(VALIDATION_FILE)
    validation = model_metadata(
        [row for row in validation_rows if row.get("dataset", "").strip() in {"TD", "GFM"}]
    )
    source = parse_source_files()

    if not source:
        raise RuntimeError("No TD/GFM summary CSV files were found.")
    if set(validation) != {"TD", "GFM"}:
        raise RuntimeError("The validation summary must contain one TD row and one GFM row.")

    model_rows = []
    projection_rows = []
    sensitivity_rows = []
    qc_rows = []

    for (dataset, scenario), files in sorted(source.items(), key=lambda item: scenario_sort_key({"dataset": item[0][0], "scenario": item[0][1]})):
        val = validation[dataset]
        chi_threshold = select_one(files, "mahalanobis_threshold_chi90.csv")
        chi_pixel = select_one(files, "chi90_pixel_summary.csv")
        p90_threshold = select_one(files, "mahalanobis_threshold_p90.csv")
        p90_pixel = select_one(files, "p90_pixel_summary.csv")
        multilevel = files.get("mahalanobis_thresholds.csv", ("", []))[1]

        projection_features = clean(value(chi_threshold, "features"))
        validation_features = clean(value(val, "environmental_factors"))
        projection_feature_list = split_features(projection_features)
        validation_feature_list = split_features(validation_features)
        extra_projection_features = [x for x in projection_feature_list if x not in validation_feature_list]
        missing_projection_features = [x for x in validation_feature_list if x not in projection_feature_list]
        feature_set_match = not extra_projection_features and not missing_projection_features

        year, ssp = scenario_parts(scenario)
        valid_pixels = as_float(value(chi_pixel, "valid_pixels"))
        kept_pixels = as_float(value(chi_pixel, "kept_pixels"))
        in_fraction = as_float(value(chi_pixel, "in_distribution_fraction"))
        extrapolated_fraction = 1.0 - in_fraction if in_fraction is not None else None

        row = {
            "dataset": dataset,
            "scenario": scenario,
            "scenario_year": year,
            "scenario_ssp": ssp,
            "response_type": clean(value(val, "response_type")),
            "model_type": clean(value(val, "model_type")),
            "validation_feature_set": validation_features,
            "validation_spatial_terms": clean(value(val, "spatial_terms")),
            "projection_feature_set": projection_features,
            "projection_feature_count": len(projection_feature_list),
            "validation_n": clean(value(val, "n")),
            "ood_threshold_training_n": clean(value(chi_threshold, "n_training_rows")),
            "primary_validation_method": clean(value(val, "best_validation_method")),
            "primary_R2": clean(value(val, "best_R2")),
            "primary_RMSE": clean(value(val, "best_RMSE")),
            "spatial_validation_method": clean(value(val, "spatial_cv_method")),
            "spatial_R2": clean(value(val, "spatial_cv_R2")),
            "spatial_RMSE": "",
            "moran_observed_I": clean(value(val, "moran_observed_I")),
            "moran_expected_I": clean(value(val, "moran_expected_I")),
            "moran_p_value": clean(value(val, "moran_p_value")),
            "moran_pass": clean(value(val, "moran_pass")),
            "threshold_profile": clean(value(chi_pixel, "profile")) or "chi90",
            "threshold_type": clean(value(chi_threshold, "threshold_type")),
            "threshold_probability": clean(value(chi_threshold, "probability")),
            "chi_square_threshold_D2": clean(value(chi_threshold, "chi_square_threshold_d2")),
            "distance_threshold": clean(value(chi_threshold, "threshold")),
            "ridge": clean(value(chi_threshold, "ridge")),
            "distance_space": clean(value(chi_threshold, "space")),
            "valid_pixels": clean(value(chi_pixel, "valid_pixels")),
            "kept_pixels": clean(value(chi_pixel, "kept_pixels")),
            "in_distribution_fraction": clean(value(chi_pixel, "in_distribution_fraction")),
            "extrapolated_fraction": fmt_number(extrapolated_fraction, 9),
            "projection_fraction_denominator": "complete environmental-predictor pixels",
            "feature_set_consistent": str(feature_set_match),
            "extra_projection_features": "; ".join(extra_projection_features),
            "missing_projection_features": "; ".join(missing_projection_features),
            "geographic_distance_in_summary": "not reported",
            "note": "Scenario projection summary uses ridge-Mahalanobis environmental distance; reconcile cropland_scenario with model-fitting description.",
        }
        projection_rows.append(row)

        # One model-level row, duplicated only once per dataset below.
        model_rows.append({
            "dataset": dataset,
            "target": clean(value(val, "target")),
            "response_type": clean(value(val, "response_type")),
            "model_type": clean(value(val, "model_type")),
            "validation_feature_set": validation_features,
            "validation_spatial_terms": clean(value(val, "spatial_terms")),
            "validation_n": clean(value(val, "n")),
            "primary_validation_method": clean(value(val, "best_validation_method")),
            "primary_R2": clean(value(val, "best_R2")),
            "primary_RMSE": clean(value(val, "best_RMSE")),
            "spatial_validation_method": clean(value(val, "spatial_cv_method")),
            "spatial_R2": clean(value(val, "spatial_cv_R2")),
            "moran_observed_I": clean(value(val, "moran_observed_I")),
            "moran_expected_I": clean(value(val, "moran_expected_I")),
            "moran_p_value": clean(value(val, "moran_p_value")),
            "moran_pass": clean(value(val, "moran_pass")),
            "cv_pass": clean(value(val, "cv_pass")),
            "recommendation_or_note": clean(value(val, "recommendation_or_note")),
        })

        # Both chi90 and empirical p90 are retained for threshold sensitivity.
        for label, threshold_row, pixel_row in (
            ("chi90", chi_threshold, chi_pixel),
            ("empirical_p90", p90_threshold, p90_pixel),
        ):
            threshold_value = as_float(value(threshold_row, "threshold"))
            fraction = as_float(value(pixel_row, "in_distribution_fraction"))
            sensitivity_rows.append({
                "dataset": dataset,
                "scenario": scenario,
                "profile": label,
                "threshold_type": clean(value(threshold_row, "threshold_type")),
                "percentile": clean(value(threshold_row, "percentile")),
                "probability": clean(value(threshold_row, "probability")),
                "df": clean(value(threshold_row, "df")),
                "chi_square_threshold_D2": clean(value(threshold_row, "chi_square_threshold_d2")),
                "distance_threshold": clean(value(threshold_row, "threshold")),
                "ridge": clean(value(threshold_row, "ridge")),
                "distance_space": clean(value(threshold_row, "space")),
                "valid_pixels": clean(value(pixel_row, "valid_pixels")),
                "kept_pixels": clean(value(pixel_row, "kept_pixels")),
                "in_distribution_fraction": clean(value(pixel_row, "in_distribution_fraction")),
                "extrapolated_fraction": fmt_number(1.0 - fraction, 9) if fraction is not None else "",
                "interpretation": "heuristic threshold sensitivity; not a formal significance test",
            })

        # QC catches the points a reviewer is most likely to ask about.
        multilevel_profiles = "; ".join(
            f"{clean(value(r, 'profile'))}:{clean(value(r, 'threshold'))}"
            for r in multilevel if clean(value(r, "profile"))
        )
        qc_rows.append({
            "dataset": dataset,
            "scenario": scenario,
            "chi90_summary_present": str(bool(chi_pixel)),
            "p90_summary_present": str(bool(p90_pixel)),
            "chi90_threshold_present": str(bool(chi_threshold)),
            "p90_threshold_present": str(bool(p90_threshold)),
            "valid_pixels": clean(value(chi_pixel, "valid_pixels")),
            "chi90_in_distribution_fraction": clean(value(chi_pixel, "in_distribution_fraction")),
            "p90_in_distribution_fraction": clean(value(p90_pixel, "in_distribution_fraction")),
            "threshold_training_n": clean(value(chi_threshold, "n_training_rows")),
            "validation_n": clean(value(val, "n")),
            "threshold_n_equals_validation_n": str(clean(value(chi_threshold, "n_training_rows")) == clean(value(val, "n"))),
            "validation_features": validation_features,
            "projection_features": projection_features,
            "feature_set_consistent": str(feature_set_match),
            "extra_projection_features": "; ".join(extra_projection_features),
            "missing_projection_features": "; ".join(missing_projection_features),
            "multilevel_thresholds_available": multilevel_profiles,
            "reviewer_action": (
                "Clarify whether cropland_scenario is an additional scenario-only predictor and why OOD n differs from validation n."
                if (not feature_set_match or clean(value(chi_threshold, "n_training_rows")) != clean(value(val, "n")))
                else "No structural discrepancy detected in the supplied summaries."
            ),
        })

    # Deduplicate the model-level table by dataset.
    unique_models = {}
    for row in model_rows:
        unique_models[row["dataset"]] = row
    model_rows = [unique_models[key] for key in ("TD", "GFM") if key in unique_models]

    # A compact joined table is convenient for a supplementary spreadsheet.
    joined_rows = []
    for row in projection_rows:
        joined_rows.append(row)

    model_fields = list(model_rows[0].keys())
    projection_fields = list(projection_rows[0].keys())
    sensitivity_fields = list(sensitivity_rows[0].keys())
    qc_fields = list(qc_rows[0].keys())

    write_csv(os.path.join(OUT_DIR, "Supplementary_Table_TD_GFM_Model_Validation.csv"), model_rows, model_fields)
    write_csv(os.path.join(OUT_DIR, "Supplementary_Table_TD_GFM_Projection_chi90.csv"), projection_rows, projection_fields)
    write_csv(os.path.join(OUT_DIR, "Supplementary_Table_TD_GFM_Threshold_Sensitivity.csv"), sensitivity_rows, sensitivity_fields)
    write_csv(os.path.join(OUT_DIR, "Supplementary_Table_TD_GFM_Validation_and_Projection.csv"), joined_rows, projection_fields)
    write_csv(os.path.join(OUT_DIR, "Supplementary_Table_TD_GFM_QC.csv"), qc_rows, qc_fields)

    readme = f"""TD/GFM supplementary tables
============================

Source directories:
  TD:  {SUMMARY_DIRS['TD']}
  GFM: {SUMMARY_DIRS['GFM']}

Primary projection table:
  Supplementary_Table_TD_GFM_Projection_chi90.csv
  Eight rows, one per dataset x cropland scenario. The primary threshold is the theoretical chi-square 0.90 threshold.

Model validation table:
  Supplementary_Table_TD_GFM_Model_Validation.csv
  Two rows, one per endpoint. Validation n and OOD-threshold n are intentionally separate.

Threshold sensitivity:
  Supplementary_Table_TD_GFM_Threshold_Sensitivity.csv
  Compares theoretical chi-square 0.90 with the empirical training-distance 90th percentile.

Reviewer-facing QC:
  Supplementary_Table_TD_GFM_QC.csv
  Explicitly flags differences between validation and projection predictor sets and sample sizes.

Interpretation:
  in_distribution_fraction = kept_pixels / valid_pixels.
  extrapolated_fraction = 1 - in_distribution_fraction.
  valid_pixels are cells with complete projection predictors; ocean/no-data cells are excluded from the denominator.
  The supplied projection summaries include cropland_scenario as an additional predictor and report n_training_rows = 781 (TD) or 787 (GFM), whereas the validation summary reports n = 251. This is not silently harmonized: the discrepancy is preserved and flagged for Methods clarification.
"""
    with open(os.path.join(OUT_DIR, "README.txt"), "w", encoding="utf-8") as handle:
        handle.write(readme)

    print(f"Wrote {len(model_rows)} model-validation rows")
    print(f"Wrote {len(projection_rows)} chi90 projection rows")
    print(f"Wrote {len(sensitivity_rows)} threshold-sensitivity rows")
    print(f"Wrote {len(qc_rows)} QC rows")
    print(f"Output directory: {OUT_DIR}")


if __name__ == "__main__":
    main()
