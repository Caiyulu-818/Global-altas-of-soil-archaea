#!/usr/bin/env python3
"""Compile reviewer-oriented guild extrapolation supplementary tables."""

from __future__ import annotations

import csv
import math
from pathlib import Path


SUMMARY_DIR = Path("/Volumes/DEEP/extra/Guild/summary")
VALIDATION_FILE = Path(
    "/Users/caiyulu/Documents/Codex/2026-06-15/files-mentioned-by-the-user-env/"
    "supplementary_tables/supplementary_model_validation_summary_guilds_TD_GFM.csv"
)
OUT_DIR = Path(__file__).resolve().parent / "supplementary_extrapolation_tables_20260801"


ENDPOINTS = {
    "g1_occurrence": {
        "guild": "G1",
        "response_type": "occurrence suitability",
        "model_label": "G1 occurrence",
        "primary_file": "g1_occurrence_chi90_threshold_and_summary.csv",
    },
    "g1_abundance": {
        "guild": "G1",
        "response_type": "abundance",
        "model_label": "G1 abundance",
        "primary_file": "g1_abundance_chi90_risk_ratio_pixel_summary.csv",
        "threshold_file": "g1_abundance_ridge_mahalanobis_threshold_chi90.csv",
    },
    "g2_abundance": {
        "guild": "G2",
        "response_type": "abundance",
        "model_label": "G2 abundance",
        "primary_file": "g2_ridge_chi90_pixel_summary.csv",
        "threshold_file": "g2_ridge_mahalanobis_threshold_chi90.csv",
    },
    "g3_occurrence": {
        "guild": "G3",
        "response_type": "occurrence suitability",
        "model_label": "G3 occurrence",
        "primary_file": "g3_occurrence_chi90_threshold_and_summary.csv",
    },
    "g3_abundance": {
        "guild": "G3",
        "response_type": "abundance",
        "model_label": "G3 abundance",
        "primary_file": "g3_abundance_chi90_threshold_and_summary.csv",
    },
    "g4_abundance": {
        "guild": "G4",
        "response_type": "abundance",
        "model_label": "G4 abundance",
        "primary_file": "g4_ridge_chi90_pixel_summary.csv",
        "threshold_file": "g4_ridge_mahalanobis_threshold_chi90.csv",
    },
    "g5_abundance": {
        "guild": "G5",
        "response_type": "abundance",
        "model_label": "G5 abundance",
        "primary_file": "g5_ridge_chi90_pixel_summary.csv",
        "threshold_file": "g5_ridge_mahalanobis_threshold_chi90.csv",
    },
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def as_float(value: str | None):
    if value in (None, "", "NA", "nan", "NaN"):
        return None
    try:
        number = float(value)
    except ValueError:
        return value
    if not math.isfinite(number):
        return None
    return number


def first_row(path: Path) -> dict[str, str]:
    rows = read_csv(path)
    if not rows:
        raise ValueError(f"Empty CSV: {path}")
    return rows[0]


def write_csv(path: Path, rows: list[dict], columns: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def normalized_set(value: str | None) -> set[str]:
    if not value:
        return set()
    return {item.strip().lower() for item in value.split(";") if item.strip()}


def load_validation_rows() -> dict[str, dict[str, str]]:
    if not VALIDATION_FILE.exists():
        return {}

    endpoint_by_model = {
        "Guild 1": "g1_occurrence",
        "Guild 2": "g2_abundance",
        "Guild 3": "g3_occurrence",
        "Guild 4": "g4_abundance",
        "Guild 5": "g5_abundance",
    }
    result = {}
    for row in read_csv(VALIDATION_FILE):
        endpoint = endpoint_by_model.get(row.get("model_name", ""))
        if endpoint:
            result[endpoint] = row
    return result


def merge_validation_and_ood(main_rows: list[dict]) -> list[dict]:
    validation = load_validation_rows()
    merged = []

    for row in main_rows:
        endpoint = row["endpoint"]
        val = validation.get(endpoint, {})
        validated_predictors = val.get("environmental_factors", "")
        ood_predictors = row.get("predictors", "")
        val_set = normalized_set(validated_predictors)
        ood_set = normalized_set(ood_predictors)

        if not val:
            predictor_status = "MISSING_VALIDATION_ENDPOINT"
            recommendation = "Add endpoint-specific out-of-fold validation before presenting this projection as a final model."
        elif val_set == ood_set:
            predictor_status = "MATCH"
            recommendation = "No environmental-predictor discrepancy detected."
        else:
            predictor_status = "MISMATCH"
            recommendation = "Use the identical final predictor set for validation, global prediction and Mahalanobis extrapolation, or explicitly justify the separation."

        spatial_terms = val.get("spatial_terms", "")
        spatial_status = "NONE_REPORTED" if not spatial_terms else "VALIDATION_ONLY_NOT_IN_OOD_SUMMARY"

        merged.append(
            {
                **row,
                "validated_model_type": val.get("model_type", ""),
                "validated_feature_mode": val.get("feature_mode", ""),
                "validated_environmental_factors": validated_predictors,
                "validated_spatial_terms": spatial_terms,
                "validation_n": val.get("n", ""),
                "best_validation_method": val.get("best_validation_method", ""),
                "best_R2": val.get("best_R2", ""),
                "best_RMSE": val.get("best_RMSE", ""),
                "best_AUC": val.get("best_AUC", ""),
                "best_TSS": val.get("best_TSS", ""),
                "spatial_cv_method": val.get("spatial_cv_method", ""),
                "spatial_cv_R2": val.get("spatial_cv_R2", ""),
                "moran_observed_I": val.get("moran_observed_I", ""),
                "moran_expected_I": val.get("moran_expected_I", ""),
                "moran_p_value": val.get("moran_p_value", ""),
                "predictor_consistency": predictor_status,
                "spatial_terms_consistency": spatial_status,
                "reviewer_action": recommendation,
            }
        )

    return merged


def threshold_from_pixel_file(endpoint: str, cfg: dict, pixel: dict) -> dict:
    threshold_file = cfg.get("threshold_file")
    if threshold_file:
        threshold = first_row(SUMMARY_DIR / threshold_file)
    else:
        threshold = pixel

    model = endpoint
    if model == "g2_abundance":
        model = "g2"
    elif model == "g4_abundance":
        model = "g4"
    elif model == "g5_abundance":
        model = "g5"

    return {
        "endpoint": endpoint,
        "guild": cfg["guild"],
        "response_type": cfg["response_type"],
        "model_label": cfg["model_label"],
        "model_type": "ranger random forest",
        "predictors": threshold.get("features", pixel.get("features", "")),
        "n_training_rows": threshold.get("n_training_rows", ""),
        "distance_method": "ridge-regularized Mahalanobis distance",
        "distance_space": threshold.get(
            "space", "standardized raw predictor space; square-root distance"
        ),
        "ridge_lambda": threshold.get("ridge", "0.001"),
        "threshold_method": threshold.get("threshold_type", "chi_square"),
        "threshold_profile": threshold.get("profile", pixel.get("profile", "chi90")),
        "threshold_probability": threshold.get(
            "probability", pixel.get("probability", "0.9")
        ),
        "degrees_of_freedom": threshold.get("df", pixel.get("df", "")),
        "chi_square_threshold_D2": threshold.get(
            "chi_square_threshold_d2", pixel.get("chi_square_threshold_d2", "")
        ),
        "distance_threshold": threshold.get(
            "threshold", threshold.get("distance_threshold", pixel.get("threshold", ""))
        ),
        "valid_pixels": pixel.get("valid_pixels", ""),
        "in_distribution_pixels": pixel.get(
            "in_distribution_pixels", pixel.get("kept_pixels", "")
        ),
        "in_distribution_fraction": pixel.get("in_distribution_fraction", ""),
        "extrapolated_fraction": (
            None
            if pixel.get("in_distribution_fraction", "") == ""
            else 1.0 - float(pixel["in_distribution_fraction"])
        ),
        "risk_ratio_definition": (
            "grid-cell ridge-Mahalanobis distance / distance threshold"
        ),
        "fraction_denominator": "valid pixels with complete model predictors",
        "land_mask_status": "not recorded in summary; report fixed land-mask definition",
        "source_file": cfg["primary_file"],
        "risk_ratio_min": threshold.get("minimum_ratio", ""),
        "risk_ratio_max": threshold.get("maximum_ratio", ""),
    }


def load_primary(endpoint: str, cfg: dict) -> dict:
    path = SUMMARY_DIR / cfg["primary_file"]
    row = first_row(path)

    # A combined threshold_and_summary file already contains all fields.
    if "distance_threshold" in row:
        return threshold_from_pixel_file(endpoint, cfg, row)

    threshold = threshold_from_pixel_file(endpoint, cfg, row)

    # Add continuous-risk extrema where available.
    if endpoint == "g1_abundance":
        risk = first_row(
            SUMMARY_DIR / "g1_abundance_chi90_risk_ratio_pixel_summary.csv"
        )
        threshold["risk_ratio_min"] = risk.get("minimum_ratio", "")
        threshold["risk_ratio_max"] = risk.get("maximum_ratio", "")

    return threshold


def load_sensitivity_rows() -> list[dict]:
    rows = []
    files = sorted(SUMMARY_DIR.glob("*.csv"))
    for path in files:
        name = path.name
        if name.startswith("._"):
            continue

        # Current, well-defined pixel summaries.
        if "pixel_summary" not in name and "threshold_and_summary" not in name:
            continue
        if "ood_pixel_summary" in name:
            continue

        data = first_row(path)
        model = data.get("model", "")
        if model == "g2":
            endpoint = "g2_abundance"
        elif model == "g4":
            endpoint = "g4_abundance"
        elif model == "g5":
            endpoint = "g5_abundance"
        elif model == "g3":
            endpoint = "g3_occurrence"
        else:
            endpoint = model

        if endpoint not in ENDPOINTS:
            continue

        cfg = ENDPOINTS[endpoint]
        profile = data.get("profile", "")
        # The standalone risk-ratio summary duplicates the chi90 pixel summary.
        if "risk_ratio" in name:
            continue
        if profile == "p90":
            profile = "empirical_p90"
        if profile == "" and data.get("threshold_type") == "chi_square":
            profile = "chi90"
        if profile == "":
            profile = "unknown"

        threshold_lookup = {
            "g2_abundance": "g2_ridge_mahalanobis_threshold_p90.csv",
            "g4_abundance": "g4_ridge_mahalanobis_threshold_p90.csv",
            "g5_abundance": "g5_ridge_mahalanobis_threshold_p90.csv",
        }
        threshold_row = {}
        threshold_name = threshold_lookup.get(endpoint)
        if threshold_name and profile == "empirical_p90":
            threshold_row = first_row(SUMMARY_DIR / threshold_name)

        probability = data.get("probability", "")
        if profile == "empirical_p90" and threshold_row.get("percentile"):
            probability = str(float(threshold_row["percentile"]) / 100.0)

        row = {
            "endpoint": endpoint,
            "guild": cfg["guild"],
            "response_type": cfg["response_type"],
            "profile": profile,
            "threshold_type": threshold_row.get(
                "threshold_type", data.get("threshold_type", "chi_square")
            ),
            "probability": probability,
            "distance_threshold": data.get(
                "threshold", data.get("distance_threshold", "")
            ),
            "valid_pixels": data.get("valid_pixels", ""),
            "in_distribution_pixels": data.get(
                "kept_pixels", data.get("in_distribution_pixels", "")
            ),
            "in_distribution_fraction": data.get("in_distribution_fraction", ""),
            "extrapolated_fraction": (
                None
                if data.get("in_distribution_fraction", "") == ""
                else 1.0 - float(data["in_distribution_fraction"])
            ),
            "source_file": name,
        }
        rows.append(row)

    # Keep one row per endpoint and threshold definition. Prefer the
    # endpoint-specific threshold_and_summary file when both versions exist.
    unique = {}
    for row in rows:
        key = (row["endpoint"], row["profile"], row["threshold_type"])
        if key not in unique or "threshold_and_summary" in row["source_file"]:
            unique[key] = row
    return list(unique.values())


def build_qc_rows(main_rows: list[dict], sensitivity_rows: list[dict]) -> list[dict]:
    present = {row["endpoint"] for row in main_rows}
    qc = []

    for endpoint, cfg in ENDPOINTS.items():
        if endpoint not in present:
            qc.append(
                {
                    "severity": "critical",
                    "endpoint": endpoint,
                    "issue": "missing primary chi90 extrapolation summary",
                    "recommendation": "Generate the endpoint-specific chi90 pixel summary before submission.",
                }
            )

    for row in main_rows:
        if not row.get("n_training_rows"):
            qc.append(
                {
                    "severity": "moderate",
                    "endpoint": row["endpoint"],
                    "issue": "training sample count is absent from the source summary",
                    "recommendation": "Add the number of aggregated training observations.",
                }
            )
        if row.get("land_mask_status", "").startswith("not recorded"):
            qc.append(
                {
                    "severity": "moderate",
                    "endpoint": row["endpoint"],
                    "issue": "land-mask definition is not recorded in the summary",
                    "recommendation": "State the fixed land mask and whether valid_pixels excludes ocean and incomplete predictors.",
                }
            )

    profiles_by_endpoint: dict[str, set[str]] = {}
    for row in sensitivity_rows:
        profiles_by_endpoint.setdefault(row["endpoint"], set()).add(row["profile"])

    for endpoint in ENDPOINTS:
        observed = profiles_by_endpoint.get(endpoint, set())
        if "empirical_p90" not in observed:
            qc.append(
                {
                    "severity": "informational",
                    "endpoint": endpoint,
                    "issue": "empirical p90 pixel coverage is not available",
                    "recommendation": "Include it only if threshold sensitivity is reported; otherwise state that chi90 is the primary threshold.",
                }
            )

    qc.append(
        {
            "severity": "critical",
            "endpoint": "all",
            "issue": "model validation metrics are not present in the Guild summary folder",
            "recommendation": "Cross-reference the separate validation table containing model formulation, CV method, R2/RMSE or AUC/TSS, and Moran's I.",
        }
    )
    qc.append(
        {
            "severity": "moderate",
            "endpoint": "all",
            "issue": "valid-pixel denominators differ among endpoints",
            "recommendation": "Report total grid pixels, land pixels, complete-predictor pixels, and define fractions as in_distribution_pixels / valid_pixels.",
        }
    )
    qc.append(
        {
            "severity": "moderate",
            "endpoint": "all",
            "issue": "old *_ridge_ood_pixel_summary.csv files do not store their threshold mapping",
            "recommendation": "Do not merge these legacy risk-class fractions into the chi90 table unless class boundaries are added.",
        }
    )
    return qc


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    main_rows = [load_primary(endpoint, cfg) for endpoint, cfg in ENDPOINTS.items()]
    main_columns = [
        "endpoint", "guild", "response_type", "model_label", "model_type",
        "predictors", "n_training_rows", "distance_method", "distance_space",
        "ridge_lambda", "threshold_method", "threshold_profile",
        "threshold_probability", "degrees_of_freedom", "chi_square_threshold_D2",
        "distance_threshold", "valid_pixels", "in_distribution_pixels",
        "in_distribution_fraction", "extrapolated_fraction", "risk_ratio_definition",
        "fraction_denominator",
        "land_mask_status", "source_file",
    ]
    write_csv(
        OUT_DIR / "Supplementary_Table_Guild_Extrapolation_chi90.csv",
        main_rows,
        main_columns,
    )

    sensitivity_rows = load_sensitivity_rows()
    sensitivity_columns = [
        "endpoint", "guild", "response_type", "profile", "threshold_type",
        "probability", "distance_threshold", "valid_pixels",
        "in_distribution_pixels", "in_distribution_fraction",
        "extrapolated_fraction", "source_file",
    ]
    write_csv(
        OUT_DIR / "Supplementary_Table_Guild_Extrapolation_sensitivity.csv",
        sensitivity_rows,
        sensitivity_columns,
    )

    merged_rows = merge_validation_and_ood(main_rows)
    merged_columns = main_columns + [
        "validated_model_type", "validated_feature_mode",
        "validated_environmental_factors", "validated_spatial_terms",
        "validation_n", "best_validation_method", "best_R2", "best_RMSE",
        "best_AUC", "best_TSS", "spatial_cv_method", "spatial_cv_R2",
        "moran_observed_I", "moran_expected_I", "moran_p_value",
        "predictor_consistency", "spatial_terms_consistency", "reviewer_action",
    ]
    write_csv(
        OUT_DIR / "Supplementary_Table_Guild_Validation_and_Extrapolation.csv",
        merged_rows,
        merged_columns,
    )

    qc_rows = build_qc_rows(main_rows, sensitivity_rows)
    for row in merged_rows:
        if row["predictor_consistency"] != "MATCH":
            qc_rows.append(
                {
                    "severity": "critical",
                    "endpoint": row["endpoint"],
                    "issue": f"validation/OOD environmental predictors: {row['predictor_consistency']}",
                    "recommendation": row["reviewer_action"],
                }
            )
        if row["spatial_terms_consistency"] != "NONE_REPORTED":
            qc_rows.append(
                {
                    "severity": "critical",
                    "endpoint": row["endpoint"],
                    "issue": "spatial terms are reported in validation but absent from OOD predictor space",
                    "recommendation": "Clarify whether spatial terms belong to the final projected model; if yes, include their global layers in prediction and OOD calculations.",
                }
            )
    write_csv(
        OUT_DIR / "Supplementary_Table_Guild_Extrapolation_QC.csv",
        qc_rows,
        ["severity", "endpoint", "issue", "recommendation"],
    )

    notes = [
        "Primary table: one row per guild endpoint using the theoretical chi-square 0.90 threshold.",
        "Distance threshold is sqrt(D2); chi_square_threshold_D2 is the squared-distance threshold.",
        "in_distribution_fraction = in_distribution_pixels / valid_pixels.",
        "valid_pixels means grid cells with complete final predictors; the land-mask definition must be added to the Methods/table note.",
        "The Guild summary folder does not contain model-validation metrics; keep those in a separate validation table.",
        "Do not include macOS ._ files or legacy *_ridge_ood_pixel_summary.csv files in the primary table without threshold definitions.",
    ]
    (OUT_DIR / "README.txt").write_text("\n".join(notes) + "\n", encoding="utf-8")

    print(f"Wrote {len(main_rows)} primary rows to {OUT_DIR}")
    print(f"Wrote {len(sensitivity_rows)} sensitivity rows")
    print(f"Wrote {len(qc_rows)} QC rows")


if __name__ == "__main__":
    main()
