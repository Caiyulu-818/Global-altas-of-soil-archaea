# Final modelling and global projection workflow

This folder is a clean, self-contained archive of the final Guild 1--5, TD and GFM workflow. All files here are copies; the original scripts and existing server outputs were not changed or deleted.

## Folder structure

```text
01_model_validation/
  Model screening, cross-validation and Moran's I assessment
02_global_prediction/
  Final ranger global prediction CSV scripts
03_guild_extrapolation/
  Guild ridge-Mahalanobis environmental extrapolation and chi-square 0.90 post-processing
04_td_gfm_extrapolation/
  TD/GFM scenario prediction and ridge-Mahalanobis extrapolation
05_supplementary_tables/
  Supplementary-table compilation scripts
```

## 1. Guild model validation

The scripts in `01_model_validation` are the analysis scripts used for model screening and validation:

- `red_abundance_model_selection_and_moran.R`: continuous abundance candidates for Guilds 1--5.
- `red_occurrence_model_selection_and_moran.R`: occurrence-suitability candidates; run with `1` or `3`.
- `red_spatial_cv_sensitivity.R`: spatial-validation sensitivity analysis.
- `td_gfm_model_selection_and_moran.R`: TD and GFM continuous-response modelling.
- `independent_rf_moran_audit.R`: independent RF/Moran audit.

The model-selection scripts compare the specified response transformations, predictor sets and validation schemes. The reported Moran's I values in these scripts are calculated from the final full-data fitted residuals. They should therefore be described as a residual spatial-autocorrelation diagnostic, not as Moran's I calculated from spatial out-of-fold residuals.

The validation scripts retain the original local input paths. Before running them on a server, update their input/output roots to the corresponding server paths. The global-prediction scripts in Section 2 already use the server paths shown below.

## 2. Final Guild global predictions

The current global environmental layers are read from:

```text
/public/home/ZhangFZ/Gisdata/resample/csv
```

Training data are read from:

```text
/public/home/ylwang/Gisdata/TC_New/lucy/env.meta.red1_scaled.csv
/public/home/ylwang/Gisdata/TC_New/lucy/env.meta.red2_scaled.csv
/public/home/ylwang/Gisdata/TC_New/lucy/env.meta.red3_scaled.csv
/public/home/ylwang/Gisdata/TC_New/lucy/env.meta.red4_scaled.csv
/public/home/ylwang/Gisdata/TC_New/lucy/env.meta.red5_scaled.csv
```

Outputs are written to:

```text
/public/home/ylwang/Gisdata/TC_New/lucy/uncertainty/g*_results/pca_multilevel1
```

The final predictor sets in the archived core are:

| Endpoint | Predictors | Response |
|---|---|---|
| G1 | `ref_1`, `ndvimean`, `soc` | occurrence suitability |
| G1 abundance | `ref_1`, `ndvimean`, `Elevation` | exploratory abundance |
| G2 | `ph`, `ref_1`, `ndvimean` | abundance |
| G3 | `ref_1`, `ndvimean` | occurrence suitability |
| G3 abundance | `ndvimean`, `ref_1` | exploratory abundance |
| G4 | `ph`, `soc`, `ref_1`, `biomass` | abundance |
| G5 | `ref_1`, `ndvimean`, `biomass`, `ph` | abundance |

The G1 abundance configuration in this archive explicitly excludes ECE because the ECE global layer was incomplete in the northern hemisphere. The generic core therefore uses the same no-ECE predictor set as the dedicated `fast_g1_abundance_noECE_prediction_only.R` script.

Run one model or several models with:

```bash
Rscript rerun_guild_multilevel_fixed.R g1
Rscript rerun_guild_multilevel_fixed.R g2
Rscript rerun_guild_multilevel_fixed.R g3 g3_abundance
Rscript rerun_guild_multilevel_fixed.R g4
Rscript rerun_guild_multilevel_fixed.R g5
Rscript rerun_guild_multilevel_fixed.R g1_abundance
```

`rerun_g1_multilevel_fixed.R` and the other Guild-specific wrappers are also provided. The prediction-only no-ECE script can be used when only the G1 abundance prediction CSV needs to be rebuilt:

```bash
Rscript fast_g1_abundance_noECE_prediction_only.R
```

By default it writes `prediction_abundancermece.csv`; use `OUTPUT_TAG=_abundance` if the standard `prediction_abundance.csv` name is required by a downstream script.

## 3. Guild environmental extrapolation

The recommended fast scripts are in `03_guild_extrapolation`:

```bash
python make_g1_ridge_mahalanobis_ood_direct_tif_cpu_threads.py
python make_g2_ridge_mahalanobis_ood_direct_tif_cpu_threads.py
python make_g3_ridge_mahalanobis_ood_direct_tif_cpu_threads.py
python make_g4_ridge_mahalanobis_ood_direct_tif_cpu_threads.py
python make_g5_ridge_mahalanobis_ood_direct_tif_cpu_threads.py
```

They read the completed prediction CSVs from `pca_multilevel1`, calculate ridge-regularized Mahalanobis distances in the standardized final-predictor space, and write GeoTIFFs directly to:

```text
/public/home/ylwang/Gisdata/TC_New/lucy/uncertainty/g*_results/ridge_mahalanobis_ood_tif_direct2
```

The default sensitivity profiles are empirical training-distance thresholds at 97.5%, 99.0% and 99.5%. The `environmental_extrapolation_risk_ratio.tif` is a continuous layer: values above 1 are beyond the reference threshold used by that script. The `environmental_extrapolation_mask_<profile>.tif` files are binary threshold masks, and `prediction_reliable_<profile>.tif` files are the prediction maps after applying those environmental masks.

The chi-square 0.90 scripts in the same folder are post-processing scripts. They reuse an existing Mahalanobis-distance TIF and do not rerun the RF model or recompute the global prediction:

```bash
python g2_chi90_from_existing_mahalanobis_20260730.py
python g4_chi90_from_existing_mahalanobis_20260730.py
python g5_chi90_from_existing_mahalanobis_20260730.py
python g1_occurrence_chi90_from_existing_mahalanobis_20260730.py
python g1_abundance_chi90_risk_ratio_fast_20260730.py
python g3_occurrence_chi90_from_existing_mahalanobis_20260730.py
python g3_abundance_chi90_from_existing_mahalanobis_20260730.py
```

These use the theoretical threshold `sqrt(qchisq(0.90, df = p))`, where `p` is the number of final predictors. This is distinct from the empirical 90th-percentile threshold.

## 4. TD/GFM scenario extrapolation

The eight `ridge_TD_*.py` and `ridge_GFM_*.py` scripts in `04_td_gfm_extrapolation` run the four scenario-specific TD/GFM extrapolation jobs. They use the TD/GFM-specific predictor sets and cropland scenario layer, and write outputs under:

```text
/public/home/ylwang/Gisdata/TC_New/lucy/cropland/scenario_extrapolation
```

The `add_TD_*.py` and `add_GFM_*.py` scripts add empirical p90 or theoretical chi-square 0.90 profiles from existing Mahalanobis-distance TIFs. They do not rerun RF prediction.

## 5. Supplementary tables

Use the two scripts in `05_supplementary_tables` to compile Guild and TD/GFM extrapolation summaries. Check their input roots before running on a different machine.

## Reproducibility notes

- `Lat` and `Lon` are used for sample coordinates and are not included as environmental predictors in the final Guild predictor sets.
- The global environmental layers and training predictors must use the same variable definitions and scaling convention.
- `prediction.csv` and `prediction_abundance.csv` are inputs to extrapolation; the extrapolation scripts do not alter them.
- A complete output must contain the expected 36,000 longitude rows before CSV-to-TIF conversion or downstream masking.
