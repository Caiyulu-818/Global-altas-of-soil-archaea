# TD/GFM Ridge-Mahalanobis 90% Profiles

These scripts add two additional 90% OOD profiles to existing TD/GFM ridge outputs.
They do not rerun RF prediction and do not recompute global Mahalanobis distance.

Required existing files for each dataset/scenario:

- `scenario_extrapolation/<TD|GFM>/<scenario>/prediction.csv`
- `scenario_extrapolation/<TD|GFM>/<scenario>/ridge_mahalanobis_ood_tif_direct/mahalanobis_distance.tif`

Output directory:

- `scenario_extrapolation/<TD|GFM>/<scenario>/ridge_mahalanobis_ood_tif_direct/`

## Profile 1: empirical p90

Threshold:

```text
90th percentile of training-sample ridge-Mahalanobis distances
```

Example:

```bash
nohup python add_TD_cropland2020_SSP245_empirical_p90.py > TD_2020_SSP245_empirical_p90.log 2>&1 &
```

Outputs:

- `environmental_extrapolation_mask_p90.tif`
- `prediction_ridge_reliable_p90.tif`
- `<dataset>_<scenario>_ridge_mahalanobis_threshold_p90.csv`
- `<dataset>_<scenario>_ridge_p90_pixel_summary.csv`

## Profile 2: chi-square 0.90

Threshold:

```text
sqrt(qchisq(0.9, df = number_of_predictors))
```

The square root is used because existing `mahalanobis_distance.tif` stores `sqrt(D^2)`.

Example:

```bash
nohup python add_TD_cropland2020_SSP245_chisq90.py > TD_2020_SSP245_chisq90.log 2>&1 &
```

Outputs:

- `environmental_extrapolation_mask_chi90.tif`
- `prediction_ridge_reliable_chi90.tif`
- `<dataset>_<scenario>_ridge_mahalanobis_threshold_chi90.csv`
- `<dataset>_<scenario>_ridge_chi90_pixel_summary.csv`

## Mask values

- `0` = in-distribution / retained
- `1` = environmental outlier / masked
- `255` = NoData

## All wrapper scripts

TD empirical p90:

```bash
nohup python add_TD_cropland2020_SSP245_empirical_p90.py > TD_2020_SSP245_empirical_p90.log 2>&1 &
nohup python add_TD_cropland2100_SSP126_empirical_p90.py > TD_2100_SSP126_empirical_p90.log 2>&1 &
nohup python add_TD_cropland2100_SSP245_empirical_p90.py > TD_2100_SSP245_empirical_p90.log 2>&1 &
nohup python add_TD_cropland2100_SSP370_empirical_p90.py > TD_2100_SSP370_empirical_p90.log 2>&1 &
```

TD chi-square 0.90:

```bash
nohup python add_TD_cropland2020_SSP245_chisq90.py > TD_2020_SSP245_chisq90.log 2>&1 &
nohup python add_TD_cropland2100_SSP126_chisq90.py > TD_2100_SSP126_chisq90.log 2>&1 &
nohup python add_TD_cropland2100_SSP245_chisq90.py > TD_2100_SSP245_chisq90.log 2>&1 &
nohup python add_TD_cropland2100_SSP370_chisq90.py > TD_2100_SSP370_chisq90.log 2>&1 &
```

GFM empirical p90:

```bash
nohup python add_GFM_cropland2020_SSP245_empirical_p90.py > GFM_2020_SSP245_empirical_p90.log 2>&1 &
nohup python add_GFM_cropland2100_SSP126_empirical_p90.py > GFM_2100_SSP126_empirical_p90.log 2>&1 &
nohup python add_GFM_cropland2100_SSP245_empirical_p90.py > GFM_2100_SSP245_empirical_p90.log 2>&1 &
nohup python add_GFM_cropland2100_SSP370_empirical_p90.py > GFM_2100_SSP370_empirical_p90.log 2>&1 &
```

GFM chi-square 0.90:

```bash
nohup python add_GFM_cropland2020_SSP245_chisq90.py > GFM_2020_SSP245_chisq90.log 2>&1 &
nohup python add_GFM_cropland2100_SSP126_chisq90.py > GFM_2100_SSP126_chisq90.log 2>&1 &
nohup python add_GFM_cropland2100_SSP245_chisq90.py > GFM_2100_SSP245_chisq90.log 2>&1 &
nohup python add_GFM_cropland2100_SSP370_chisq90.py > GFM_2100_SSP370_chisq90.log 2>&1 &
```
