# Fast Ridge-Mahalanobis Direct-to-TIF for G1-G5

This version is faster than the previous direct-to-TIF script because it can:

1. Read `pca_multilevel1_tif/prediction*.tif` directly and skip the expensive `prediction.csv -> memmap` step.
2. Read environmental predictor values from GeoTIFF rows instead of parsing huge environmental CSV rows.

It does not rerun RF prediction.

## Input

Prediction:

- Preferred fast input:
  - `g*_results/pca_multilevel1_tif/prediction.tif`
  - `g*_results/pca_multilevel1_tif/prediction_abundance.tif`
- Fallback input:
  - `g*_results/pca_multilevel1/prediction.csv`
  - `g*_results/pca_multilevel1/prediction_abundance.csv`

Environmental predictors:

- `/public/home/ZhangFZ/Gisdata/resample/tif/New_*.tif`

## Output

Outputs are written to:

```text
g*_results/ridge_mahalanobis_ood_tif_direct2_fast/
```

This avoids overwriting previous results in `ridge_mahalanobis_ood_tif_direct1` or `ridge_mahalanobis_ood_tif_direct2`.

## Run

```bash
python make_g4_ridge_mahalanobis_ood_direct_tif_fast.py > g4_ridge_fast.log 2>&1
```

Other guilds:

```bash
python make_g1_ridge_mahalanobis_ood_direct_tif_fast.py > g1_ridge_fast.log 2>&1
python make_g2_ridge_mahalanobis_ood_direct_tif_fast.py > g2_ridge_fast.log 2>&1
python make_g3_ridge_mahalanobis_ood_direct_tif_fast.py > g3_ridge_fast.log 2>&1
python make_g4_ridge_mahalanobis_ood_direct_tif_fast.py > g4_ridge_fast.log 2>&1
python make_g5_ridge_mahalanobis_ood_direct_tif_fast.py > g5_ridge_fast.log 2>&1
```

## Important switches

Inside `make_g1_g5_ridge_mahalanobis_ood_direct_tif_fast.py`:

```python
USE_PREDICTION_TIF_IF_AVAILABLE = True
USE_FEATURE_TIFS_FOR_MD = True
```

If you are not sure that predictor TIFs are on the same scale as the prediction CSV inputs, set:

```python
USE_FEATURE_TIFS_FOR_MD = False
```

That falls back to the original global environmental CSV parsing path.
