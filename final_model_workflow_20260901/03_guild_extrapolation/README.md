# CPU-Threaded Ridge-Mahalanobis Direct-to-TIF for G1-G5

This is the archived final Guild extrapolation workflow. The G1 abundance
configuration uses `ref_1`, `ndvimean` and `Elevation`; ECE is intentionally
not used because the former ECE layer was incomplete in the northern
hemisphere.

This version keeps the same ridge-Mahalanobis calculation and output logic as your direct-to-TIF script, but adds GDAL thread settings for GeoTIFF writing/compression.

## What CPU_THREADS changes

Inside `make_g1_g5_ridge_mahalanobis_ood_direct_tif_cpu_threads.py`:

```python
CPU_THREADS = 8
GDAL_CACHEMAX_MB = 2048
```

You can increase `CPU_THREADS` to match the CPU resources requested on the node, for example:

```python
CPU_THREADS = 16
```

or:

```python
CPU_THREADS = 28
```

This mainly helps the GeoTIFF writing/compression step. The huge CSV parsing and prediction CSV-to-memmap step are still mostly I/O-bound and single-process.

## Best way to use more CPU

Run different guild scripts on different nodes or as separate jobs:

```bash
python make_g1_ridge_mahalanobis_ood_direct_tif_cpu_threads.py > g1_ridge_cpu.log 2>&1
python make_g2_ridge_mahalanobis_ood_direct_tif_cpu_threads.py > g2_ridge_cpu.log 2>&1
python make_g3_ridge_mahalanobis_ood_direct_tif_cpu_threads.py > g3_ridge_cpu.log 2>&1
python make_g4_ridge_mahalanobis_ood_direct_tif_cpu_threads.py > g4_ridge_cpu.log 2>&1
python make_g5_ridge_mahalanobis_ood_direct_tif_cpu_threads.py > g5_ridge_cpu.log 2>&1
```

G1 and G3 wrappers also run the abundance versions.

## Output

The output directory is unchanged from the attached script:

```text
g*_results/ridge_mahalanobis_ood_tif_direct2/
```

Existing files with the same names in that directory may be overwritten by GDAL when rerun.

## Chi-square 0.90 post-processing

The `*_chi90_*.py` scripts reuse an existing Mahalanobis-distance TIF and add
the theoretical `sqrt(qchisq(0.90, df = p))` threshold. They do not rerun the
RF model or modify the source prediction CSV. Their outputs are written to
separate chi-square-90 directories configured inside each script.
