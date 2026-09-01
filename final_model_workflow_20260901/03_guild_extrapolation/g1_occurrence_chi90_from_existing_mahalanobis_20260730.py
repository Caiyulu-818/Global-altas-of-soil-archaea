#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""G1 occurrence chi90 risk products from an existing Mahalanobis TIFF."""

import os
import chi90_risk_ratio_from_existing_md_core as core

RESULT_DIR = "/public/home/ylwang/Gisdata/TC_New/lucy/uncertainty/g1_results"

core.run({
    "model": "g1_occurrence",
    "abundance": False,
    "features": ["ref_1", "ndvimean", "soc"],
    "df": 3,
    "chi_square_d2": 6.251388631170325,  # qchisq(0.90, df = 3)
    "md_tif": os.path.join(
        RESULT_DIR, "ridge_mahalanobis_ood_tif_direct1", "mahalanobis_distance.tif"
    ),
    "output_dir": os.path.join(
        RESULT_DIR, "chi90_from_existing_mahalanobis_occurrence_20260730"
    ),
})
