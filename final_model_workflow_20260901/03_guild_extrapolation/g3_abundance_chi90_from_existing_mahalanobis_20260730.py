#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""G3-abundance chi90 risk products from an existing Mahalanobis TIFF."""

import os
import chi90_risk_ratio_from_existing_md_core as core

RESULT_DIR = "/public/home/ylwang/Gisdata/TC_New/lucy/uncertainty/g3_results"

core.run({
    "model": "g3_abundance",
    "abundance": True,
    "features": ["ref_1", "ndvimean"],
    "df": 2,
    "chi_square_d2": 4.605170185988092,  # qchisq(0.90, df = 2)
    "md_tif": os.path.join(
        RESULT_DIR, "ridge_mahalanobis_ood_tif_direct1", "mahalanobis_distance_abundance.tif"
    ),
    "output_dir": os.path.join(
        RESULT_DIR, "chi90_from_existing_mahalanobis_abundance_20260730"
    ),
})
