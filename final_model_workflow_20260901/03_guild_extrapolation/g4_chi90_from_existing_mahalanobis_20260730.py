#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""G4 chi90 risk products from an existing Mahalanobis TIFF."""

import os
import chi90_risk_ratio_from_existing_md_core as core

RESULT_DIR = "/public/home/ylwang/Gisdata/TC_New/lucy/uncertainty/g4_results"

core.run({
    "model": "g4",
    "abundance": False,
    "features": ["ph", "soc", "ref_1", "biomass"],
    "df": 4,
    "chi_square_d2": 7.779440339734858,  # qchisq(0.90, df = 4)
    "md_tif": os.path.join(
        RESULT_DIR, "ridge_mahalanobis_ood_tif_direct1", "mahalanobis_distance.tif"
    ),
    "output_dir": os.path.join(RESULT_DIR, "chi90_from_existing_mahalanobis_20260730"),
})
