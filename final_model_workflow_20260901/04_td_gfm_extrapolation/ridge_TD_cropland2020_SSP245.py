#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import td_gfm_ridge_mahalanobis_ood_direct_tif as core

core.DATASET_OVERRIDE = "TD"
core.SCENARIO_OVERRIDE = "cropland2020_SSP245"

if __name__ == "__main__":
    core.main()
