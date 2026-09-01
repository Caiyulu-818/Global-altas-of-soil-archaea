#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import td_gfm_add_ridge_90_profile_core as core
core.DATASET_OVERRIDE = "GFM"
core.SCENARIO_OVERRIDE = "cropland2100_SSP126"
core.PROFILE_MODE = "chisq90"
if __name__ == "__main__":
    core.main()
