#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import td_gfm_add_ridge_90_profile_core as core
core.DATASET_OVERRIDE = "TD"
core.SCENARIO_OVERRIDE = "cropland2020_SSP245"
core.PROFILE_MODE = "empirical_p90"
if __name__ == "__main__":
    core.main()
