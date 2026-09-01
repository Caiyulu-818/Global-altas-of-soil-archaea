#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import make_g1_g5_ridge_mahalanobis_ood_direct_tif_fast as core

core.MODELS_TO_RUN = ["g5"]

if __name__ == "__main__":
    core.main()
