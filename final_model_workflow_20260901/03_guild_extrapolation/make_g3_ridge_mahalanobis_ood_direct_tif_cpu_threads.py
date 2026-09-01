#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import make_g1_g5_ridge_mahalanobis_ood_direct_tif_cpu_threads as core

core.MODELS_TO_RUN = ["g3", "g3_abundance"]

if __name__ == "__main__":
    core.main()
