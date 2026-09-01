#!/usr/bin/env Rscript

# Entry point for G2. The shared core script contains the updated 20260628
# global environmental CSV paths and writes to pca_multilevel1.
guild_id <- "g2"
args_all <- commandArgs(FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
source(file.path(script_dir, "guild_global_prediction_extrapolation_multilevel_csv.R"), chdir = TRUE)
