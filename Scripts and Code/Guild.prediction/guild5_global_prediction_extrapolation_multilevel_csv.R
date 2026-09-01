#!/usr/bin/env Rscript

# Fixed version: the sourced core script uses scaled environmental predictors
# but replaces standardized Lat/Lon with raw env.meta.red*.csv coordinates
# for nearest-sample geographic distance.
guild_id <- "g5"
args_all <- commandArgs(FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
source(file.path(script_dir, "guild_global_prediction_extrapolation_multilevel_csv.R"), chdir = TRUE)
