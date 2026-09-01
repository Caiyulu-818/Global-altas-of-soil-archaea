#!/usr/bin/env Rscript

# Extra continuous abundance prediction for G3.
# Main G3 map remains presence/absence; this script writes abundance-specific
# outputs such as prediction_abundance.csv and prediction_abundance_reliable_*.csv.
guild_id <- "g3_abundance"
args_all <- commandArgs(FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
source(file.path(script_dir, "guild_global_prediction_extrapolation_multilevel_csv.R"), chdir = TRUE)
