#!/usr/bin/env Rscript

args_all <- commandArgs(FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
run_target_override <- "ALL"
source(file.path(script_dir, "td_gfm_global_scenario_extrapolation_csv.R"), chdir = TRUE, local = FALSE)
