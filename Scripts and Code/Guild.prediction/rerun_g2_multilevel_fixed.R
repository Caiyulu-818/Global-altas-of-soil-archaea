#!/usr/bin/env Rscript

# Rerun G2 output with updated global CSVs.
# Server:
#   nohup Rscript rerun_g2_multilevel_fixed.R > g2_rerun_fixed_20260628.log 2>&1 &

args_all <- commandArgs(FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
script <- file.path(script_dir, "rerun_guild_multilevel_fixed.R")
status <- system2("Rscript", c(script, "g2"))
if (!identical(status, 0L)) stop("g2 rerun failed with exit status ", status)
