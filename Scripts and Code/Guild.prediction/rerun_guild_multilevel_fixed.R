#!/usr/bin/env Rscript

# Rerun fixed multi-level global predictions for guild models.
#
# This 20260628 version uses the resample/csv global environmental layers and
# writes all CSV outputs to:
#   /public/home/ylwang/Gisdata/TC_New/lucy/uncertainty/g*_results/pca_multilevel1
#
# Usage examples on server:
#   Rscript rerun_guild_multilevel_fixed.R
#   Rscript rerun_guild_multilevel_fixed.R g1 g1_abundance
#   Rscript rerun_guild_multilevel_fixed.R g5
#
# Set CLEAN_OLD_OUTPUTS=FALSE only if you deliberately want to append/check
# without clearing the output directory first. For normal reruns, leave it TRUE.

base_dir <- "/public/home/ylwang/Gisdata/TC_New/lucy"
uncertainty_dir <- file.path(base_dir, "uncertainty")
output_dirname <- "pca_multilevel1"
expected_csv_rows <- 36000L
clean_old_outputs <- toupper(Sys.getenv("CLEAN_OLD_OUTPUTS", unset = "TRUE")) != "FALSE"

default_models <- c("g1", "g1_abundance", "g2", "g3", "g3_abundance", "g4", "g5")
args <- commandArgs(trailingOnly = TRUE)
models_to_run <- if (length(args) > 0) unlist(strsplit(paste(args, collapse = ","), ",")) else default_models
models_to_run <- trimws(tolower(models_to_run))
models_to_run <- models_to_run[nzchar(models_to_run)]

valid_models <- default_models
bad_models <- setdiff(models_to_run, valid_models)
if (length(bad_models) > 0) {
  stop("Unknown model(s): ", paste(bad_models, collapse = ", "),
       "\nValid models are: ", paste(valid_models, collapse = ", "))
}

model_to_guild <- function(model) sub("_abundance$", "", model)
prediction_file_name <- function(model) {
  if (grepl("_abundance$", model)) "prediction_abundance.csv" else "prediction.csv"
}

args_all <- commandArgs(FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
code_dir_env <- Sys.getenv("GUILD_CODE_DIR", unset = "")
if (nzchar(code_dir_env)) {
  script_dir <- normalizePath(code_dir_env, mustWork = FALSE)
}

core_script <- file.path(script_dir, "guild_global_prediction_extrapolation_multilevel_csv.R")
if (!file.exists(core_script)) {
  stop(
    "Missing core script: ", core_script, "\n",
    "Put guild_global_prediction_extrapolation_multilevel_csv.R in the same directory as this rerun script, ",
    "or run with GUILD_CODE_DIR=/path/to/scripts."
  )
}

count_lines <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  out <- system2("wc", c("-l", path), stdout = TRUE, stderr = TRUE)
  as.integer(strsplit(trimws(out[1]), "\\s+")[[1]][1])
}

if (clean_old_outputs) {
  guilds_to_clean <- unique(vapply(models_to_run, model_to_guild, FUN.VALUE = character(1)))
  for (guild in guilds_to_clean) {
    out_dir <- file.path(uncertainty_dir, paste0(guild, "_results"), output_dirname)
    message("[", guild, "] cleaning old CSV outputs in: ", out_dir)
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    unlink(file.path(out_dir, "*.csv"))
  }
}

for (model in models_to_run) {
  guild <- model_to_guild(model)
  out_dir <- file.path(uncertainty_dir, paste0(guild, "_results"), output_dirname)
  pred_csv <- file.path(out_dir, prediction_file_name(model))

  message("\n==============================")
  message("Running model: ", model)
  message("Output dir   : ", out_dir)
  message("==============================")

  status <- system2("Rscript", c(core_script, model))
  if (!identical(status, 0L)) {
    stop(model, " failed with exit status ", status)
  }

  n_rows <- count_lines(pred_csv)
  message("[", model, "] ", basename(pred_csv), " rows: ", n_rows)
  if (is.na(n_rows)) {
    stop(model, " did not create expected prediction CSV: ", pred_csv)
  }
  if (n_rows != expected_csv_rows) {
    stop(model, " prediction CSV incomplete: ", n_rows, " rows, expected ", expected_csv_rows)
  }
}

message("\nAll requested guild models finished successfully.")
