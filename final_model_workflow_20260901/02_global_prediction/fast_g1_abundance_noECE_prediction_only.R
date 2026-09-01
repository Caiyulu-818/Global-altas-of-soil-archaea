#!/usr/bin/env Rscript

# ============================================================
# Fast G1 abundance global prediction only, no ECE
# ============================================================
#
# Purpose:
#   Rebuild only the G1 abundance prediction CSV after removing ECE.
#   This script deliberately skips Mahalanobis, nearest-sample distance,
#   risk masks, and reliable-prediction CSVs. Use the Python ridge/OOD
#   scripts after this if you need reliable TIFs.
#
# Main output:
#   /public/home/ylwang/Gisdata/TC_New/lucy/uncertainty/g1_results/pca_multilevel1/
#     prediction_abundancermece.csv
#
# Speed controls can be overridden from shell:
#   NTREE_MAIN=1000 PREDICT_THREADS=16 COL_BLOCK_SIZE=20 Rscript fast_g1_abundance_noECE_prediction_only.R
#
# If you want to overwrite the standard abundance name instead:
#   OUTPUT_TAG=_abundance Rscript fast_g1_abundance_noECE_prediction_only.R

suppressPackageStartupMessages({
  library(ranger)
  library(data.table)
})

set.seed(123)

# -----------------------------
# 0. User settings
# -----------------------------
base_dir <- "/public/home/ylwang/Gisdata/TC_New/lucy"
uncertainty_dir <- file.path(base_dir, "uncertainty")
env_csv_root <- "/public/home/ZhangFZ/Gisdata/resample/csv"

training_file <- file.path(base_dir, "env.meta.red1_scaled.csv")
coord_file <- file.path(base_dir, "env.meta.red1.csv")
path_prefix <- file.path(uncertainty_dir, "g1_results", "pca_multilevel1")
if (!dir.exists(path_prefix)) dir.create(path_prefix, recursive = TRUE)

features <- c("ref_1", "ndvimean", "Elevation")
target <- "Guild.1"
output_tag <- Sys.getenv("OUTPUT_TAG", "_abundancermece")

global_files <- c(
  Elevation = file.path(env_csv_root, "New_dem.csv"),
  ref_1    = file.path(env_csv_root, "New_Nadir_Reflectance_1.csv"),
  ndvimean = file.path(env_csv_root, "New_NDVImean.csv")
)
global_files <- global_files[features]

ntree_main <- as.integer(Sys.getenv("NTREE_MAIN", "1000"))
predict_threads <- as.integer(Sys.getenv("PREDICT_THREADS", "16"))
col_block_size <- as.integer(Sys.getenv("COL_BLOCK_SIZE", "20"))
model_mtry <- as.integer(Sys.getenv("MTRY", "3"))

setDTthreads(max(1L, predict_threads))

tagged_name <- function(prefix, suffix = ".csv") {
  paste0(prefix, output_tag, suffix)
}

append_csv_row <- function(x, file) {
  fwrite(data.table(t(x)), file, append = TRUE, col.names = FALSE, na = "NA")
}

looks_like_real_lonlat <- function(lat, lon) {
  all(is.finite(lat), is.finite(lon)) &&
    min(lat, na.rm = TRUE) >= -90 &&
    max(lat, na.rm = TRUE) <= 90 &&
    min(lon, na.rm = TRUE) >= -180 &&
    max(lon, na.rm = TRUE) <= 180 &&
    diff(range(lat, na.rm = TRUE)) > 10 &&
    diff(range(lon, na.rm = TRUE)) > 10
}

# -----------------------------
# 1. Load training data
# -----------------------------
cat("Loading training data...\n")
env.meta <- fread(training_file)
if (!"biome" %in% names(env.meta)) env.meta[, biome := factor(1)]

if (!looks_like_real_lonlat(env.meta$Lat, env.meta$Lon)) {
  cat("Lat/Lon in scaled training file look standardized; replacing with raw coordinates from:\n")
  cat(coord_file, "\n")
  coord.meta <- fread(coord_file)
  if (!all(c("Lat", "Lon") %in% names(coord.meta))) {
    stop("Raw coordinate file does not contain Lat/Lon: ", coord_file)
  }
  if (nrow(coord.meta) != nrow(env.meta)) {
    stop("Cannot safely replace Lat/Lon because row counts differ.")
  }
  env.meta[, Lat := coord.meta$Lat]
  env.meta[, Lon := coord.meta$Lon]
}

needed_cols <- c(features, target, "Lat", "Lon", "biome")
missing_cols <- setdiff(needed_cols, names(env.meta))
if (length(missing_cols) > 0) {
  stop("Missing columns in training data: ", paste(missing_cols, collapse = ", "))
}

train_df <- na.omit(as.data.frame(env.meta[, ..needed_cols]))
train_df$biome <- as.factor(train_df$biome)

cat("Training rows:", nrow(train_df), "\n")
cat("Guild: g1_abundance_noECE\n")
cat("Features:", paste(features, collapse = ", "), "\n")
cat("Output tag:", output_tag, "\n")
cat("ntree_main:", ntree_main, "\n")
cat("predict_threads:", predict_threads, "\n")
cat("col_block_size:", col_block_size, "\n")

# -----------------------------
# 2. Train RF_shift_log1p model
# -----------------------------
x_train <- train_df[, features, drop = FALSE]
y_train <- train_df[[target]]
baseline <- min(y_train, na.rm = TRUE)
shifted <- y_train - baseline

cat("Training main ranger model...\n")
main_rf <- ranger(
  .y ~ .,
  data = data.frame(.y = log1p(shifted), x_train, check.names = FALSE),
  num.trees = ntree_main,
  mtry = max(1L, min(ncol(x_train), model_mtry)),
  min.node.size = 5,
  importance = "permutation",
  num.threads = predict_threads,
  seed = 123
)

importance_file <- file.path(path_prefix, "g1_abundance_noECE_ranger_importance.csv")
fwrite(
  data.frame(
    variable = names(main_rf$variable.importance),
    importance = as.numeric(main_rf$variable.importance)
  ),
  importance_file
)

predict_abundance <- function(newdata) {
  pr <- predict(main_rf, data = newdata, num.threads = predict_threads)$predictions
  as.numeric(expm1(pr) + baseline)
}

# -----------------------------
# 3. Load global feature CSVs
# -----------------------------
cat("Loading global environmental CSV matrices...\n")
missing_layers <- setdiff(features, names(global_files))
if (length(missing_layers) > 0) {
  stop("Missing global file paths for features: ", paste(missing_layers, collapse = ", "))
}

global_layers <- lapply(features, function(v) {
  cat("  reading", v, "from", global_files[[v]], "\n")
  fread(global_files[[v]])
})
names(global_layers) <- features

dims <- lapply(global_layers, dim)
dim_text <- vapply(dims, paste, collapse = "x", FUN.VALUE = character(1))
if (length(unique(dim_text)) != 1) {
  stop("Global layer dimensions do not match: ", paste(names(dim_text), dim_text, collapse = "; "))
}

n_lat <- nrow(global_layers[[1]])
n_lon <- ncol(global_layers[[1]])
cat("Global grid:", n_lon, "columns x", n_lat, "rows =",
    format(n_lon * n_lat, scientific = FALSE), "pixels\n")

metadata_file <- file.path(path_prefix, "g1_abundance_noECE_fast_prediction_settings.csv")
fwrite(
  data.frame(
    n_lon = n_lon,
    n_lat = n_lat,
    features = paste(features, collapse = ";"),
    target = target,
    model_kind = "RF_shift_log1p",
    baseline = baseline,
    ntree_main = ntree_main,
    mtry = max(1L, min(length(features), model_mtry)),
    predict_threads = predict_threads,
    col_block_size = col_block_size,
    output_tag = output_tag,
    output_dir = path_prefix
  ),
  metadata_file
)

# -----------------------------
# 4. Prediction-only global loop
# -----------------------------
file_prediction <- file.path(path_prefix, tagged_name("prediction"))
if (file.exists(file_prediction)) {
  cat("Removing only existing prediction file:\n")
  cat(file_prediction, "\n")
  invisible(file.remove(file_prediction))
}

cat("Starting fast prediction-only CSV writing...\n")
start_time <- Sys.time()

for (i0 in seq(1L, n_lon, by = col_block_size)) {
  i1 <- min(i0 + col_block_size - 1L, n_lon)
  idx <- i0:i1
  n_block <- length(idx)

  mats <- lapply(features, function(v) {
    as.matrix(global_layers[[v]][, idx, with = FALSE])
  })
  names(mats) <- features

  env.block <- as.data.frame(
    setNames(lapply(mats, as.vector), features),
    check.names = FALSE
  )

  valid <- complete.cases(env.block)
  pred_vec <- rep(NA_real_, nrow(env.block))
  if (any(valid)) {
    pred_vec[valid] <- predict_abundance(env.block[valid, , drop = FALSE])
  }

  pred_mat <- matrix(pred_vec, nrow = n_lat, ncol = n_block)
  for (j in seq_len(n_block)) {
    append_csv_row(pred_mat[, j], file_prediction)
  }

  if (i1 %% 100L == 0L || i1 == n_lon) {
    elapsed_h <- as.numeric(difftime(Sys.time(), start_time, units = "hours"))
    done <- i1 / n_lon
    eta_h <- elapsed_h / max(done, 1e-9) - elapsed_h
    cat(
      "Progress:", round(done * 100, 2), "%",
      "| longitude columns:", i1, "/", n_lon,
      "| elapsed h:", round(elapsed_h, 2),
      "| ETA h:", round(eta_h, 2), "\n"
    )
  }
}

cat("Done. Final prediction CSV:\n")
cat(file_prediction, "\n")
cat("Settings:\n")
cat(metadata_file, "\n")
cat("Importance:\n")
cat(importance_file, "\n")
