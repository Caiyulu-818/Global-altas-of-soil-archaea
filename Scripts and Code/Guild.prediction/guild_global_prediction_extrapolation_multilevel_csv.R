#!/usr/bin/env Rscript

# ============================================================
# Guild global prediction + multi-level extrapolation risk
# CSV-only version: no terra, no raster, no gdal
# ============================================================
#
# Required R packages:
#   ranger
#   data.table
#   parallel
#   RANN
#
# Final common CSV outputs:
#   prediction.csv
#   mahalanobis_distance.csv
#   nearest_sample_distance.csv
#
# Final profile-specific CSV outputs:
#   environmental_extrapolation_mask_<profile>.csv
#   combined_extrapolation_risk_<profile>.csv
#   reliable_prediction_mask_<profile>.csv
#   prediction_reliable_<profile>.csv
#
# Profiles are strict / moderate / lenient, so one global run gives
# sensitivity outputs for different extrapolation-risk thresholds.
#
# Output orientation:
#   This follows your original script:
#   - one output row per longitude column
#   - each output row contains all latitude values
#   Therefore output dimensions are usually:
#   36000 rows x 18000 columns
#
# For local raster conversion/visualization, you may need to transpose the CSV
# depending on your plotting workflow.

suppressPackageStartupMessages({
  library(ranger)
  library(data.table)
  library(parallel)
})

if (!requireNamespace("RANN", quietly = TRUE)) {
  stop("Package RANN is required. Install on the server with: install.packages('RANN')")
}

set.seed(123)

# -----------------------------
# 0. User settings
# -----------------------------
if (!exists("guild_id")) {
  args <- commandArgs(trailingOnly = TRUE)
  guild_id <- if (length(args) >= 1) args[[1]] else "g2"
}
guild_id <- tolower(guild_id)

base_dir <- "/public/home/ylwang/Gisdata/TC_New/lucy"
uncertainty_dir <- file.path(base_dir, "uncertainty")
env_csv_root <- "/public/home/ZhangFZ/Gisdata/resample/csv"
training_files <- c(
  g1 = file.path(base_dir, "env.meta.red1_scaled.csv"),
  g2 = file.path(base_dir, "env.meta.red2_scaled.csv"),
  g3 = file.path(base_dir, "env.meta.red3_scaled.csv"),
  g4 = file.path(base_dir, "env.meta.red4_scaled.csv"),
  g5 = file.path(base_dir, "env.meta.red5_scaled.csv")
)
coord_files <- c(
  g1 = file.path(base_dir, "env.meta.red1.csv"),
  g2 = file.path(base_dir, "env.meta.red2.csv"),
  g3 = file.path(base_dir, "env.meta.red3.csv"),
  g4 = file.path(base_dir, "env.meta.red4.csv"),
  g5 = file.path(base_dir, "env.meta.red5.csv")
)

global_file_catalog <- c(
  Elevation = file.path(env_csv_root, "New_dem.csv"),
  evenness  = file.path(env_csv_root, "New_evenness_01_05_1km_uint16.csv"),
  northness = file.path(env_csv_root, "New_northness.csv"),
  ref_1     = file.path(env_csv_root, "New_Nadir_Reflectance_1.csv"),
  ECE       = file.path(env_csv_root, "New_Mean_of_annual_predicted_ECe_1980_2018.csv"),
  biomass   = file.path(env_csv_root, "New_soilbiomass_1km.csv"),
  ph        = file.path(env_csv_root, "New_Average_phh2o2.csv"),
  soc       = file.path(env_csv_root, "New_Average_soc2.csv"),
  ndvimean  = file.path(env_csv_root, "New_NDVImean.csv")
)

guild_configs <- list(
  # Scaled-limited PA alternative. Full-env PA models need extra global layer paths.
  g1 = list(
    output_subdir = "g1_results",
    training_file = training_files[["g1"]],
    coord_file = coord_files[["g1"]],
    output_tag = "",
    features = c("ref_1", "ndvimean", "soc"),
    target = "Guild.1",
    response_type = "presence_absence",
    model_kind = "RF_unweighted",
    mtry = NA_integer_
  ),
  # Extra continuous abundance map for G1. The main G1 map remains PA because
  # abundance CV performance was weak; this output is for visualization/supplement.
  g1_abundance = list(
    output_subdir = "g1_results",
    training_file = training_files[["g1"]],
    coord_file = coord_files[["g1"]],
    output_tag = "_abundance",
    features = c("ref_1", "ndvimean", "ECE", "Elevation"),
    target = "Guild.1",
    response_type = "abundance",
    model_kind = "RF_shift_log1p",
    mtry = 3L
  ),
  # Same as your current G2 script.
  g2 = list(
    output_subdir = "g2_results",
    training_file = training_files[["g2"]],
    coord_file = coord_files[["g2"]],
    output_tag = "",
    features = c("ph", "ref_1", "ndvimean"),
    target = "Guild.2",
    response_type = "abundance",
    model_kind = "RF_raw",
    mtry = NA_integer_
  ),
  # Scaled-limited PA alternative. Full-env PA model needs extra global layer paths.
  g3 = list(
    output_subdir = "g3_results",
    training_file = training_files[["g3"]],
    coord_file = coord_files[["g3"]],
    output_tag = "",
    features = c("ref_1", "ndvimean"),
    target = "Guild.3",
    response_type = "presence_absence",
    model_kind = "RF_unweighted",
    mtry = NA_integer_
  ),
  # Extra continuous abundance map for G3. The main G3 map remains PA because
  # abundance spatial CV was not reliable; this output is for visualization/supplement.
  g3_abundance = list(
    output_subdir = "g3_results",
    training_file = training_files[["g3"]],
    coord_file = coord_files[["g3"]],
    output_tag = "_abundance",
    features = c("ndvimean", "ref_1"),
    target = "Guild.3",
    response_type = "abundance",
    model_kind = "RF_hurdle_shift_log1p",
    mtry = 1L
  ),
  g4 = list(
    output_subdir = "g4_results",
    training_file = training_files[["g4"]],
    coord_file = coord_files[["g4"]],
    output_tag = "",
    features = c("ph", "soc", "ref_1", "biomass"),
    target = "Guild.4",
    response_type = "abundance",
    model_kind = "RF_shift_log1p",
    mtry = 4L
  ),
  g5 = list(
    output_subdir = "g5_results",
    training_file = training_files[["g5"]],
    coord_file = coord_files[["g5"]],
    output_tag = "",
    features = c("ref_1", "ndvimean", "biomass", "ph"),
    target = "Guild.5",
    response_type = "abundance",
    model_kind = "RF_hurdle_shift_log1p",
    mtry = 4L
  )
)

if (!guild_id %in% names(guild_configs)) {
  stop("guild_id must be one of: ", paste(names(guild_configs), collapse = ", "))
}
cfg_guild <- guild_configs[[guild_id]]

path_prefix <- file.path(uncertainty_dir, cfg_guild$output_subdir, "pca_multilevel1")
if (!dir.exists(path_prefix)) dir.create(path_prefix, recursive = TRUE)

training_file <- cfg_guild$training_file
coord_file <- cfg_guild$coord_file
features <- cfg_guild$features
target <- cfg_guild$target
response_type <- cfg_guild$response_type
model_kind <- cfg_guild$model_kind
output_tag <- cfg_guild$output_tag
if (is.null(output_tag)) output_tag <- ""
model_mtry <- cfg_guild$mtry
if (is.na(model_mtry)) model_mtry <- max(1, floor(sqrt(length(features))))
global_files <- global_file_catalog[features]
if (any(is.na(global_files))) {
  stop("Missing global CSV path(s) for feature(s): ", paste(features[is.na(global_files)], collapse = ", "))
}

tagged_name <- function(prefix, suffix = ".csv") {
  paste0(prefix, output_tag, suffix)
}

tagged_profile_name <- function(prefix, profile, suffix = ".csv") {
  paste0(prefix, output_tag, "_", profile, suffix)
}

ntree_main <- 2000
predict_threads <- 8

# FALSE is recommended for the first global run.
# TRUE means prediction.csv is the mean of bootstrap models, much slower.
use_bootstrap_prediction <- FALSE
n_boot <- 100
ntree_boot <- 500
train_cores <- 28

# Multiple extrapolation profiles.
# strict   = conservative map: fewer pixels treated as reliable.
# moderate = recommended main map, using the 99% empirical training-distance envelope.
# lenient  = sensitivity map: more pixels retained as reliable.
risk_profiles <- data.table(
  profile = c("strict", "moderate", "lenient"),
  mahalanobis_quantile = c(0.975, 0.990, 0.995),
  geo_cap_quantile = c(0.975, 0.990, 0.995),
  geo_cap_multiplier = c(1.5, 2.0, 2.5),
  geo_cap_min_km = c(1000, 1500, 2500),
  env_weight = c(2, 2, 1),
  geo_weight = c(1, 1, 1),
  combined_risk_cutoff = c(0.60, 0.75, 0.90)
)

# Recommended for extrapolation assessment:
# TRUE  = PCA-Mahalanobis distance. PCA is fit on training environmental space,
#         PCs explaining pca_variance of training variance are retained, and
#         Mahalanobis distance is calculated in retained PC space.
# FALSE = raw-predictor Mahalanobis distance.
use_pca_mahalanobis <- TRUE
pca_variance <- 0.95

# Choose geographic caps from the training-point nearest-neighbour
# distance distribution. This makes the geographic threshold data-driven
# and produces one cap for each strict/moderate/lenient profile.
use_auto_geo_cap <- TRUE

# Assumes rows in each input CSV are north-to-south:
# row 1 near 90N, final row near 90S.
lat_descending <- TRUE
if (!lat_descending) {
  stop("This script assumes north-to-south CSV rows. Flip rows first if needed.")
}

# -----------------------------
# 1. Helper functions
# -----------------------------
append_csv_row <- function(x, file) {
  fwrite(data.table(t(x)), file, append = TRUE, col.names = FALSE, na = "NA")
}

lonlat_to_xyz <- function(lon, lat) {
  lonr <- lon * pi / 180
  latr <- lat * pi / 180
  cbind(cos(latr) * cos(lonr), cos(latr) * sin(lonr), sin(latr))
}

earth_radius_km <- 6371.0088

# -----------------------------
# 2. Load training data
# -----------------------------
cat("Loading training data...\n")
env.meta <- fread(training_file)
if (!"biome" %in% names(env.meta)) env.meta[, biome := factor(1)]

looks_like_real_lonlat <- function(lat, lon) {
  all(is.finite(lat), is.finite(lon)) &&
    min(lat, na.rm = TRUE) >= -90 &&
    max(lat, na.rm = TRUE) <= 90 &&
    min(lon, na.rm = TRUE) >= -180 &&
    max(lon, na.rm = TRUE) <= 180 &&
    diff(range(lat, na.rm = TRUE)) > 10 &&
    diff(range(lon, na.rm = TRUE)) > 10
}

if (!looks_like_real_lonlat(env.meta$Lat, env.meta$Lon)) {
  cat("Lat/Lon in scaled training file look standardized; replacing them with raw coordinates from:\n")
  cat(coord_file, "\n")
  coord.meta <- fread(coord_file)
  if (!all(c("Lat", "Lon") %in% names(coord.meta))) {
    stop("Raw coordinate file does not contain Lat/Lon: ", coord_file)
  }
  if (nrow(coord.meta) != nrow(env.meta)) {
    stop(
      "Cannot safely replace Lat/Lon because row counts differ. scaled n=",
      nrow(env.meta), "; raw n=", nrow(coord.meta), "; file=", coord_file
    )
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

if (response_type == "presence_absence") {
  train_df$.y <- factor(
    ifelse(train_df[[target]] > min(train_df[[target]], na.rm = TRUE) + 1e-10, "present", "absent"),
    levels = c("absent", "present")
  )
} else {
  train_df$.y <- train_df[[target]]
}

cat("Training rows:", nrow(train_df), "\n")
cat("Guild:", guild_id, "\n")
cat("Features:", paste(features, collapse = ", "), "\n")
cat("Response type:", response_type, "\n")
cat("Model kind:", model_kind, "\n")
cat("mtry:", model_mtry, "\n")

# -----------------------------
# 3. Fit ranger model
# -----------------------------
fit_model_object <- function(d, trees, threads = 1, importance = "none", seed = 123) {
  x <- d[, features, drop = FALSE]
  y <- d[[target]]
  mtry <- max(1, min(ncol(x), model_mtry))

  if (response_type == "presence_absence") {
    yy <- factor(
      ifelse(y > min(y, na.rm = TRUE) + 1e-10, "present", "absent"),
      levels = c("absent", "present")
    )
    train <- data.frame(.y = yy, x, check.names = FALSE)
    class_weights <- NULL
    if (identical(model_kind, "RF_balanced")) {
      tab <- table(yy)
      class_weights <- as.numeric(sum(tab) / (length(tab) * tab))
      names(class_weights) <- names(tab)
    }
    fit <- ranger(
      .y ~ .,
      data = train,
      num.trees = trees,
      mtry = mtry,
      min.node.size = 5,
      probability = TRUE,
      class.weights = class_weights,
      importance = importance,
      num.threads = threads,
      seed = seed
    )
    return(list(kind = "presence_absence", fit = fit))
  }

  if (model_kind == "RF_raw") {
    fit <- ranger(
      .y ~ .,
      data = data.frame(.y = y, x, check.names = FALSE),
      num.trees = trees,
      mtry = mtry,
      min.node.size = 5,
      importance = importance,
      num.threads = threads,
      seed = seed
    )
    return(list(kind = "RF_raw", fit = fit))
  }

  baseline <- min(y, na.rm = TRUE)
  shifted <- y - baseline

  if (model_kind == "RF_shift_log1p") {
    fit <- ranger(
      .y ~ .,
      data = data.frame(.y = log1p(shifted), x, check.names = FALSE),
      num.trees = trees,
      mtry = mtry,
      min.node.size = 5,
      importance = importance,
      num.threads = threads,
      seed = seed
    )
    return(list(kind = "RF_shift_log1p", fit = fit, baseline = baseline))
  }

  if (model_kind == "RF_hurdle_shift_log1p") {
    tol <- sqrt(.Machine$double.eps)
    present <- factor(ifelse(y > baseline + tol, "present", "absent"), levels = c("absent", "present"))
    class_fit <- ranger(
      .y ~ .,
      data = data.frame(.y = present, x, check.names = FALSE),
      num.trees = trees,
      mtry = mtry,
      min.node.size = 5,
      probability = TRUE,
      importance = importance,
      num.threads = threads,
      seed = seed
    )
    pos <- y > baseline + tol
    if (sum(pos) < 8) {
      reg_fit <- NULL
      mean_pos_excess <- mean(shifted[pos], na.rm = TRUE)
    } else {
      reg_fit <- ranger(
        .y ~ .,
        data = data.frame(.y = log1p(shifted[pos]), x[pos, , drop = FALSE], check.names = FALSE),
        num.trees = trees,
        mtry = mtry,
        min.node.size = 3,
        importance = importance,
        num.threads = threads,
        seed = seed + 17
      )
      mean_pos_excess <- NA_real_
    }
    return(list(
      kind = "RF_hurdle_shift_log1p",
      class_fit = class_fit,
      reg_fit = reg_fit,
      mean_pos_excess = mean_pos_excess,
      baseline = baseline
    ))
  }

  stop("Unsupported model_kind: ", model_kind)
}

predict_one_model <- function(model, newdata) {
  if (model$kind == "presence_absence") {
    pr <- predict(model$fit, data = newdata, num.threads = predict_threads)$predictions
    return(as.numeric(pr[, "present"]))
  }
  if (model$kind == "RF_raw") {
    return(as.numeric(predict(model$fit, data = newdata, num.threads = predict_threads)$predictions))
  }
  if (model$kind == "RF_shift_log1p") {
    pr <- predict(model$fit, data = newdata, num.threads = predict_threads)$predictions
    return(as.numeric(expm1(pr) + model$baseline))
  }
  if (model$kind == "RF_hurdle_shift_log1p") {
    prob <- predict(model$class_fit, data = newdata, num.threads = predict_threads)$predictions[, "present"]
    if (is.null(model$reg_fit)) {
      pos_excess <- rep(model$mean_pos_excess, nrow(newdata))
    } else {
      pos_excess <- expm1(predict(model$reg_fit, data = newdata, num.threads = predict_threads)$predictions)
    }
    return(as.numeric(model$baseline + prob * pos_excess))
  }
  stop("Unsupported fitted model kind: ", model$kind)
}

extract_importance <- function(model) {
  imp <- NULL
  if (!is.null(model$fit)) imp <- model$fit$variable.importance
  if (is.null(imp) && !is.null(model$class_fit)) imp <- model$class_fit$variable.importance
  if (is.null(imp)) imp <- setNames(rep(NA_real_, length(features)), features)
  imp
}

cat("Training main ranger model...\n")
main_rf <- fit_model_object(train_df, trees = ntree_main, threads = predict_threads, importance = "permutation", seed = 123)

fwrite(
  data.frame(variable = names(extract_importance(main_rf)), importance = as.numeric(extract_importance(main_rf))),
  file.path(path_prefix, paste0(guild_id, "_ranger_importance.csv"))
)

boot_models <- NULL
if (use_bootstrap_prediction) {
  cat("Training", n_boot, "bootstrap ranger models...\n")
  boot_models <- mclapply(seq_len(n_boot), function(b) {
    set.seed(123 + b)
    boot_indices <- unlist(lapply(unique(train_df$biome), function(lv) {
      idx <- which(train_df$biome == lv)
      sample(idx, size = length(idx), replace = TRUE)
    }))
    fit_model_object(train_df[boot_indices, , drop = FALSE], trees = ntree_boot, threads = 1, importance = "none", seed = 123 + b)
  }, mc.cores = train_cores)
} else {
  cat("Bootstrap prediction disabled; prediction.csv will use the main model.\n")
}

# -----------------------------
# 4. Environmental extrapolation parameters
# -----------------------------
cat("Preparing environmental extrapolation parameters...\n")
X_train <- as.matrix(train_df[, features, drop = FALSE])
X_train_df <- as.data.frame(X_train)

make_mahalanobis_fun <- function(scores_train) {
  scores_train <- as.matrix(scores_train)
  center <- colMeans(scores_train, na.rm = TRUE)
  cov_mat <- cov(scores_train, use = "complete.obs")
  if (is.null(dim(cov_mat))) cov_mat <- matrix(cov_mat, nrow = 1)
  ridge <- max(diag(cov_mat), na.rm = TRUE) * 1e-6
  cov_inv <- solve(cov_mat + diag(ridge, ncol(cov_mat)))
  function(scores_new) {
    scores_new <- as.matrix(scores_new)
    z <- sweep(scores_new, 2, center, "-")
    rowSums((z %*% cov_inv) * z)
  }
}

if (use_pca_mahalanobis) {
  # Variables are assumed already scaled consistently with training data.
  # Therefore PCA centers but does not re-scale by default.
  pca_fit <- prcomp(X_train_df, center = TRUE, scale. = FALSE)
  variance_explained <- cumsum(pca_fit$sdev^2) / sum(pca_fit$sdev^2)
  n_pc <- which(variance_explained >= pca_variance)[1]
  pc_keep <- seq_len(n_pc)
  train_scores <- predict(pca_fit, newdata = X_train_df)[, pc_keep, drop = FALSE]
  mahalanobis_in_space <- make_mahalanobis_fun(train_scores)
  environmental_distance <- function(x) {
    scores <- predict(pca_fit, newdata = as.data.frame(x))[, pc_keep, drop = FALSE]
    mahalanobis_in_space(scores)
  }
  extrapolation_method <- "PCA_Mahalanobis"
  extrapolation_df <- n_pc
  cat("Using PCA-Mahalanobis with", n_pc, "PC axes; cumulative variance =",
      round(variance_explained[n_pc], 4), "\n")
} else {
  mahalanobis_in_space <- make_mahalanobis_fun(X_train)
  environmental_distance <- function(x) mahalanobis_in_space(x)
  extrapolation_method <- "Raw_predictor_Mahalanobis"
  extrapolation_df <- length(features)
  cat("Using raw-predictor Mahalanobis with", length(features), "variables.\n")
}

md_train <- environmental_distance(X_train_df)
risk_profiles[, mahalanobis_threshold := as.numeric(
  quantile(md_train, mahalanobis_quantile, na.rm = TRUE)
), by = profile]

fwrite(
  rbind(
    data.frame(
      profile = risk_profiles$profile,
      threshold_type = "training_empirical",
      probability = risk_profiles$mahalanobis_quantile,
      method = extrapolation_method,
      df = extrapolation_df,
      n_original_features = length(features),
      pca_variance_target = ifelse(use_pca_mahalanobis, pca_variance, NA_real_),
      threshold = risk_profiles$mahalanobis_threshold
    ),
    data.frame(
      profile = "reference",
      threshold_type = "chi_square",
      probability = c(0.95, 0.975, 0.99),
      method = extrapolation_method,
      df = extrapolation_df,
      n_original_features = length(features),
      pca_variance_target = ifelse(use_pca_mahalanobis, pca_variance, NA_real_),
      threshold = qchisq(c(0.95, 0.975, 0.99), df = extrapolation_df)
    )
  ),
  file.path(path_prefix, paste0(guild_id, "_mahalanobis_thresholds_multilevel.csv"))
)
cat("Mahalanobis empirical thresholds:\n")
print(risk_profiles[, .(profile, mahalanobis_quantile, mahalanobis_threshold)])

# -----------------------------
# 5. Nearest-sample distance setup
# -----------------------------
# Use unique coordinates so duplicated samples at the same location do not
# create zero nearest-neighbour distances and artificially shrink geo_cap_km.
train_coords <- unique(train_df[, c("Lon", "Lat")])
train_xyz <- lonlat_to_xyz(train_coords$Lon, train_coords$Lat)

if (nrow(train_coords) >= 2) {
  train_nn <- RANN::nn2(data = train_xyz, query = train_xyz, k = 2)
  train_nn_chord <- train_nn$nn.dists[, 2]
  train_nn_km <- 2 * earth_radius_km * asin(pmin(train_nn_chord / 2, 1))
  train_nn_probs <- c(0.5, 0.75, 0.9, 0.95, 0.975, 0.99)
  train_nn_summary <- data.frame(
    probability = train_nn_probs,
    nearest_neighbor_distance_km = as.numeric(quantile(train_nn_km, train_nn_probs, na.rm = TRUE)),
    n_unique_training_locations = nrow(train_coords)
  )
  fwrite(
    train_nn_summary,
    file.path(path_prefix, paste0(guild_id, "_training_nearest_neighbor_distance_summary.csv"))
  )

  risk_profiles[, geo_cap_km := pmax(
    geo_cap_min_km,
    as.numeric(quantile(train_nn_km, geo_cap_quantile, na.rm = TRUE)) * geo_cap_multiplier
  ), by = profile]

  fwrite(
    risk_profiles[, .(
      profile,
      use_auto_geo_cap,
      geo_cap_km,
      geo_cap_quantile,
      geo_cap_multiplier,
      geo_cap_min_km,
      mahalanobis_quantile,
      mahalanobis_threshold,
      env_weight,
      geo_weight,
      combined_risk_cutoff,
      n_unique_training_locations = nrow(train_coords)
    )],
    file.path(path_prefix, paste0(guild_id, "_extrapolation_profile_settings.csv"))
  )
  cat("Geographic caps used:\n")
  print(risk_profiles[, .(profile, geo_cap_quantile, geo_cap_multiplier, geo_cap_km)])
} else {
  stop("Fewer than two unique training coordinates; cannot calculate nearest-sample distance profiles.")
}

nearest_dist_km <- function(lon, lat) {
  qxyz <- lonlat_to_xyz(lon, lat)
  nn <- RANN::nn2(data = train_xyz, query = qxyz, k = 1)
  chord <- nn$nn.dists[, 1]
  2 * earth_radius_km * asin(pmin(chord / 2, 1))
}

# -----------------------------
# 6. Load global environmental CSV matrices
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
cat("Global grid:", n_lon, "columns x", n_lat, "rows =", format(n_lon * n_lat, scientific = FALSE), "pixels\n")

lon_vec <- -180 + (seq_len(n_lon) - 0.5) * 360 / n_lon
lat_vec <- 90 - (seq_len(n_lat) - 0.5) * 180 / n_lat

fwrite(
  data.frame(
    n_lon = n_lon,
    n_lat = n_lat,
    xmin = -180,
    xmax = 180,
    ymin = -90,
    ymax = 90,
    crs = "EPSG:4326",
    orientation = "rows_are_longitude_columns__columns_are_latitude_rows",
    row1_lon = lon_vec[1],
    last_row_lon = lon_vec[n_lon],
    col1_lat = lat_vec[1],
    last_col_lat = lat_vec[n_lat]
  ),
  file.path(path_prefix, paste0(guild_id, "_csv_grid_metadata.csv"))
)

# -----------------------------
# 7. Prepare output CSV files
# -----------------------------
file_prediction <- file.path(path_prefix, tagged_name("prediction"))
file_md <- file.path(path_prefix, tagged_name("mahalanobis_distance"))
file_nearest <- file.path(path_prefix, tagged_name("nearest_sample_distance"))

profile_files <- lapply(risk_profiles$profile, function(p) {
  list(
    env_mask = file.path(path_prefix, tagged_profile_name("environmental_extrapolation_mask", p)),
    risk = file.path(path_prefix, tagged_profile_name("combined_extrapolation_risk", p)),
    reliable = file.path(path_prefix, tagged_profile_name("reliable_prediction_mask", p)),
    prediction_reliable = file.path(path_prefix, paste0("prediction", output_tag, "_reliable_", p, ".csv"))
  )
})
names(profile_files) <- risk_profiles$profile

out_files <- c(
  file_prediction,
  file_md,
  file_nearest,
  unlist(profile_files, use.names = FALSE)
)
invisible(file.remove(out_files[file.exists(out_files)]))

profile_summary <- data.table(
  profile = risk_profiles$profile,
  valid_pixels = 0,
  environmental_ood_pixels = 0,
  reliable_pixels = 0,
  combined_risk_sum = 0
)

# -----------------------------
# 8. Global prediction loop
# -----------------------------
cat("Starting global prediction and CSV writing...\n")
start_time <- Sys.time()

for (i in seq_len(n_lon)) {
  env.curr <- as.data.frame(lapply(features, function(v) global_layers[[v]][[i]]))
  names(env.curr) <- features

  valid <- complete.cases(env.curr)

  prediction <- rep(NA_real_, n_lat)
  m_dist <- rep(NA_real_, n_lat)
  nearest_km <- rep(NA_real_, n_lat)
  profile_rows <- lapply(profile_files, function(x) {
    list(
      env_mask = rep(NA_integer_, n_lat),              # 1 = extrapolated/outside; 0 = represented; NA = no data
      risk = rep(NA_real_, n_lat),                     # 0 low risk; 1 high risk; NA = no data
      reliable = rep(NA_integer_, n_lat),              # 1 reliable; 0 unreliable; NA = no data
      prediction_reliable = rep(NA_real_, n_lat)       # prediction, but unreliable pixels are NA
    )
  })

  if (any(valid)) {
    x_valid <- env.curr[valid, , drop = FALSE]

    # Environmental extrapolation: PCA-Mahalanobis or raw Mahalanobis distance.
    m_dist[valid] <- environmental_distance(x_valid)

    # Geographic extrapolation: nearest training sample distance.
    lon_col <- rep(lon_vec[i], n_lat)
    nearest_km[valid] <- nearest_dist_km(lon_col[valid], lat_vec[valid])

    # Prediction.
    if (use_bootstrap_prediction) {
      all_preds <- sapply(boot_models, predict_one_model, newdata = x_valid)
      if (is.null(dim(all_preds))) all_preds <- matrix(all_preds, ncol = 1)
      prediction[valid] <- rowMeans(all_preds, na.rm = TRUE)
    } else {
      prediction[valid] <- predict_one_model(main_rf, x_valid)
    }

    # Profile-specific extrapolation masks and reliable predictions.
    for (p in risk_profiles$profile) {
      prof <- risk_profiles[profile == p]
      env_mask_valid <- as.integer(m_dist[valid] > prof$mahalanobis_threshold)
      env_risk_valid <- pmin(m_dist[valid] / prof$mahalanobis_threshold, 1)
      geo_risk_valid <- pmin(nearest_km[valid] / prof$geo_cap_km, 1)
      risk_valid <- (prof$env_weight * env_risk_valid + prof$geo_weight * geo_risk_valid) /
        (prof$env_weight + prof$geo_weight)
      reliable_valid <- as.integer(
        env_mask_valid == 0 &
          nearest_km[valid] <= prof$geo_cap_km &
          risk_valid <= prof$combined_risk_cutoff
      )

      profile_rows[[p]]$env_mask[valid] <- env_mask_valid
      profile_rows[[p]]$risk[valid] <- risk_valid
      profile_rows[[p]]$reliable[valid] <- reliable_valid
      profile_rows[[p]]$prediction_reliable[valid] <- ifelse(reliable_valid == 1, prediction[valid], NA_real_)

      profile_summary[profile == p, `:=`(
        valid_pixels = valid_pixels + sum(valid),
        environmental_ood_pixels = environmental_ood_pixels + sum(env_mask_valid == 1, na.rm = TRUE),
        reliable_pixels = reliable_pixels + sum(reliable_valid == 1, na.rm = TRUE),
        combined_risk_sum = combined_risk_sum + sum(risk_valid, na.rm = TRUE)
      )]
    }
  }

  append_csv_row(prediction, file_prediction)
  append_csv_row(m_dist, file_md)
  append_csv_row(nearest_km, file_nearest)
  for (p in names(profile_files)) {
    append_csv_row(profile_rows[[p]]$env_mask, profile_files[[p]]$env_mask)
    append_csv_row(profile_rows[[p]]$risk, profile_files[[p]]$risk)
    append_csv_row(profile_rows[[p]]$reliable, profile_files[[p]]$reliable)
    append_csv_row(profile_rows[[p]]$prediction_reliable, profile_files[[p]]$prediction_reliable)
  }

  if (i %% 100 == 0 || i == n_lon) {
    elapsed_h <- as.numeric(difftime(Sys.time(), start_time, units = "hours"))
    done <- i / n_lon
    eta_h <- elapsed_h / max(done, 1e-9) - elapsed_h
    cat(
      "Progress:", round(done * 100, 2), "%",
      "| longitude columns:", i, "/", n_lon,
      "| elapsed h:", round(elapsed_h, 2),
      "| ETA h:", round(eta_h, 2), "\n"
    )
  }
}

cat("Done. Final CSV files:\n")
cat(file_prediction, "\n")
cat(file_md, "\n")
cat(file_nearest, "\n")
for (p in names(profile_files)) {
  cat(profile_files[[p]]$env_mask, "\n")
  cat(profile_files[[p]]$risk, "\n")
  cat(profile_files[[p]]$reliable, "\n")
  cat(profile_files[[p]]$prediction_reliable, "\n")
}

profile_summary <- merge(
  profile_summary,
  risk_profiles,
  by = "profile",
  all.x = TRUE,
  sort = FALSE
)
profile_summary[, environmental_ood_fraction := environmental_ood_pixels / valid_pixels]
profile_summary[, reliable_fraction := reliable_pixels / valid_pixels]
profile_summary[, mean_combined_risk := combined_risk_sum / valid_pixels]
summary_file <- file.path(path_prefix, paste0(guild_id, "_extrapolation_profile_pixel_summary.csv"))
fwrite(profile_summary, summary_file)
cat(summary_file, "\n")
