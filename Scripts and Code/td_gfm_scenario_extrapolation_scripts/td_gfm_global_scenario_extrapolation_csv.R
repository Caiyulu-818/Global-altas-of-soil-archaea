#!/usr/bin/env Rscript

# ============================================================
# TD / GFM global scenario prediction + extrapolation risk
# CSV-only version for server runs: no terra, no raster, no gdal
# ============================================================
#
# Packages required:
#   ranger
#   data.table
#   RANN
#
# Usage examples:
#   Rscript td_gfm_global_scenario_extrapolation_csv.R TD
#   Rscript td_gfm_global_scenario_extrapolation_csv.R GFM
#   Rscript td_gfm_global_scenario_extrapolation_csv.R all
#
# Outputs:
#   /public/home/ylwang/Gisdata/TC_New/lucy/cropland/scenario_extrapolation/<TD|GFM>/<scenario>/
#     prediction.csv
#     mahalanobis_distance.csv
#     nearest_sample_distance.csv
#     environmental_extrapolation_mask_<strict|moderate|lenient>.csv
#     combined_extrapolation_risk_<strict|moderate|lenient>.csv
#     reliable_prediction_mask_<strict|moderate|lenient>.csv
#     prediction_reliable_<strict|moderate|lenient>.csv
#
# Meaning:
#   environmental_extrapolation_mask: 1 = environmental OOD, 0 = represented
#   reliable_prediction_mask: 1 = reliable, 0 = mask out
#   prediction_reliable: prediction already masked by reliable_prediction_mask

suppressPackageStartupMessages({
  library(data.table)
  library(ranger)
})

if (!requireNamespace("RANN", quietly = TRUE)) {
  stop("Package RANN is required. Install on the server with install.packages('RANN').")
}

set.seed(20260624)

# -----------------------------
# 0. User settings
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
run_target <- if (exists("run_target_override", inherits = TRUE)) {
  toupper(get("run_target_override", inherits = TRUE))
} else if (length(args) >= 1) {
  toupper(args[[1]])
} else {
  "ALL"
}
if (!run_target %in% c("TD", "GFM", "ALL")) {
  stop("First argument must be TD, GFM, or all.")
}

base_dir <- "/public/home/ylwang/Gisdata/TC_New/lucy"
training_dir <- file.path(base_dir, "cropland")
env_csv_root <- "/public/home/ZhangFZ/Gisdata/resample/csv"
out_root <- file.path(training_dir, "scenario_extrapolation")

# If TRUE, cropland_scenario is included as a model predictor so SSP layers
# actually change the predictions. This is the setting for future scenario maps.
# If FALSE, the exact screened environmental-only models are used; then SSP
# scenarios do not change prediction, only an optional OOD diagnosis would differ.
use_cropland_as_predictor <- TRUE

ntree_main <- 2000
predict_threads <- 8

# Multi-level extrapolation settings. Strict is safest for main maps;
# moderate/lenient are useful sensitivity levels.
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

use_pca_mahalanobis <- TRUE
pca_variance <- 0.95

# CSV orientation follows your original global scripts:
# global layer data.tables are n_lat rows x n_lon columns;
# outputs are written as n_lon rows x n_lat columns.
lat_descending <- TRUE
if (!lat_descending) stop("This script assumes north-to-south CSV rows.")

# -----------------------------
# 1. Model settings from previous model testing
# -----------------------------
# TD: best validated model type was RF_raw. Core screened predictors:
#   New_NDVImean; New_wc2.1_30s_bio_12; aridity;
#   Soil.pH.x.10.in.H2O; New_wc2.1_30s_bio_1
#
# GFM: best validated model type was RF_log1p. Core screened predictors:
#   aridity; New_NDVImean; New_northness
#
# We add cropland_scenario when use_cropland_as_predictor = TRUE so the
# four SSP layers can drive different future predictions.
dataset_configs <- list(
  TD = list(
    training_file = file.path(training_dir, "TD.env.meta.csv"),
    target = "Shannon.TD",
    model_kind = "RF_raw",
    features_static = c(
      "New_NDVImean",
      "New_wc2.1_30s_bio_12",
      "aridity",
      "Soil.pH.x.10.in.H2O",
      "New_wc2.1_30s_bio_1"
    ),
    training_cropland_col = "cropland2020_SSP245"
  ),
  GFM = list(
    training_file = file.path(training_dir, "GFM.env.meta.with_lonlat.csv"),
    target = "GFM",
    model_kind = "RF_log1p",
    features_static = c(
      "aridity",
      "New_NDVImean",
      "New_northness"
    ),
    training_cropland_col = "cropland2020_SSP245"
  )
)

# Global environmental layer catalog.
global_file_catalog <- c(
  aridity = file.path(env_csv_root, "New_ai_et0.csv"),
  Cation.exchange.capacity.of.soil = file.path(env_csv_root, "New_Average_cec2.csv"),
  Coarse.fragments.volumetric = file.path(env_csv_root, "New_Average_cfvo2.csv"),
  Soil.pH.x.10.in.H2O = file.path(env_csv_root, "New_Average_phh2o2.csv"),
  Sand.content = file.path(env_csv_root, "New_Average_sand2.csv"),
  Silt.content = file.path(env_csv_root, "New_Average_silt2.csv"),
  Soil.organic.carbon.density = file.path(env_csv_root, "New_Average_ocd2.csv"),
  New_eastness = file.path(env_csv_root, "New_eastness.csv"),
  New_northness = file.path(env_csv_root, "New_northness.csv"),
  New_NDVImean = file.path(env_csv_root, "New_NDVImean.csv"),
  New_Nppmean = file.path(env_csv_root, "New_Nppmean.csv"),
  New_wc2.1_30s_bio_1 = file.path(env_csv_root, "New_wc2.1_30s_bio_1.csv"),
  New_wc2.1_30s_bio_12 = file.path(env_csv_root, "New_wc2.1_30s_bio_12.csv"),
  New_Mean_of_annual_predicted_ECe_1980_2018 = file.path(env_csv_root, "New_Mean_of_annual_predicted_ECe_1980_2018.csv"),
  New_slope = file.path(env_csv_root, "New_slope.csv"),
  New_soilbiomass_1km = file.path(env_csv_root, "New_soilbiomass_1km.csv")
)

cropland_root <- file.path(training_dir)
cropland_scenarios <- data.table(
  scenario = c(
    "cropland2100_SSP126",
    "cropland2100_SSP370",
    "cropland2020_SSP245",
    "cropland2100_SSP245"
  ),
  file = c(
    file.path(cropland_root, "New_globalCropland_2100CE_SSP126_LandOnly.csv", "New_globalCropland_2100CE_SSP126_LandOnly.csv"),
    file.path(cropland_root, "New_globalCropland_2100CE_SSP370_LandOnly.csv", "New_globalCropland_2100CE_SSP370_LandOnly.csv"),
    file.path(cropland_root, "New_globalCropland_2020CE_SSP245_LandOnly.csv", "New_globalCropland_2020CE_SSP245_LandOnly.csv"),
    file.path(cropland_root, "New_globalCropland_2100CE_SSP245_LandOnly.csv", "New_globalCropland_2100CE_SSP245_LandOnly.csv")
  )
)

# -----------------------------
# 2. Helper functions
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

looks_like_real_lonlat <- function(lat, lon) {
  ok <- is.finite(lat) & is.finite(lon)
  if (sum(ok) < 2) return(FALSE)
  min(lat[ok], na.rm = TRUE) >= -90 &&
    max(lat[ok], na.rm = TRUE) <= 90 &&
    min(lon[ok], na.rm = TRUE) >= -180 &&
    max(lon[ok], na.rm = TRUE) <= 180 &&
    diff(range(lat[ok], na.rm = TRUE)) > 10 &&
    diff(range(lon[ok], na.rm = TRUE)) > 10
}

safe_r2 <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  if (length(obs) < 2 || sum((obs - mean(obs))^2) == 0) return(NA_real_)
  1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
}

metrics_reg <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  data.frame(
    R2 = safe_r2(obs, pred),
    RMSE = sqrt(mean((obs - pred)^2)),
    MAE = mean(abs(obs - pred)),
    n = length(obs)
  )
}

make_mahalanobis_fun <- function(train_matrix) {
  train_matrix <- as.matrix(train_matrix)
  center <- colMeans(train_matrix, na.rm = TRUE)
  cov_mat <- cov(train_matrix, use = "complete.obs")
  ridge <- max(diag(cov_mat), na.rm = TRUE) * 1e-6
  cov_inv <- solve(cov_mat + diag(ridge, ncol(cov_mat)))
  function(x) {
    x <- as.matrix(x)
    z <- sweep(x, 2, center, "-")
    rowSums((z %*% cov_inv) * z)
  }
}

fit_model <- function(train_df, target, features, model_kind) {
  x <- train_df[, features, drop = FALSE]
  y <- train_df[[target]]
  mtry <- max(1, floor(sqrt(length(features))))

  if (model_kind == "RF_raw") {
    fit <- ranger(
      .y ~ .,
      data = data.frame(.y = y, x, check.names = FALSE),
      num.trees = ntree_main,
      mtry = mtry,
      min.node.size = 5,
      importance = "permutation",
      num.threads = predict_threads,
      seed = 20260624
    )
    return(list(kind = model_kind, fit = fit))
  }

  if (model_kind == "RF_log1p") {
    baseline <- min(y, na.rm = TRUE)
    shifted <- y - baseline
    fit <- ranger(
      .y ~ .,
      data = data.frame(.y = log1p(shifted), x, check.names = FALSE),
      num.trees = ntree_main,
      mtry = mtry,
      min.node.size = 5,
      importance = "permutation",
      num.threads = predict_threads,
      seed = 20260624
    )
    return(list(kind = model_kind, fit = fit, baseline = baseline))
  }

  stop("Unsupported model_kind: ", model_kind)
}

predict_model <- function(model, newdata) {
  p <- predict(model$fit, data = newdata, num.threads = predict_threads)$predictions
  if (model$kind == "RF_log1p") {
    p <- expm1(p) + model$baseline
    p <- pmax(p, model$baseline)
  }
  as.numeric(p)
}

remove_old_outputs <- function(files) {
  old <- files[file.exists(files)]
  if (length(old) > 0) invisible(file.remove(old))
}

load_global_layers <- function(features) {
  missing <- setdiff(features, names(global_file_catalog))
  if (length(missing) > 0) {
    stop("No global CSV path defined for feature(s): ", paste(missing, collapse = ", "))
  }
  layers <- lapply(features, function(v) {
    path <- global_file_catalog[[v]]
    if (!file.exists(path)) stop("Missing global layer: ", path)
    cat("  reading", v, "from", path, "\n")
    fread(path)
  })
  names(layers) <- features
  dims <- vapply(layers, function(x) paste(dim(x), collapse = "x"), character(1))
  if (length(unique(dims)) != 1) {
    stop("Global layer dimensions do not match: ", paste(names(dims), dims, collapse = "; "))
  }
  layers
}

# -----------------------------
# 3. Dataset run
# -----------------------------
run_dataset <- function(dataset_name) {
  cfg <- dataset_configs[[dataset_name]]
  if (is.null(cfg)) stop("Unknown dataset: ", dataset_name)

  cat("\n==============================\n")
  cat("Dataset:", dataset_name, "\n")
  cat("==============================\n")

  if (!file.exists(cfg$training_file)) {
    stop("Missing training file: ", cfg$training_file)
  }

  train0 <- fread(cfg$training_file)
  if (!all(c("Lat", "Lon") %in% names(train0))) {
    stop("Training file must contain real Lat and Lon columns: ", cfg$training_file)
  }
  if (!looks_like_real_lonlat(train0$Lat, train0$Lon)) {
    stop("Lat/Lon do not look like real coordinates in: ", cfg$training_file)
  }
  missing_coords <- sum(!is.finite(train0$Lat) | !is.finite(train0$Lon))
  if (missing_coords > 0) {
    cat("Rows with missing Lat/Lon will be removed before model fitting:", missing_coords, "\n")
  }

  features <- cfg$features_static
  if (use_cropland_as_predictor) features <- c(features, "cropland_scenario")

  needed <- c(cfg$features_static, cfg$target, "Lat", "Lon")
  if (use_cropland_as_predictor) needed <- c(needed, cfg$training_cropland_col)
  missing <- setdiff(needed, names(train0))
  if (length(missing) > 0) {
    stop("Missing columns in training file: ", paste(missing, collapse = ", "))
  }

  train0[, cropland_scenario := get(cfg$training_cropland_col)]
  train_df <- na.omit(as.data.frame(train0[, c(features, cfg$target, "Lat", "Lon"), with = FALSE]))
  cat("Training rows after NA removal:", nrow(train_df), "\n")
  cat("Model kind:", cfg$model_kind, "\n")
  cat("Features:", paste(features, collapse = ", "), "\n")

  model <- fit_model(train_df, cfg$target, features, cfg$model_kind)
  train_pred <- predict_model(model, train_df[, features, drop = FALSE])
  train_metrics <- metrics_reg(train_df[[cfg$target]], train_pred)

  dataset_out <- file.path(out_root, dataset_name)
  dir.create(dataset_out, recursive = TRUE, showWarnings = FALSE)

  importance <- data.frame(
    variable = names(model$fit$variable.importance),
    importance = as.numeric(model$fit$variable.importance),
    stringsAsFactors = FALSE
  )
  fwrite(importance, file.path(dataset_out, paste0(dataset_name, "_ranger_importance.csv")))
  fwrite(train_metrics, file.path(dataset_out, paste0(dataset_name, "_training_fit_metrics.csv")))

  # Environmental extrapolation is calculated in model predictor space,
  # including cropland_scenario when use_cropland_as_predictor = TRUE.
  X_train <- train_df[, features, drop = FALSE]
  if (use_pca_mahalanobis) {
    pca <- prcomp(X_train, center = TRUE, scale. = TRUE)
    var_exp <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
    n_pc <- which(var_exp >= pca_variance)[1]
    train_pca <- pca$x[, seq_len(n_pc), drop = FALSE]
    md_fun <- make_mahalanobis_fun(train_pca)
    environmental_distance <- function(x) {
      px <- predict(pca, newdata = x)[, seq_len(n_pc), drop = FALSE]
      md_fun(px)
    }
    extrapolation_method <- "PCA_Mahalanobis"
    extrapolation_df <- n_pc
    cat("Using PCA-Mahalanobis with", n_pc, "PC axes; cumulative variance =", round(var_exp[n_pc], 4), "\n")
  } else {
    environmental_distance <- make_mahalanobis_fun(X_train)
    extrapolation_method <- "Raw_predictor_Mahalanobis"
    extrapolation_df <- ncol(X_train)
  }

  md_train <- environmental_distance(X_train)
  risk_profiles[, mahalanobis_threshold := as.numeric(
    quantile(md_train, mahalanobis_quantile, na.rm = TRUE)
  ), by = profile]

  train_coords <- unique(train_df[, c("Lon", "Lat")])
  train_xyz <- lonlat_to_xyz(train_coords$Lon, train_coords$Lat)
  if (nrow(train_coords) < 2) stop("Need at least two unique training coordinates.")

  train_nn <- RANN::nn2(data = train_xyz, query = train_xyz, k = 2)
  train_nn_km <- 2 * earth_radius_km * asin(pmin(train_nn$nn.dists[, 2] / 2, 1))
  risk_profiles[, geo_cap_km := pmax(
    geo_cap_min_km,
    as.numeric(quantile(train_nn_km, geo_cap_quantile, na.rm = TRUE)) * geo_cap_multiplier
  ), by = profile]

  nearest_dist_km <- function(lon, lat) {
    qxyz <- lonlat_to_xyz(lon, lat)
    nn <- RANN::nn2(data = train_xyz, query = qxyz, k = 1)
    2 * earth_radius_km * asin(pmin(nn$nn.dists[, 1] / 2, 1))
  }

  fwrite(
    data.frame(
      dataset = dataset_name,
      use_cropland_as_predictor = use_cropland_as_predictor,
      extrapolation_method = extrapolation_method,
      extrapolation_df = extrapolation_df,
      n_training_rows = nrow(train_df),
      n_unique_training_locations = nrow(train_coords),
      feature = features
    ),
    file.path(dataset_out, paste0(dataset_name, "_model_projection_settings.csv"))
  )
  fwrite(
    risk_profiles[, .(
      profile, mahalanobis_quantile, mahalanobis_threshold,
      geo_cap_quantile, geo_cap_multiplier, geo_cap_min_km, geo_cap_km,
      env_weight, geo_weight, combined_risk_cutoff
    )],
    file.path(dataset_out, paste0(dataset_name, "_extrapolation_profile_settings.csv"))
  )

  cat("Loading static global layers...\n")
  static_layers <- load_global_layers(cfg$features_static)
  n_lat <- nrow(static_layers[[1]])
  n_lon <- ncol(static_layers[[1]])
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
    file.path(dataset_out, paste0(dataset_name, "_csv_grid_metadata.csv"))
  )

  for (s in seq_len(nrow(cropland_scenarios))) {
    scen <- cropland_scenarios$scenario[s]
    cropland_file <- cropland_scenarios$file[s]
    if (!file.exists(cropland_file)) {
      stop("Missing cropland scenario layer: ", cropland_file)
    }

    scenario_out <- file.path(dataset_out, scen)
    dir.create(scenario_out, recursive = TRUE, showWarnings = FALSE)
    cat("\n--- Scenario:", scen, "---\n")
    cat("Cropland layer:", cropland_file, "\n")
    cropland_layer <- fread(cropland_file)

    if (!identical(dim(cropland_layer), dim(static_layers[[1]]))) {
      stop("Cropland layer dimensions do not match static layers for scenario: ", scen)
    }

    file_prediction <- file.path(scenario_out, "prediction.csv")
    file_md <- file.path(scenario_out, "mahalanobis_distance.csv")
    file_nearest <- file.path(scenario_out, "nearest_sample_distance.csv")

    profile_files <- lapply(risk_profiles$profile, function(p) {
      list(
        env_mask = file.path(scenario_out, paste0("environmental_extrapolation_mask_", p, ".csv")),
        risk = file.path(scenario_out, paste0("combined_extrapolation_risk_", p, ".csv")),
        reliable = file.path(scenario_out, paste0("reliable_prediction_mask_", p, ".csv")),
        prediction_reliable = file.path(scenario_out, paste0("prediction_reliable_", p, ".csv"))
      )
    })
    names(profile_files) <- risk_profiles$profile

    out_files <- c(file_prediction, file_md, file_nearest, unlist(profile_files, use.names = FALSE))
    remove_old_outputs(out_files)

    profile_summary <- data.table(
      profile = risk_profiles$profile,
      valid_pixels = 0,
      environmental_ood_pixels = 0,
      reliable_pixels = 0,
      combined_risk_sum = 0
    )

    start_time <- Sys.time()
    for (i in seq_len(n_lon)) {
      env.curr <- as.data.frame(lapply(cfg$features_static, function(v) static_layers[[v]][[i]]))
      names(env.curr) <- cfg$features_static
      if (use_cropland_as_predictor) {
        env.curr$cropland_scenario <- cropland_layer[[i]]
      }

      valid <- complete.cases(env.curr)
      prediction <- rep(NA_real_, n_lat)
      m_dist <- rep(NA_real_, n_lat)
      nearest_km <- rep(NA_real_, n_lat)

      profile_rows <- lapply(profile_files, function(x) {
        list(
          env_mask = rep(NA_integer_, n_lat),
          risk = rep(NA_real_, n_lat),
          reliable = rep(NA_integer_, n_lat),
          prediction_reliable = rep(NA_real_, n_lat)
        )
      })

      if (any(valid)) {
        x_valid <- env.curr[valid, features, drop = FALSE]
        m_dist[valid] <- environmental_distance(x_valid)

        lon_col <- rep(lon_vec[i], n_lat)
        nearest_km[valid] <- nearest_dist_km(lon_col[valid], lat_vec[valid])
        prediction[valid] <- predict_model(model, x_valid)

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
          dataset_name, scen,
          "| progress:", round(done * 100, 2), "%",
          "| longitude columns:", i, "/", n_lon,
          "| elapsed h:", round(elapsed_h, 2),
          "| ETA h:", round(eta_h, 2), "\n"
        )
      }
    }

    profile_summary <- merge(profile_summary, risk_profiles, by = "profile", all.x = TRUE, sort = FALSE)
    profile_summary[, environmental_ood_fraction := environmental_ood_pixels / valid_pixels]
    profile_summary[, reliable_fraction := reliable_pixels / valid_pixels]
    profile_summary[, mean_combined_risk := combined_risk_sum / valid_pixels]
    fwrite(profile_summary, file.path(scenario_out, paste0(dataset_name, "_", scen, "_pixel_summary.csv")))

    cat("Scenario done:", scenario_out, "\n")
  }
}

targets <- if (run_target == "ALL") c("TD", "GFM") else run_target
for (x in targets) run_dataset(x)

cat("\nAll requested TD/GFM scenario extrapolation runs finished.\n")
