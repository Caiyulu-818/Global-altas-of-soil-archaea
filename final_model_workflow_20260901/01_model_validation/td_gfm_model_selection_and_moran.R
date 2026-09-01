#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ranger)
  library(caret)
  library(spdep)
  library(CAST)
})

set.seed(20260618)

td_file <- "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/TD.env.meta.csv"
gfm_file <- "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/GFM.env.meta.csv"
out_dir <- "model_outputs_td_gfm"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cfg <- list(
  seed = 20260618,
  candidate_trees = 300,
  final_trees = 1000,
  k_folds = 10,
  mem_count = 5,
  moran_k = 8,
  correlation_cutoff = 0.90,
  aggregate_same_coords = TRUE,
  allow_future_cropland = FALSE
)

mode_value <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
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
  rmse <- sqrt(mean((obs - pred)^2))
  obs_sd <- sd(obs)
  obs_range <- diff(range(obs))
  data.frame(
    R2 = safe_r2(obs, pred),
    RMSE = rmse,
    NRMSE_SD = ifelse(is.finite(obs_sd) && obs_sd > 0, rmse / obs_sd, NA_real_),
    NRMSE_range = ifelse(is.finite(obs_range) && obs_range > 0, rmse / obs_range, NA_real_),
    MAE = mean(abs(obs - pred)),
    n = length(obs),
    stringsAsFactors = FALSE
  )
}

make_split_list <- function(test_list, n, prefix) {
  lapply(seq_along(test_list), function(i) {
    test <- sort(unique(test_list[[i]]))
    list(name = paste0(prefix, i), train = setdiff(seq_len(n), test), test = test)
  })
}

random_folds <- function(y, k = 10, repeats = 1) {
  folds <- createMultiFolds(y, k = k, times = repeats)
  lapply(names(folds), function(nm) {
    train <- sort(unique(folds[[nm]]))
    list(name = nm, train = train, test = setdiff(seq_along(y), train))
  })
}

stratified_biome_folds <- function(y, biome, k = 10) {
  if (length(unique(biome)) > 1) {
    strata <- biome
  } else {
    strata <- cut(y, breaks = unique(quantile(y, probs = seq(0, 1, length.out = min(k, 5) + 1), na.rm = TRUE)), include.lowest = TRUE)
  }
  test_list <- createFolds(strata, k = k, list = TRUE, returnTrain = FALSE)
  make_split_list(test_list, length(y), "Biome")
}

loocv_folds <- function(n) {
  lapply(seq_len(n), function(i) list(name = paste0("LOO", i), train = setdiff(seq_len(n), i), test = i))
}

spatial_kmeans <- function(coords, k = 10) {
  set.seed(cfg$seed)
  km <- kmeans(scale(coords), centers = min(k, nrow(coords)), nstart = 100)
  make_split_list(split(seq_len(nrow(coords)), km$cluster), nrow(coords), paste0("SK", k, "_"))
}

spatial_grid <- function(coords, nx = 3, ny = 3) {
  lat_bin <- cut(coords$Lat, breaks = nx, include.lowest = TRUE, labels = FALSE)
  lon_bin <- cut(coords$Lon, breaks = ny, include.lowest = TRUE, labels = FALSE)
  make_split_list(split(seq_len(nrow(coords)), paste(lat_bin, lon_bin, sep = "_")), nrow(coords), paste0("Grid", nx, "x", ny, "_"))
}

knndm_feature <- function(coords, k = 10) {
  tryCatch({
    coord_id <- paste(round(coords$Lat, 6), round(coords$Lon, 6), sep = ",")
    tpoints_all <- as.data.frame(coords[, c("Lat", "Lon")])
    tpoints <- tpoints_all[!duplicated(coord_id), , drop = FALSE]
    tpoint_ids <- unique(coord_id)
    gx <- seq(min(tpoints$Lat), max(tpoints$Lat), length.out = 45)
    gy <- seq(min(tpoints$Lon), max(tpoints$Lon), length.out = 45)
    predpoints <- expand.grid(Lat = gx, Lon = gy)
    random_assignment <- FALSE
    knn <- withCallingHandlers(
      CAST::knndm(tpoints = tpoints, predpoints = predpoints, space = "feature", k = k, maxp = 0.8),
      message = function(m) {
        if (grepl("random CV assignment", conditionMessage(m), fixed = TRUE)) {
          random_assignment <<- TRUE
          invokeRestart("muffleMessage")
        }
      }
    )
    if (random_assignment) stop("kNNDM returned random assignment")
    list(
      label = "kNNDM_10fold",
      note = paste0("CAST::knndm in Lat/Lon feature space; W=", signif(knn$W, 6)),
      splits = lapply(seq_along(knn$indx_test), function(i) {
        train_ids <- tpoint_ids[sort(unique(knn$indx_train[[i]]))]
        test_ids <- tpoint_ids[sort(unique(knn$indx_test[[i]]))]
        list(
          name = paste0("kNNDM", i),
          train = which(coord_id %in% train_ids),
          test = which(coord_id %in% test_ids)
        )
      })
    )
  }, error = function(e) {
    list(label = "spatial_kmeans_10fold", note = paste0("kNNDM fallback: ", conditionMessage(e)), splits = spatial_kmeans(coords, k))
  })
}

make_mem <- function(coords, k = 8, n_mem = 5) {
  xy <- as.matrix(coords[, c("Lon", "Lat")])
  if (any(duplicated(as.data.frame(xy)))) {
    set.seed(cfg$seed)
    xy <- xy + matrix(rnorm(length(xy), 0, 1e-9), ncol = 2)
  }
  k <- min(k, nrow(xy) - 2)
  knn <- knearneigh(xy, k = k, longlat = FALSE)
  W <- listw2mat(nb2listw(knn2nb(knn), style = "B", zero.policy = TRUE))
  W <- (W + t(W)) / 2
  n <- nrow(W)
  C <- diag(n) - matrix(1 / n, n, n)
  eig <- eigen(C %*% W %*% C, symmetric = TRUE)
  keep <- head(which(eig$values > 1e-8), n_mem)
  mem <- as.data.frame(eig$vectors[, keep, drop = FALSE])
  names(mem) <- paste0("MEM", seq_len(ncol(mem)))
  mem
}

calc_moran <- function(residuals, coords, k = 8) {
  xy <- as.matrix(coords[, c("Lon", "Lat")])
  ok <- is.finite(residuals) & is.finite(xy[, 1]) & is.finite(xy[, 2])
  xy <- xy[ok, , drop = FALSE]
  residuals <- residuals[ok]
  if (nrow(xy) < 5) return(data.frame(moran_observed_I = NA_real_, moran_expected_I = NA_real_, moran_p_value = NA_real_))
  k <- min(k, nrow(xy) - 2)
  if (any(duplicated(as.data.frame(xy)))) {
    set.seed(cfg$seed)
    xy <- xy + matrix(rnorm(length(xy), 0, 1e-9), ncol = 2)
  }
  lw <- nb2listw(knn2nb(knearneigh(xy, k = k, longlat = FALSE)), style = "W", zero.policy = TRUE)
  mt <- moran.test(residuals, lw, randomisation = TRUE, zero.policy = TRUE)
  data.frame(
    moran_observed_I = unname(mt$estimate[["Moran I statistic"]]),
    moran_expected_I = unname(mt$estimate[["Expectation"]]),
    moran_p_value = mt$p.value,
    stringsAsFactors = FALSE
  )
}

median_impute <- function(d) {
  for (nm in names(d)) {
    if (anyNA(d[[nm]])) d[[nm]][is.na(d[[nm]])] <- median(d[[nm]], na.rm = TRUE)
  }
  d
}

aggregate_by_coord <- function(d, target, env_vars, spatial_vars = character(0)) {
  d$coord_id <- paste(round(d$Lat, 6), round(d$Lon, 6), sep = ",")
  vars_mean <- unique(c(target, env_vars, spatial_vars, "Lat", "Lon"))
  rows <- lapply(split(d, d$coord_id), function(z) {
    out <- as.data.frame(lapply(z[, vars_mean, drop = FALSE], function(x) mean(x, na.rm = TRUE)), check.names = FALSE)
    out$biome <- mode_value(z$biome)
    out$n_samples_at_coord <- nrow(z)
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out$biome[is.na(out$biome)] <- "Unknown"
  out
}

keep_individual_samples <- function(d, target, env_vars, spatial_vars = character(0)) {
  keep <- unique(c(target, env_vars, spatial_vars, "Lat", "Lon", "biome", "site"))
  keep <- keep[keep %in% names(d)]
  out <- d[, keep, drop = FALSE]
  out$biome[is.na(out$biome) | out$biome == ""] <- "Unknown"
  out <- out[is.finite(out[[target]]) & is.finite(out$Lat) & is.finite(out$Lon), , drop = FALSE]
  out$n_samples_at_coord <- ave(seq_len(nrow(out)), paste(round(out$Lat, 6), round(out$Lon, 6), sep = ","), FUN = length)
  rownames(out) <- NULL
  out
}

load_td <- function() {
  d <- read.csv(td_file, check.names = FALSE)
  names(d)[names(d) == ""] <- "row_id"
  target <- "Shannon.TD"
  env_vars <- setdiff(names(d)[vapply(d, is.numeric, logical(1))], c("row_id", "X", target, "Lat", "Lon"))
  if (!cfg$allow_future_cropland) env_vars <- env_vars[!grepl("^cropland2100", env_vars)]
  d$biome[is.na(d$biome)] <- "Unknown"
  if (cfg$aggregate_same_coords) {
    agg <- aggregate_by_coord(d, target, env_vars)
  } else {
    agg <- keep_individual_samples(d, target, env_vars)
  }
  list(label = "TD", target = target, data = agg, env_vars = env_vars, spatial_vars = character(0), spatial_source = "MEM")
}

load_gfm <- function() {
  g <- read.csv(gfm_file, check.names = FALSE)
  names(g)[names(g) == ""] <- "row_id"
  td_coords <- read.csv(td_file, check.names = FALSE)[, c("site", "Lat", "Lon"), drop = FALSE]
  td_coords <- td_coords[!duplicated(td_coords$site), , drop = FALSE]
  g <- merge(g, td_coords, by = "site", all.x = FALSE, sort = FALSE)
  target <- "GFM"
  spatial_vars <- grep("^PCNM", names(g), value = TRUE)
  env_vars <- setdiff(names(g)[vapply(g, is.numeric, logical(1))], c("row_id", "X", target, spatial_vars, "Lat", "Lon"))
  if (!cfg$allow_future_cropland) env_vars <- env_vars[!grepl("^cropland2100", env_vars)]
  g$biome[is.na(g$biome)] <- "Unknown"
  if (cfg$aggregate_same_coords) {
    agg <- aggregate_by_coord(g, target, env_vars, spatial_vars)
  } else {
    agg <- keep_individual_samples(g, target, env_vars, spatial_vars)
  }
  list(label = "GFM", target = target, data = agg, env_vars = env_vars, spatial_vars = spatial_vars, spatial_source = "PCNM")
}

prepare_dataset <- function(src) {
  d <- src$data
  x_env <- median_impute(d[, src$env_vars, drop = FALSE])
  nzv <- nearZeroVar(x_env, saveMetrics = TRUE)
  x_env <- x_env[, rownames(nzv)[!nzv$nzv], drop = FALSE]
  if (ncol(x_env) > 1) {
    cors <- cor(x_env, use = "pairwise.complete.obs")
    drop_cor <- findCorrelation(cors, cutoff = cfg$correlation_cutoff, names = TRUE, exact = TRUE)
    x_env <- x_env[, setdiff(names(x_env), drop_cor), drop = FALSE]
  }
  d[, names(x_env)] <- x_env
  list(
    label = src$label,
    target = src$target,
    data = d,
    x_env = x_env,
    env_vars = names(x_env),
    spatial_vars = src$spatial_vars[src$spatial_vars %in% names(d)],
    spatial_source = src$spatial_source,
    y = d[[src$target]],
    coords = d[, c("Lat", "Lon"), drop = FALSE],
    biome = factor(d$biome)
  )
}

rank_vars <- function(ds) {
  fit <- ranger(
    .y ~ .,
    data = data.frame(.y = ds$y, ds$x_env, check.names = FALSE),
    num.trees = cfg$candidate_trees,
    importance = "permutation",
    mtry = max(1, floor(sqrt(ncol(ds$x_env)))),
    seed = cfg$seed
  )
  sort(fit$variable.importance, decreasing = TRUE)
}

build_x <- function(ds, vars, feature_mode) {
  x <- ds$x_env[, vars, drop = FALSE]
  spatial_terms <- character(0)
  if (grepl("biome", feature_mode)) {
    mm <- as.data.frame(model.matrix(~ ds$biome - 1))
    names(mm) <- make.names(gsub("^ds\\$biome", "biome_", names(mm)), unique = TRUE)
    x <- cbind(x, mm)
    spatial_terms <- c(spatial_terms, "biome")
  }
  if (grepl("coord", feature_mode)) {
    coord <- data.frame(
      coord_Lat = ds$coords$Lat,
      coord_Lon = ds$coords$Lon,
      coord_Lat2 = ds$coords$Lat^2,
      coord_Lon2 = ds$coords$Lon^2,
      coord_LatLon = ds$coords$Lat * ds$coords$Lon
    )
    x <- cbind(x, coord)
    spatial_terms <- c(spatial_terms, names(coord))
  }
  if (grepl("spatial", feature_mode)) {
    if (ds$spatial_source == "PCNM" && length(ds$spatial_vars) > 0) {
      sp <- ds$data[, ds$spatial_vars, drop = FALSE]
      x <- cbind(x, sp)
      spatial_terms <- c(spatial_terms, ds$spatial_vars)
    } else {
      mem <- make_mem(ds$coords, cfg$moran_k, cfg$mem_count)
      x <- cbind(x, mem)
      spatial_terms <- c(spatial_terms, names(mem))
    }
  }
  attr(x, "spatial_terms") <- setdiff(spatial_terms, "biome")
  as.data.frame(x, check.names = FALSE)
}

fit_predict <- function(train_x, train_y, test_x, model_kind, trees) {
  if (model_kind == "RF_log1p") {
    y_fit <- log1p(train_y)
  } else {
    y_fit <- train_y
  }
  fit <- ranger(
    .y ~ .,
    data = data.frame(.y = y_fit, train_x, check.names = FALSE),
    num.trees = trees,
    mtry = max(1, floor(sqrt(ncol(train_x)))),
    min.node.size = 5,
    seed = cfg$seed
  )
  pred <- as.numeric(predict(fit, data = test_x)$predictions)
  if (model_kind == "RF_log1p") pred <- pmax(expm1(pred), 0)
  pred
}

full_predict <- function(x, y, model_kind, trees = cfg$final_trees) {
  fit_predict(x, y, x, model_kind, trees)
}

fit_full_model_object <- function(x, y, model_kind, trees = cfg$final_trees, importance = "permutation") {
  if (model_kind == "RF_log1p") {
    y_fit <- log1p(y)
  } else {
    y_fit <- y
  }
  ranger(
    .y ~ .,
    data = data.frame(.y = y_fit, x, check.names = FALSE),
    num.trees = trees,
    mtry = max(1, floor(sqrt(ncol(x)))),
    min.node.size = 5,
    importance = importance,
    seed = cfg$seed
  )
}

predict_full_model_object <- function(model, x, model_kind) {
  pred <- as.numeric(predict(model, data = x)$predictions)
  if (model_kind == "RF_log1p") pred <- pmax(expm1(pred), 0)
  pred
}

cv_predict <- function(x, y, splits, model_kind, trees) {
  rows <- vector("list", length(splits))
  for (i in seq_along(splits)) {
    sp <- splits[[i]]
    pred <- fit_predict(x[sp$train, , drop = FALSE], y[sp$train], x[sp$test, , drop = FALSE], model_kind, trees)
    rows[[i]] <- data.frame(row_id = sp$test, obs = y[sp$test], pred = pred)
  }
  do.call(rbind, rows)
}

score_candidates <- function(d) {
  d$score <- d$mean_R2 + 0.75 * d$best_spatial_R2 - 0.25 * d$mean_RMSE - 0.25 * d$best_spatial_RMSE
  d[order(-d$score, -d$best_spatial_R2, d$best_spatial_RMSE, d$n_env_factors), ]
}

run_one <- function(src) {
  ds <- prepare_dataset(src)
  n_unique_coords <- nrow(unique(ds$coords))
  message("Dataset ", ds$label, ": n=", length(ds$y), ", unique coords=", n_unique_coords, ", env candidates=", length(ds$env_vars))
  fwrite(
    data.frame(
      dataset = ds$label,
      original_target = ds$target,
      n_model_rows = length(ds$y),
      n_unique_coords = n_unique_coords,
      n_env_after_filter = length(ds$env_vars),
      aggregate_same_coords = cfg$aggregate_same_coords,
      allow_future_cropland = cfg$allow_future_cropland
    ),
    file.path(out_dir, paste0(ds$label, "_preprocessing_summary.csv"))
  )
  imp <- rank_vars(ds)
  fwrite(data.frame(variable = names(imp), importance = as.numeric(imp)), file.path(out_dir, paste0(ds$label, "_predictor_importance.csv")))
  fwrite(data.frame(variable = names(imp), importance = as.numeric(imp)), file.path(out_dir, paste0(ds$label, "_screening_importance.csv")))

  ranked <- names(imp)
  sizes <- unique(pmin(length(ranked), c(3, 5, 8, 12, 16, 20, length(ranked))))
  sizes <- sizes[sizes >= 1]
  feature_modes <- c("env_only", "env_plus_biome", "env_plus_coord", "env_plus_biome_coord", "env_plus_spatial", "env_plus_biome_spatial")

  cv_methods <- list(
    random_10fold = random_folds(ds$y, cfg$k_folds, repeats = 1),
    biome_stratified_10fold = stratified_biome_folds(ds$y, ds$biome, cfg$k_folds),
    spatial_kmeans_10fold = spatial_kmeans(ds$coords, cfg$k_folds),
    spatial_grid_3x3 = spatial_grid(ds$coords, 3, 3)
  )
  knn <- knndm_feature(ds$coords, cfg$k_folds)
  cv_methods[[knn$label]] <- knn$splits

  candidate_rows <- list()
  candidate_cv <- list()
  cid <- 1
  for (model_kind in c("RF_raw", "RF_log1p")) {
    for (feature_mode in feature_modes) {
      for (sz in sizes) {
        vars <- ranked[seq_len(sz)]
        x <- build_x(ds, vars, feature_mode)
        spatial_terms <- attr(x, "spatial_terms")
        cv_rows <- list()
        for (method_name in names(cv_methods)) {
          pred <- cv_predict(x, ds$y, cv_methods[[method_name]], model_kind, cfg$candidate_trees)
          met <- metrics_reg(pred$obs, pred$pred)
          cv_rows[[method_name]] <- data.frame(
            dataset = ds$label,
            candidate_key = paste0("cand", cid),
            model_kind = model_kind,
            feature_mode = feature_mode,
            n_env_factors = length(vars),
            environmental_factors = paste(vars, collapse = "; "),
            spatial_terms = paste(spatial_terms, collapse = "; "),
          cv_method = method_name,
          R2 = met$R2,
          RMSE = met$RMSE,
          NRMSE_SD = met$NRMSE_SD,
          NRMSE_range = met$NRMSE_range,
          MAE = met$MAE,
            stringsAsFactors = FALSE
          )
        }
        cv_eval <- do.call(rbind, cv_rows)
        pred_full <- full_predict(x, ds$y, model_kind, cfg$candidate_trees)
        moran <- calc_moran(ds$y - pred_full, ds$coords, cfg$moran_k)
        spatial <- cv_eval[cv_eval$cv_method %in% c("kNNDM_10fold", "spatial_kmeans_10fold", "spatial_grid_3x3"), ]
        best_spatial <- spatial[order(-spatial$R2, spatial$RMSE), ][1, ]
        candidate_rows[[cid]] <- data.frame(
          dataset = ds$label,
          candidate_key = paste0("cand", cid),
          model_kind = model_kind,
          feature_mode = feature_mode,
          n_env_factors = length(vars),
          environmental_factors = paste(vars, collapse = "; "),
          spatial_terms = paste(spatial_terms, collapse = "; "),
          mean_R2 = mean(cv_eval$R2, na.rm = TRUE),
          mean_RMSE = mean(cv_eval$RMSE, na.rm = TRUE),
          best_spatial_method = best_spatial$cv_method,
          best_spatial_R2 = best_spatial$R2,
          best_spatial_RMSE = best_spatial$RMSE,
          moran_observed_I = moran$moran_observed_I,
          moran_expected_I = moran$moran_expected_I,
          moran_p_value = moran$moran_p_value,
          moran_pass = is.na(moran$moran_p_value) || moran$moran_p_value > 0.05,
          stringsAsFactors = FALSE
        )
        candidate_cv[[cid]] <- cv_eval
        cid <- cid + 1
      }
    }
  }
  candidates <- score_candidates(do.call(rbind, candidate_rows))
  fwrite(candidates, file.path(out_dir, paste0(ds$label, "_candidate_summary.csv")))
  fwrite(do.call(rbind, candidate_cv), file.path(out_dir, paste0(ds$label, "_candidate_cv_long.csv")))

  pass <- candidates[candidates$moran_pass & candidates$best_spatial_R2 > 0, ]
  if (nrow(pass) == 0) pass <- candidates[candidates$moran_pass, ]
  if (nrow(pass) == 0) pass <- candidates
  chosen <- pass[order(-pass$score, -pass$best_spatial_R2, pass$best_spatial_RMSE, pass$n_env_factors), ][1, ]

  vars_chosen <- trimws(strsplit(chosen$environmental_factors, ";", fixed = TRUE)[[1]])
  x_chosen <- build_x(ds, vars_chosen, chosen$feature_mode)
  final_fit <- fit_full_model_object(x_chosen, ds$y, chosen$model_kind, cfg$final_trees, importance = "permutation")
  final_imp <- sort(final_fit$variable.importance, decreasing = TRUE)
  final_spatial_terms <- attr(x_chosen, "spatial_terms")
  final_imp_df <- data.frame(
    variable = names(final_imp),
    importance = as.numeric(final_imp),
    role = ifelse(
      names(final_imp) %in% vars_chosen,
      "environmental",
      ifelse(names(final_imp) %in% final_spatial_terms, "spatial", "biome")
    ),
    stringsAsFactors = FALSE
  )
  fwrite(final_imp_df, file.path(out_dir, paste0(ds$label, "_final_model_importance.csv")))

  final_methods <- list(
    LOOCV = loocv_folds(length(ds$y)),
    random_10fold_x5 = random_folds(ds$y, cfg$k_folds, repeats = 5),
    biome_stratified_10fold = stratified_biome_folds(ds$y, ds$biome, cfg$k_folds),
    spatial_kmeans_10fold = spatial_kmeans(ds$coords, cfg$k_folds),
    spatial_grid_3x3 = spatial_grid(ds$coords, 3, 3)
  )
  final_knn <- knndm_feature(ds$coords, cfg$k_folds)
  final_methods[[final_knn$label]] <- final_knn$splits

  final_cv <- lapply(names(final_methods), function(method_name) {
    pred <- cv_predict(x_chosen, ds$y, final_methods[[method_name]], chosen$model_kind, cfg$final_trees)
    met <- metrics_reg(pred$obs, pred$pred)
    data.frame(dataset = ds$label, model_kind = chosen$model_kind, feature_mode = chosen$feature_mode, cv_method = method_name, R2 = met$R2, RMSE = met$RMSE, NRMSE_SD = met$NRMSE_SD, NRMSE_range = met$NRMSE_range, MAE = met$MAE, stringsAsFactors = FALSE)
  })
  final_cv <- do.call(rbind, final_cv)
  fwrite(final_cv, file.path(out_dir, paste0(ds$label, "_final_cv.csv")))

  pred_full <- predict_full_model_object(final_fit, x_chosen, chosen$model_kind)
  full_met <- metrics_reg(ds$y, pred_full)
  moran <- calc_moran(ds$y - pred_full, ds$coords, cfg$moran_k)
  spatial_final <- final_cv[final_cv$cv_method %in% c("kNNDM_10fold", "spatial_kmeans_10fold", "spatial_grid_3x3"), ]
  best_spatial <- spatial_final[order(-spatial_final$R2, spatial_final$RMSE), ][1, ]
  best_overall <- final_cv[order(-final_cv$R2, final_cv$RMSE), ][1, ]
  data.frame(
    dataset = ds$label,
    target = ds$target,
    n = length(ds$y),
    model_type = chosen$model_kind,
    feature_mode = chosen$feature_mode,
    environmental_factors = chosen$environmental_factors,
    spatial_terms = chosen$spatial_terms,
    n_environmental_factors = chosen$n_env_factors,
    full_fit_R2 = full_met$R2,
    full_fit_RMSE = full_met$RMSE,
    full_fit_NRMSE_SD = full_met$NRMSE_SD,
    full_fit_NRMSE_range = full_met$NRMSE_range,
    LOOCV_R2 = final_cv$R2[final_cv$cv_method == "LOOCV"][1],
    LOOCV_RMSE = final_cv$RMSE[final_cv$cv_method == "LOOCV"][1],
    LOOCV_NRMSE_SD = final_cv$NRMSE_SD[final_cv$cv_method == "LOOCV"][1],
    LOOCV_NRMSE_range = final_cv$NRMSE_range[final_cv$cv_method == "LOOCV"][1],
    random_10fold_R2 = final_cv$R2[final_cv$cv_method == "random_10fold_x5"][1],
    random_10fold_RMSE = final_cv$RMSE[final_cv$cv_method == "random_10fold_x5"][1],
    random_10fold_NRMSE_SD = final_cv$NRMSE_SD[final_cv$cv_method == "random_10fold_x5"][1],
    random_10fold_NRMSE_range = final_cv$NRMSE_range[final_cv$cv_method == "random_10fold_x5"][1],
    biome_10fold_R2 = final_cv$R2[final_cv$cv_method == "biome_stratified_10fold"][1],
    biome_10fold_RMSE = final_cv$RMSE[final_cv$cv_method == "biome_stratified_10fold"][1],
    biome_10fold_NRMSE_SD = final_cv$NRMSE_SD[final_cv$cv_method == "biome_stratified_10fold"][1],
    biome_10fold_NRMSE_range = final_cv$NRMSE_range[final_cv$cv_method == "biome_stratified_10fold"][1],
    best_cv_method = best_overall$cv_method,
    best_cv_R2 = best_overall$R2,
    best_cv_RMSE = best_overall$RMSE,
    best_cv_NRMSE_SD = best_overall$NRMSE_SD,
    best_cv_NRMSE_range = best_overall$NRMSE_range,
    best_spatial_method = best_spatial$cv_method,
    best_spatial_R2 = best_spatial$R2,
    best_spatial_RMSE = best_spatial$RMSE,
    best_spatial_NRMSE_SD = best_spatial$NRMSE_SD,
    best_spatial_NRMSE_range = best_spatial$NRMSE_range,
    moran_observed_I = moran$moran_observed_I,
    moran_expected_I = moran$moran_expected_I,
    moran_p_value = moran$moran_p_value,
    moran_pass = is.na(moran$moran_p_value) || moran$moran_p_value > 0.05,
    cv_pass = best_spatial$R2 > 0 && (is.na(moran$moran_p_value) || moran$moran_p_value > 0.05),
    stringsAsFactors = FALSE
  )
}

sources <- list(load_td(), load_gfm())
best <- do.call(rbind, lapply(sources, run_one))
fwrite(best, file.path(out_dir, "TD_GFM_best_model_summary.csv"))
rounded <- best
num <- vapply(rounded, is.numeric, logical(1))
rounded[num] <- lapply(rounded[num], function(x) round(x, 6))
fwrite(rounded, file.path(out_dir, "TD_GFM_best_model_summary_rounded.csv"))
print(rounded, row.names = FALSE)
