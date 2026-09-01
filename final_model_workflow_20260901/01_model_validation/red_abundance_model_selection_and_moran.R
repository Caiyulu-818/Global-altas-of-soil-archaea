#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(caret)
  library(ranger)
  library(spdep)
  library(CAST)
})

set.seed(20260616)

base_dir <- "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F4/map_prediction"
files <- file.path(base_dir, sprintf("env.meta.red%d_scaled.csv", 1:5))
groups <- paste0("red", 1:5)

output_dir <- file.path(getwd(), "model_outputs_red_optimized")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cfg <- list(
  seed = 20260616,
  k_folds = 10,
  screening_trees = 150,
  final_cv_trees = 500,
  final_fit_trees = 2000,
  moran_k = 8,
  mem_count = 5,
  top_refine = 8
)

rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
mae <- function(obs, pred) mean(abs(obs - pred), na.rm = TRUE)
pred_r2 <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  denom <- sum((obs - mean(obs))^2)
  if (denom <= 0) return(NA_real_)
  1 - sum((obs - pred)^2) / denom
}
metric_row <- function(obs, pred) {
  data.frame(
    R2 = pred_r2(obs, pred),
    RMSE = rmse(obs, pred),
    MAE = mae(obs, pred),
    n = sum(is.finite(obs) & is.finite(pred)),
    stringsAsFactors = FALSE
  )
}

safe <- function(x) make.names(x, unique = TRUE)

make_split_list <- function(test_list, n, prefix) {
  lapply(seq_along(test_list), function(i) {
    test <- sort(unique(test_list[[i]]))
    list(name = paste0(prefix, i), train = setdiff(seq_len(n), test), test = test)
  })
}

make_random_folds <- function(y, k = 10, repeats = 1) {
  train_list <- createMultiFolds(y, k = k, times = repeats)
  lapply(names(train_list), function(nm) {
    train <- sort(unique(train_list[[nm]]))
    list(name = nm, train = train, test = setdiff(seq_along(y), train))
  })
}

make_biome_folds <- function(biome, y, k = 10) {
  if (!is.null(biome) && length(unique(biome)) > 1) {
    test_list <- createFolds(as.factor(biome), k = k, list = TRUE, returnTrain = FALSE)
  } else {
    test_list <- createFolds(y, k = k, list = TRUE, returnTrain = FALSE)
  }
  make_split_list(test_list, length(y), "Biome")
}

make_loocv <- function(n) {
  lapply(seq_len(n), function(i) list(name = paste0("LOO", i), train = setdiff(seq_len(n), i), test = i))
}

make_knndm_feature_folds <- function(coords, k = 10) {
  out <- tryCatch({
    tpoints <- as.data.frame(coords)
    names(tpoints) <- c("Lat", "Lon")
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
    if (random_assignment) stop("CAST::knndm returned a random assignment.")
    splits <- lapply(seq_along(knn$indx_test), function(i) {
      list(
        name = paste0("kNNDM", i),
        train = sort(unique(knn$indx_train[[i]])),
        test = sort(unique(knn$indx_test[[i]]))
      )
    })
    list(splits = splits, label = "kNNDM_10fold", note = paste0("CAST::knndm in scaled Lat/Lon feature space; W=", signif(knn$W, 6)))
  }, error = function(e) {
    set.seed(cfg$seed)
    km <- kmeans(scale(coords), centers = k, nstart = 50)
    test_list <- split(seq_len(nrow(coords)), km$cluster)
    list(
      splits = make_split_list(test_list, nrow(coords), "SpatialCluster"),
      label = "spatial_cluster_10fold",
      note = paste0("kNNDM fallback: ", conditionMessage(e))
    )
  })
  out
}

calc_moran <- function(residuals, coords, k = 8) {
  coords <- as.matrix(coords)
  ok <- is.finite(residuals) & is.finite(coords[, 1]) & is.finite(coords[, 2])
  coords <- coords[ok, , drop = FALSE]
  residuals <- residuals[ok]
  if (nrow(coords) < k + 2) k <- max(1, nrow(coords) - 2)
  if (any(duplicated(as.data.frame(coords)))) {
    set.seed(cfg$seed)
    coords <- coords + matrix(rnorm(length(coords), 0, 1e-9), ncol = 2)
  }
  knn <- knearneigh(coords, k = k, longlat = FALSE)
  nb <- knn2nb(knn)
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  mt <- moran.test(residuals, lw, randomisation = TRUE, zero.policy = TRUE)
  data.frame(
    moran_observed_I = unname(mt$estimate[["Moran I statistic"]]),
    moran_expected_I = unname(mt$estimate[["Expectation"]]),
    moran_p_value = mt$p.value,
    moran_weight = paste0(k, "-nearest neighbors, row-standardized"),
    stringsAsFactors = FALSE
  )
}

make_mem <- function(coords, k = 8, n_mem = 5) {
  coords <- as.matrix(coords)
  if (any(duplicated(as.data.frame(coords)))) {
    set.seed(cfg$seed)
    coords <- coords + matrix(rnorm(length(coords), 0, 1e-9), ncol = 2)
  }
  knn <- knearneigh(coords, k = k, longlat = FALSE)
  nb <- knn2nb(knn)
  W <- listw2mat(nb2listw(nb, style = "B", zero.policy = TRUE))
  W <- (W + t(W)) / 2
  n <- nrow(W)
  C <- diag(n) - matrix(1 / n, n, n)
  eig <- eigen(C %*% W %*% C, symmetric = TRUE)
  keep <- head(which(eig$values > 1e-8), n_mem)
  mem <- as.data.frame(eig$vectors[, keep, drop = FALSE])
  names(mem) <- paste0("MEM", seq_len(ncol(mem)))
  mem
}

fit_predict_once <- function(train_x, train_y, test_x, model_kind, mtry, num_trees) {
  mtry <- max(1, min(ncol(train_x), mtry))

  if (model_kind == "RF_raw") {
    train <- data.frame(.y = train_y, train_x, check.names = FALSE)
    fit <- ranger(.y ~ ., data = train, num.trees = num_trees, mtry = mtry, min.node.size = 5, seed = cfg$seed)
    return(as.numeric(predict(fit, data = test_x)$predictions))
  }

  baseline <- min(train_y, na.rm = TRUE)
  shifted <- train_y - baseline

  if (model_kind == "RF_shift_log1p") {
    train <- data.frame(.y = log1p(shifted), train_x, check.names = FALSE)
    fit <- ranger(.y ~ ., data = train, num.trees = num_trees, mtry = mtry, min.node.size = 5, seed = cfg$seed)
    return(as.numeric(expm1(predict(fit, data = test_x)$predictions) + baseline))
  }

  if (model_kind == "RF_hurdle_shift_log1p") {
    tol <- sqrt(.Machine$double.eps)
    present <- factor(ifelse(train_y > baseline + tol, "present", "absent"), levels = c("absent", "present"))
    class_train <- data.frame(.y = present, train_x, check.names = FALSE)
    class_fit <- ranger(.y ~ ., data = class_train, num.trees = num_trees, mtry = mtry, min.node.size = 5, probability = TRUE, seed = cfg$seed)
    prob <- predict(class_fit, data = test_x)$predictions[, "present"]

    pos <- train_y > baseline + tol
    if (sum(pos) < 8) {
      pos_excess <- rep(mean(shifted[pos], na.rm = TRUE), nrow(test_x))
    } else {
      reg_train <- data.frame(.y = log1p(shifted[pos]), train_x[pos, , drop = FALSE], check.names = FALSE)
      reg_fit <- ranger(.y ~ ., data = reg_train, num.trees = num_trees, mtry = mtry, min.node.size = 3, seed = cfg$seed + 17)
      pos_excess <- expm1(predict(reg_fit, data = test_x)$predictions)
    }
    return(as.numeric(baseline + prob * pos_excess))
  }

  stop("Unknown model kind: ", model_kind)
}

cv_predict <- function(x, y, splits, model_kind, mtry, num_trees) {
  rows <- vector("list", length(splits))
  for (i in seq_along(splits)) {
    sp <- splits[[i]]
    set.seed(cfg$seed + i)
    pred <- fit_predict_once(x[sp$train, , drop = FALSE], y[sp$train], x[sp$test, , drop = FALSE], model_kind, mtry, num_trees)
    rows[[i]] <- data.frame(row_id = sp$test, fold = sp$name, obs = y[sp$test], pred = pred)
  }
  do.call(rbind, rows)
}

fit_full_pred <- function(x, y, model_kind, mtry, num_trees) {
  fit_predict_once(x, y, x, model_kind, mtry, num_trees)
}

rank_predictors <- function(x, y) {
  train <- data.frame(.y = y, x, check.names = FALSE)
  fit <- ranger(.y ~ ., data = train, num.trees = cfg$final_fit_trees, mtry = max(1, floor(sqrt(ncol(x)))), importance = "permutation", seed = cfg$seed)
  sort(fit$variable.importance, decreasing = TRUE)
}

score_table <- function(d) {
  key <- paste(d$model_kind, d$feature_mode, d$candidate_id, sep = "__")
  parts <- split(d, key)
  rows <- lapply(parts, function(z) {
    first <- z[1, ]
    dist <- z[z$cv_method %in% c("kNNDM_10fold", "spatial_cluster_10fold"), ][1, ]
    data.frame(
      key = paste(first$model_kind, first$feature_mode, first$candidate_id, sep = "__"),
      model_kind = first$model_kind,
      feature_mode = first$feature_mode,
      candidate_id = first$candidate_id,
      n_environmental_factors = first$n_environmental_factors,
      environmental_factors = first$environmental_factors,
      environmental_factors_safe = first$environmental_factors_safe,
      mean_R2 = mean(z$R2, na.rm = TRUE),
      mean_RMSE = mean(z$RMSE, na.rm = TRUE),
      min_R2 = min(z$R2, na.rm = TRUE),
      distance_R2 = dist$R2,
      distance_RMSE = dist$RMSE,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$rank_mean_R2 <- rank(-out$mean_R2, ties.method = "min")
  out$rank_mean_RMSE <- rank(out$mean_RMSE, ties.method = "min")
  out$rank_distance_R2 <- rank(-out$distance_R2, ties.method = "min")
  out$rank_min_R2 <- rank(-out$min_R2, ties.method = "min")
  out$score <- out$rank_mean_R2 + out$rank_mean_RMSE + 1.5 * out$rank_distance_R2 + 0.5 * out$rank_min_R2
  out[order(out$score, -out$mean_R2, out$mean_RMSE), ]
}

build_x <- function(env_x, biome_x, mem_x, coord_x, vars_safe, feature_mode) {
  out <- env_x[, vars_safe, drop = FALSE]
  if (grepl("biome", feature_mode)) out <- cbind(out, biome_x)
  if (grepl("coord", feature_mode)) out <- cbind(out, coord_x)
  if (grepl("MEM", feature_mode)) out <- cbind(out, mem_x)
  as.data.frame(out, check.names = FALSE)
}

all_best <- list()
all_cv <- list()

for (idx in seq_along(files)) {
  group <- groups[idx]
  file <- files[idx]
  message("Running ", group)
  dat <- read.csv(file, check.names = FALSE)
  target <- grep("^Guild[.]", names(dat), value = TRUE)
  coord_cols <- c("Lat", "Lon")
  env_cols <- setdiff(names(dat)[vapply(dat, is.numeric, logical(1))], c(target, coord_cols))
  env_raw <- dat[, env_cols, drop = FALSE]
  names(env_raw) <- safe(names(env_raw))
  name_map <- setNames(env_cols, names(env_raw))
  y <- dat[[target]]
  coords <- dat[, coord_cols, drop = FALSE]

  biome_dummy <- if ("biome" %in% names(dat)) {
    mm <- model.matrix(~ factor(dat$biome) - 1)
    mm <- as.data.frame(mm)
    names(mm) <- safe(gsub("^factor\\(dat\\$biome\\)", "biome_", names(mm)))
    mm
  } else {
    data.frame()
  }
  mem <- make_mem(coords, k = cfg$moran_k, n_mem = cfg$mem_count)
  coord_features <- data.frame(
    coord_Lat = coords$Lat,
    coord_Lon = coords$Lon,
    coord_Lat2 = coords$Lat^2,
    coord_Lon2 = coords$Lon^2,
    coord_LatLon = coords$Lat * coords$Lon
  )

  imp <- rank_predictors(env_raw, y)
  rank_df <- data.frame(variable_safe = names(imp), variable = name_map[names(imp)], importance = as.numeric(imp), stringsAsFactors = FALSE)
  write.csv(rank_df, file.path(output_dir, paste0(group, "_predictor_importance.csv")), row.names = FALSE)
  ranked <- names(imp)
  sizes <- unique(pmin(length(ranked), c(2, 3, 4, 5, 6, length(ranked))))
  sizes <- sizes[sizes >= 1]

  random_splits <- make_random_folds(y, k = cfg$k_folds, repeats = 1)
  biome_splits <- make_biome_folds(dat$biome, y, k = cfg$k_folds)
  distance_cv <- make_knndm_feature_folds(coords, k = cfg$k_folds)
  distance_splits <- distance_cv$splits
  screening_methods <- list(random_10fold = random_splits, biome_stratified_10fold = biome_splits)
  screening_methods[[distance_cv$label]] <- distance_splits

  screening_rows <- list()
  row_id <- 1
  feature_modes <- c(
    "env_only", "env_plus_biome", "env_plus_MEM", "env_plus_biome_MEM",
    "env_plus_coord", "env_plus_biome_coord", "env_plus_coord_MEM", "env_plus_biome_coord_MEM"
  )
  for (model_kind in c("RF_raw", "RF_shift_log1p", "RF_hurdle_shift_log1p")) {
    for (feature_mode in feature_modes) {
      if (ncol(biome_dummy) == 0 && grepl("biome", feature_mode)) next
      for (sz in sizes) {
        vars_safe <- ranked[seq_len(sz)]
        x <- build_x(env_raw, biome_dummy, mem, coord_features, vars_safe, feature_mode)
        mtry <- max(1, floor(sqrt(ncol(x))))
        for (method_name in names(screening_methods)) {
          pred <- cv_predict(x, y, screening_methods[[method_name]], model_kind, mtry, cfg$screening_trees)
          met <- metric_row(pred$obs, pred$pred)
          screening_rows[[row_id]] <- data.frame(
            group = group,
            target = target,
            model_kind = model_kind,
            feature_mode = feature_mode,
            candidate_id = paste0("top", sz),
            n_environmental_factors = sz,
            environmental_factors_safe = paste(vars_safe, collapse = "; "),
            environmental_factors = paste(name_map[vars_safe], collapse = "; "),
            cv_method = method_name,
            mtry = mtry,
            R2 = met$R2,
            RMSE = met$RMSE,
            MAE = met$MAE,
            stringsAsFactors = FALSE
          )
          row_id <- row_id + 1
        }
      }
    }
  }
  screening <- do.call(rbind, screening_rows)
  write.csv(screening, file.path(output_dir, paste0(group, "_screening_cv_results.csv")), row.names = FALSE)
  screening_score <- score_table(screening)
  write.csv(screening_score, file.path(output_dir, paste0(group, "_screening_scores.csv")), row.names = FALSE)

  refine <- head(screening_score, cfg$top_refine)
  candidate_rows <- list()
  candidate_cv <- list()
  cand_id <- 1
  for (i in seq_len(nrow(refine))) {
    row <- refine[i, ]
    vars_safe <- strsplit(row$environmental_factors_safe, "; ", fixed = TRUE)[[1]]
    x <- build_x(env_raw, biome_dummy, mem, coord_features, vars_safe, row$feature_mode)
    mtry_grid <- unique(pmax(1, pmin(ncol(x), c(1, floor(sqrt(ncol(x))), floor(ncol(x) / 3), ncol(x)))))
    for (mtry in mtry_grid) {
      full_pred <- fit_full_pred(x, y, row$model_kind, mtry, cfg$final_fit_trees)
      full_met <- metric_row(y, full_pred)
      moran <- calc_moran(y - full_pred, coords, k = cfg$moran_k)

      quick_methods <- list(random_10fold = random_splits, biome_stratified_10fold = biome_splits)
      quick_methods[[distance_cv$label]] <- distance_splits
      cv_rows <- list()
      for (method_name in names(quick_methods)) {
        pred <- cv_predict(x, y, quick_methods[[method_name]], row$model_kind, mtry, cfg$final_cv_trees)
        met <- metric_row(pred$obs, pred$pred)
        cv_rows[[method_name]] <- data.frame(
          candidate_key = paste0("cand", cand_id),
          cv_method = method_name,
          R2 = met$R2,
          RMSE = met$RMSE,
          MAE = met$MAE,
          stringsAsFactors = FALSE
        )
      }
      cv_eval <- do.call(rbind, cv_rows)
      dist_r2 <- cv_eval$R2[cv_eval$cv_method == distance_cv$label][1]
      dist_rmse <- cv_eval$RMSE[cv_eval$cv_method == distance_cv$label][1]
      candidate_rows[[cand_id]] <- data.frame(
        group = group,
        target = target,
        candidate_key = paste0("cand", cand_id),
        model_kind = row$model_kind,
        feature_mode = row$feature_mode,
        n_environmental_factors = row$n_environmental_factors,
        environmental_factors = row$environmental_factors,
        environmental_factors_safe = row$environmental_factors_safe,
        spatial_terms = paste(c(
          if (grepl("coord", row$feature_mode)) names(coord_features) else character(0),
          if (grepl("MEM", row$feature_mode)) names(mem) else character(0)
        ), collapse = "; "),
        mtry = mtry,
        full_fit_R2 = full_met$R2,
        full_fit_RMSE = full_met$RMSE,
        mean_cv_R2 = mean(cv_eval$R2, na.rm = TRUE),
        mean_cv_RMSE = mean(cv_eval$RMSE, na.rm = TRUE),
        min_cv_R2 = min(cv_eval$R2, na.rm = TRUE),
        distance_cv_method = distance_cv$label,
        distance_cv_note = distance_cv$note,
        distance_cv_R2 = dist_r2,
        distance_cv_RMSE = dist_rmse,
        moran_observed_I = moran$moran_observed_I,
        moran_expected_I = moran$moran_expected_I,
        moran_p_value = moran$moran_p_value,
        moran_pass = moran$moran_p_value > 0.05,
        stringsAsFactors = FALSE
      )
      candidate_cv[[cand_id]] <- cv_eval
      cand_id <- cand_id + 1
    }
  }

  candidates <- do.call(rbind, candidate_rows)
  candidates$cv_pass <- candidates$mean_cv_R2 > 0 & candidates$distance_cv_R2 > 0
  candidates$rank_mean_R2 <- rank(-candidates$mean_cv_R2, ties.method = "min")
  candidates$rank_distance_R2 <- rank(-candidates$distance_cv_R2, ties.method = "min")
  candidates$rank_mean_RMSE <- rank(candidates$mean_cv_RMSE, ties.method = "min")
  candidates$rank_abs_I <- rank(abs(candidates$moran_observed_I), ties.method = "min")
  candidates$final_score <- candidates$rank_mean_R2 + 1.5 * candidates$rank_distance_R2 + candidates$rank_mean_RMSE + candidates$rank_abs_I
  write.csv(candidates, file.path(output_dir, paste0(group, "_candidate_comparison.csv")), row.names = FALSE)

  pass <- candidates[candidates$moran_pass & candidates$cv_pass, , drop = FALSE]
  if (nrow(pass) == 0) {
    pass <- candidates[candidates$moran_pass, , drop = FALSE]
  }
  if (nrow(pass) == 0) {
    pass <- candidates
  }
  chosen <- pass[order(pass$final_score, -pass$mean_cv_R2, pass$mean_cv_RMSE), ][1, ]
  vars_safe <- strsplit(chosen$environmental_factors_safe, "; ", fixed = TRUE)[[1]]
  x_chosen <- build_x(env_raw, biome_dummy, mem, coord_features, vars_safe, chosen$feature_mode)

  final_methods <- list(
    LOOCV = make_loocv(length(y)),
    random_10fold_x5 = make_random_folds(y, k = cfg$k_folds, repeats = 5),
    biome_stratified_10fold = biome_splits
  )
  final_methods[[distance_cv$label]] <- distance_splits

  final_cv_rows <- list()
  final_pred_rows <- list()
  for (method_name in names(final_methods)) {
    message("  ", group, " final CV: ", method_name)
    pred <- cv_predict(x_chosen, y, final_methods[[method_name]], chosen$model_kind, chosen$mtry, cfg$final_cv_trees)
    met <- metric_row(pred$obs, pred$pred)
    final_cv_rows[[method_name]] <- data.frame(
      group = group,
      target = target,
      model_kind = chosen$model_kind,
      feature_mode = chosen$feature_mode,
      cv_method = method_name,
      R2 = met$R2,
      RMSE = met$RMSE,
      MAE = met$MAE,
      mtry = chosen$mtry,
      stringsAsFactors = FALSE
    )
    pred$group <- group
    pred$target <- target
    pred$cv_method <- method_name
    final_pred_rows[[method_name]] <- pred
  }
  final_cv <- do.call(rbind, final_cv_rows)
  final_pred <- fit_full_pred(x_chosen, y, chosen$model_kind, chosen$mtry, cfg$final_fit_trees)
  full_met <- metric_row(y, final_pred)
  final_moran <- calc_moran(y - final_pred, coords, k = cfg$moran_k)
  best_cv <- final_cv[order(-final_cv$R2, final_cv$RMSE), ][1, ]

  best <- data.frame(
    group = group,
    target = target,
    model_type = chosen$model_kind,
    feature_mode = chosen$feature_mode,
    environmental_factors = chosen$environmental_factors,
    spatial_terms = chosen$spatial_terms,
    n_environmental_factors = chosen$n_environmental_factors,
    mtry = chosen$mtry,
    full_fit_R2 = full_met$R2,
    full_fit_RMSE = full_met$RMSE,
    best_cv_method = best_cv$cv_method,
    best_cv_R2 = best_cv$R2,
    best_cv_RMSE = best_cv$RMSE,
    random_10fold_R2 = final_cv$R2[final_cv$cv_method == "random_10fold_x5"][1],
    random_10fold_RMSE = final_cv$RMSE[final_cv$cv_method == "random_10fold_x5"][1],
    biome_10fold_R2 = final_cv$R2[final_cv$cv_method == "biome_stratified_10fold"][1],
    biome_10fold_RMSE = final_cv$RMSE[final_cv$cv_method == "biome_stratified_10fold"][1],
    LOOCV_R2 = final_cv$R2[final_cv$cv_method == "LOOCV"][1],
    LOOCV_RMSE = final_cv$RMSE[final_cv$cv_method == "LOOCV"][1],
    distance_cv_method = distance_cv$label,
    distance_cv_R2 = final_cv$R2[final_cv$cv_method == distance_cv$label][1],
    distance_cv_RMSE = final_cv$RMSE[final_cv$cv_method == distance_cv$label][1],
    distance_cv_note = distance_cv$note,
    moran_observed_I = final_moran$moran_observed_I,
    moran_expected_I = final_moran$moran_expected_I,
    moran_p_value = final_moran$moran_p_value,
    moran_pass = final_moran$moran_p_value > 0.05,
    cv_pass = mean(final_cv$R2, na.rm = TRUE) > 0 & final_cv$R2[final_cv$cv_method == distance_cv$label][1] > 0,
    stringsAsFactors = FALSE
  )

  write.csv(final_cv, file.path(output_dir, paste0(group, "_final_cv_results.csv")), row.names = FALSE)
  write.csv(do.call(rbind, final_pred_rows), file.path(output_dir, paste0(group, "_final_cv_predictions.csv")), row.names = FALSE)
  all_best[[group]] <- best
  all_cv[[group]] <- final_cv
}

best_summary <- do.call(rbind, all_best)
final_cv_long <- do.call(rbind, all_cv)
write.csv(best_summary, file.path(output_dir, "red_optimized_best_model_summary.csv"), row.names = FALSE)
write.csv(final_cv_long, file.path(output_dir, "red_optimized_final_cv_long.csv"), row.names = FALSE)

rounded <- best_summary
num <- vapply(rounded, is.numeric, logical(1))
rounded[num] <- lapply(rounded[num], function(x) round(x, 6))
write.csv(rounded, file.path(output_dir, "red_optimized_best_model_summary_rounded.csv"), row.names = FALSE)

message("Done: ", file.path(output_dir, "red_optimized_best_model_summary.csv"))
