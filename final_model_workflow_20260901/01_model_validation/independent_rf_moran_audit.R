#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(caret)
  library(randomForest)
  library(spdep)
  library(sf)
  library(CAST)
})

set.seed(20260615)

base_dir <- "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F4/map_prediction"
files <- file.path(base_dir, sprintf("env.meta.red%d_scaled.csv", 1:5))
groups <- paste0("red", 1:5)

output_dir <- file.path(getwd(), "model_outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
old_outputs <- list.files(
  output_dir,
  pattern = "^(red[1-5]_.*|rf_group_model_.*\\.(csv|md))$",
  full.names = TRUE
)
if (length(old_outputs) > 0) {
  unlink(old_outputs)
}

cfg <- list(
  seed = 20260615,
  screening_ntree = 200,
  tuning_ntree = 400,
  final_cv_ntree = 500,
  final_fit_ntree = 2000,
  screening_repeats = 2,
  tuning_repeats = 5,
  k_folds = 10,
  n_loocv_refine = 3,
  moran_k = 8
)

rmse <- function(obs, pred) {
  sqrt(mean((obs - pred)^2, na.rm = TRUE))
}

mae <- function(obs, pred) {
  mean(abs(obs - pred), na.rm = TRUE)
}

pred_r2 <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  denom <- sum((obs - mean(obs))^2)
  if (denom == 0) {
    return(NA_real_)
  }
  1 - sum((obs - pred)^2) / denom
}

metric_row <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  data.frame(
    RMSE = rmse(obs, pred),
    R2 = pred_r2(obs, pred),
    MAE = mae(obs, pred),
    n_predictions = length(obs),
    stringsAsFactors = FALSE
  )
}

make_split_list <- function(test_list, n) {
  lapply(seq_along(test_list), function(i) {
    test_idx <- sort(unique(test_list[[i]]))
    list(
      name = paste0("Fold", i),
      train = setdiff(seq_len(n), test_idx),
      test = test_idx
    )
  })
}

make_repeated_folds <- function(y, k = 10, repeats = 2) {
  train_list <- createMultiFolds(y, k = k, times = repeats)
  lapply(names(train_list), function(nm) {
    train_idx <- sort(unique(train_list[[nm]]))
    list(
      name = nm,
      train = train_idx,
      test = setdiff(seq_along(y), train_idx)
    )
  })
}

make_stratified_biome_folds <- function(df, y, k = 10) {
  if ("biome" %in% names(df) && length(unique(df$biome)) > 1) {
    test_list <- createFolds(as.factor(df$biome), k = k, list = TRUE, returnTrain = FALSE)
  } else {
    test_list <- createFolds(y, k = k, list = TRUE, returnTrain = FALSE)
  }
  make_split_list(test_list, nrow(df))
}

make_spatial_cluster_folds <- function(df, coord_cols, k = 10) {
  coords <- as.matrix(df[, coord_cols, drop = FALSE])
  storage.mode(coords) <- "double"
  set.seed(cfg$seed)
  km <- kmeans(coords, centers = k, nstart = 50)
  test_list <- split(seq_len(nrow(df)), km$cluster)
  make_split_list(test_list, nrow(df))
}

make_knndm_folds <- function(df, coord_cols, k = 10) {
  tryCatch({
    random_assignment_warning <- FALSE
    knn <- withCallingHandlers({
      df_sf <- st_as_sf(df, coords = coord_cols, crs = 4326, remove = FALSE)
      CAST::knndm(tpoints = df_sf, predpoints = df_sf, k = k)
    }, message = function(m) {
      if (grepl("random CV assignment", conditionMessage(m), fixed = TRUE)) {
        random_assignment_warning <<- TRUE
        invokeRestart("muffleMessage")
      }
    }, warning = function(w) {
      if (grepl("random CV assignment", conditionMessage(w), fixed = TRUE)) {
        random_assignment_warning <<- TRUE
        invokeRestart("muffleWarning")
      }
    })
    if (random_assignment_warning) {
      stop("CAST::knndm returned a random CV assignment instead of spatially separated folds.")
    }
    if (is.null(knn$indx_train) || is.null(knn$indx_test)) {
      stop("CAST::knndm did not return indx_train/indx_test.")
    }
    splits <- lapply(seq_along(knn$indx_test), function(i) {
      list(
        name = paste0("kNNDM", i),
        train = sort(unique(knn$indx_train[[i]])),
        test = sort(unique(knn$indx_test[[i]]))
      )
    })
    list(
      splits = splits,
      label = "spatial_kNNDM_10fold",
      note = "CAST::knndm spatial folds"
    )
  }, error = function(e) {
    message("  CAST::knndm failed; using coordinate-cluster spatial folds instead: ", conditionMessage(e))
    list(
      splits = make_spatial_cluster_folds(df, coord_cols, k = k),
      label = "spatial_cluster_10fold",
      note = paste0("CAST::knndm fallback: ", conditionMessage(e))
    )
  })
}

make_loocv_folds <- function(n) {
  lapply(seq_len(n), function(i) {
    list(
      name = paste0("LOO", i),
      train = setdiff(seq_len(n), i),
      test = i
    )
  })
}

rf_cv_predict <- function(df, target_col, vars, splits, mtry, ntree, seed_offset = 0) {
  y <- df[[target_col]]
  out <- vector("list", length(splits))
  mtry <- max(1, min(length(vars), mtry))

  for (i in seq_along(splits)) {
    sp <- splits[[i]]
    set.seed(cfg$seed + seed_offset + i)
    fit <- randomForest(
      x = df[sp$train, vars, drop = FALSE],
      y = y[sp$train],
      ntree = ntree,
      mtry = mtry,
      importance = FALSE
    )
    pred <- predict(fit, newdata = df[sp$test, vars, drop = FALSE])
    out[[i]] <- data.frame(
      row_id = sp$test,
      fold = sp$name,
      obs = y[sp$test],
      pred = as.numeric(pred),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, out)
}

score_candidate_table <- function(cv_results) {
  parts <- split(cv_results, cv_results$candidate_id)
  rows <- lapply(parts, function(d) {
    first <- d[1, , drop = FALSE]
    data.frame(
      candidate_id = first$candidate_id,
      n_predictors = first$n_predictors,
      predictors = first$predictors,
      mean_RMSE = mean(d$RMSE, na.rm = TRUE),
      mean_R2 = mean(d$R2, na.rm = TRUE),
      min_R2 = min(d$R2, na.rm = TRUE),
      max_RMSE = max(d$RMSE, na.rm = TRUE),
      n_cv_methods = length(unique(d$cv_method)),
      stringsAsFactors = FALSE
    )
  })
  score <- do.call(rbind, rows)
  score$rank_RMSE <- rank(score$mean_RMSE, ties.method = "min")
  score$rank_R2 <- rank(-score$mean_R2, ties.method = "min")
  score$rank_min_R2 <- rank(-score$min_R2, ties.method = "min")
  score$selection_score <- score$rank_RMSE + score$rank_R2 + 0.5 * score$rank_min_R2
  score[order(score$selection_score, -score$mean_R2, score$mean_RMSE, score$n_predictors), ]
}

calc_moran <- function(residuals, coords, k = 8) {
  coords <- as.matrix(coords)
  storage.mode(coords) <- "double"
  complete <- is.finite(residuals) & is.finite(coords[, 1]) & is.finite(coords[, 2])
  residuals <- residuals[complete]
  coords <- coords[complete, , drop = FALSE]

  if (nrow(coords) < (k + 2)) {
    k <- max(1, nrow(coords) - 2)
  }

  if (any(duplicated(as.data.frame(coords)))) {
    set.seed(cfg$seed)
    coords <- coords + matrix(rnorm(length(coords), mean = 0, sd = 1e-9), ncol = 2)
  }

  knn <- knearneigh(coords, k = k, longlat = FALSE)
  nb <- knn2nb(knn)
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  test <- moran.test(residuals, lw, randomisation = TRUE, zero.policy = TRUE)

  data.frame(
    moran_weight = paste0(k, "-nearest neighbors, row-standardized"),
    moran_observed_I = unname(test$estimate[["Moran I statistic"]]),
    moran_expected_I = unname(test$estimate[["Expectation"]]),
    moran_p_value = test$p.value,
    moran_alternative = test$alternative,
    stringsAsFactors = FALSE
  )
}

run_group <- function(file, group_name) {
  message("Running ", group_name, " from ", basename(file))
  dat <- read.csv(file, check.names = FALSE)

  target_col <- grep("^Guild\\.", names(dat), value = TRUE)
  if (length(target_col) != 1) {
    stop("Expected exactly one Guild.* column in ", file)
  }

  lat_col <- grep("^lat$", names(dat), ignore.case = TRUE, value = TRUE)[1]
  lon_col <- grep("^lon$", names(dat), ignore.case = TRUE, value = TRUE)[1]
  coord_cols <- c(lon_col, lat_col)
  if (any(is.na(coord_cols))) {
    stop("Could not find Lat/Lon columns in ", file)
  }

  excluded_always <- c(target_col, lat_col, lon_col)
  candidate_cols <- setdiff(names(dat), excluded_always)
  numeric_predictors <- candidate_cols[vapply(dat[candidate_cols], is.numeric, logical(1))]
  excluded_non_numeric <- setdiff(candidate_cols, numeric_predictors)

  keep_cols <- unique(c(numeric_predictors, target_col, lat_col, lon_col, "biome"))
  keep_cols <- keep_cols[keep_cols %in% names(dat)]
  df <- dat[, keep_cols, drop = FALSE]
  df <- df[complete.cases(df), , drop = FALSE]

  nzv <- nearZeroVar(df[, numeric_predictors, drop = FALSE], saveMetrics = TRUE)
  predictor_pool <- rownames(nzv)[!nzv$nzv]
  if (length(predictor_pool) < 1) {
    stop("No usable numeric environmental predictors after near-zero variance filtering for ", group_name)
  }

  y <- df[[target_col]]
  x <- df[, predictor_pool, drop = FALSE]

  message("  Ranking predictors by full random forest permutation importance")
  set.seed(cfg$seed)
  rank_rf <- randomForest(
    x = x,
    y = y,
    ntree = cfg$final_fit_ntree,
    mtry = max(1, floor(ncol(x) / 3)),
    importance = TRUE
  )
  rank_importance <- importance(rank_rf, type = 1)
  rank_importance <- data.frame(
    variable = rownames(rank_importance),
    permutation_importance = as.numeric(rank_importance[, 1]),
    group = group_name,
    target = target_col,
    stringsAsFactors = FALSE
  )
  rank_importance <- rank_importance[order(rank_importance$permutation_importance, decreasing = TRUE), ]
  ranked_vars <- rank_importance$variable

  repeated_splits <- make_repeated_folds(y, k = cfg$k_folds, repeats = cfg$screening_repeats)
  biome_splits <- make_stratified_biome_folds(df, y, k = cfg$k_folds)
  spatial_cv <- make_knndm_folds(df, coord_cols, k = cfg$k_folds)
  spatial_splits <- spatial_cv$splits
  spatial_method_name <- spatial_cv$label
  spatial_cv_note <- spatial_cv$note
  screening_methods <- list(
    repeated_10fold = repeated_splits,
    biome_stratified_10fold = biome_splits
  )
  screening_methods[[spatial_method_name]] <- spatial_splits

  message("  Screening top-k predictor combinations with repeated, biome-stratified, and spatial CV")
  screening_rows <- list()
  row_id <- 1
  for (k_vars in seq_along(ranked_vars)) {
    vars <- ranked_vars[seq_len(k_vars)]
    screen_mtry <- max(1, floor(length(vars) / 3))
    candidate_id <- paste0("top", k_vars)
    for (method_name in names(screening_methods)) {
      pred_df <- rf_cv_predict(
        df = df,
        target_col = target_col,
        vars = vars,
        splits = screening_methods[[method_name]],
        mtry = screen_mtry,
        ntree = cfg$screening_ntree,
        seed_offset = k_vars * 1000
      )
      met <- metric_row(pred_df$obs, pred_df$pred)
      screening_rows[[row_id]] <- data.frame(
        group = group_name,
        target = target_col,
        stage = "screening",
        candidate_id = candidate_id,
        n_predictors = length(vars),
        predictors = paste(vars, collapse = "; "),
        cv_method = method_name,
        mtry = screen_mtry,
        ntree = cfg$screening_ntree,
        RMSE = met$RMSE,
        R2 = met$R2,
        MAE = met$MAE,
        n_predictions = met$n_predictions,
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1
    }
  }

  screening_results <- do.call(rbind, screening_rows)
  coarse_score <- score_candidate_table(screening_results)
  refine_ids <- head(coarse_score$candidate_id, cfg$n_loocv_refine)

  message("  Refining top ", length(refine_ids), " combinations with LOOCV")
  loocv_splits <- make_loocv_folds(nrow(df))
  loocv_rows <- list()
  for (i in seq_along(refine_ids)) {
    cid <- refine_ids[i]
    vars <- strsplit(coarse_score$predictors[coarse_score$candidate_id == cid][1], "; ", fixed = TRUE)[[1]]
    loocv_mtry <- max(1, floor(length(vars) / 3))
    pred_df <- rf_cv_predict(
      df = df,
      target_col = target_col,
      vars = vars,
      splits = loocv_splits,
      mtry = loocv_mtry,
      ntree = cfg$final_cv_ntree,
      seed_offset = 9000 + i * 1000
    )
    met <- metric_row(pred_df$obs, pred_df$pred)
    loocv_rows[[i]] <- data.frame(
      group = group_name,
      target = target_col,
      stage = "loocv_refinement",
      candidate_id = cid,
      n_predictors = length(vars),
      predictors = paste(vars, collapse = "; "),
      cv_method = "LOOCV",
      mtry = loocv_mtry,
      ntree = cfg$final_cv_ntree,
      RMSE = met$RMSE,
      R2 = met$R2,
      MAE = met$MAE,
      n_predictions = met$n_predictions,
      stringsAsFactors = FALSE
    )
  }
  loocv_results <- do.call(rbind, loocv_rows)

  all_selection_results <- rbind(screening_results, loocv_results)
  refined_results <- all_selection_results[all_selection_results$candidate_id %in% refine_ids, , drop = FALSE]
  refined_score <- score_candidate_table(refined_results)
  selected_vars <- strsplit(refined_score$predictors[1], "; ", fixed = TRUE)[[1]]

  message("  Selected predictors: ", paste(selected_vars, collapse = ", "))
  message("  Tuning final RF mtry by repeated 10-fold CV")
  tuning_splits <- make_repeated_folds(y, k = cfg$k_folds, repeats = cfg$tuning_repeats)
  tuning_rows <- list()
  for (m in seq_len(length(selected_vars))) {
    pred_df <- rf_cv_predict(
      df = df,
      target_col = target_col,
      vars = selected_vars,
      splits = tuning_splits,
      mtry = m,
      ntree = cfg$tuning_ntree,
      seed_offset = 20000 + m * 1000
    )
    met <- metric_row(pred_df$obs, pred_df$pred)
    tuning_rows[[m]] <- data.frame(
      group = group_name,
      target = target_col,
      mtry = m,
      ntree = cfg$tuning_ntree,
      RMSE = met$RMSE,
      R2 = met$R2,
      MAE = met$MAE,
      n_predictions = met$n_predictions,
      stringsAsFactors = FALSE
    )
  }
  tuning_results <- do.call(rbind, tuning_rows)
  tuning_results <- tuning_results[order(-tuning_results$R2, tuning_results$RMSE, tuning_results$mtry), ]
  best_mtry <- tuning_results$mtry[1]

  final_methods <- list(
    LOOCV = loocv_splits,
    repeated_10fold_x5 = tuning_splits,
    biome_stratified_10fold = biome_splits
  )
  final_methods[[spatial_method_name]] <- spatial_splits

  message("  Evaluating final selected model with LOOCV, repeated 10-fold, biome-stratified 10-fold, and spatial CV")
  final_cv_rows <- list()
  final_pred_rows <- list()
  for (method_name in names(final_methods)) {
    pred_df <- rf_cv_predict(
      df = df,
      target_col = target_col,
      vars = selected_vars,
      splits = final_methods[[method_name]],
      mtry = best_mtry,
      ntree = cfg$final_cv_ntree,
      seed_offset = 30000 + length(final_cv_rows) * 1000
    )
    pred_df$group <- group_name
    pred_df$target <- target_col
    pred_df$cv_method <- method_name
    pred_df$mtry <- best_mtry
    pred_df$ntree <- cfg$final_cv_ntree
    final_pred_rows[[method_name]] <- pred_df

    met <- metric_row(pred_df$obs, pred_df$pred)
    final_cv_rows[[method_name]] <- data.frame(
      group = group_name,
      target = target_col,
      model_type = "Random Forest regression",
      cv_method = method_name,
      mtry = best_mtry,
      ntree = cfg$final_cv_ntree,
      RMSE = met$RMSE,
      R2 = met$R2,
      MAE = met$MAE,
      n_predictions = met$n_predictions,
      predictors = paste(selected_vars, collapse = "; "),
      stringsAsFactors = FALSE
    )
  }
  final_cv_results <- do.call(rbind, final_cv_rows)
  final_cv_predictions <- do.call(rbind, final_pred_rows)

  message("  Fitting final full-data RF and calculating Moran's I")
  set.seed(cfg$seed)
  final_fit <- randomForest(
    x = df[, selected_vars, drop = FALSE],
    y = y,
    ntree = cfg$final_fit_ntree,
    mtry = best_mtry,
    importance = TRUE
  )
  final_pred <- predict(final_fit, newdata = df[, selected_vars, drop = FALSE])
  final_resid <- y - final_pred
  final_fit_metrics <- metric_row(y, final_pred)
  moran <- calc_moran(final_resid, df[, coord_cols, drop = FALSE], k = cfg$moran_k)

  final_importance <- importance(final_fit, type = 1)
  final_importance <- data.frame(
    group = group_name,
    target = target_col,
    variable = rownames(final_importance),
    permutation_importance = as.numeric(final_importance[, 1]),
    stringsAsFactors = FALSE
  )
  final_importance <- final_importance[order(final_importance$permutation_importance, decreasing = TRUE), ]

  write.csv(rank_importance, file.path(output_dir, paste0(group_name, "_initial_rf_importance_ranking.csv")), row.names = FALSE)
  write.csv(all_selection_results, file.path(output_dir, paste0(group_name, "_candidate_combination_cv_results.csv")), row.names = FALSE)
  write.csv(coarse_score, file.path(output_dir, paste0(group_name, "_coarse_candidate_scores.csv")), row.names = FALSE)
  write.csv(refined_score, file.path(output_dir, paste0(group_name, "_refined_candidate_scores.csv")), row.names = FALSE)
  write.csv(tuning_results, file.path(output_dir, paste0(group_name, "_final_mtry_tuning_results.csv")), row.names = FALSE)
  write.csv(final_cv_results, file.path(output_dir, paste0(group_name, "_final_cv_results.csv")), row.names = FALSE)
  write.csv(final_cv_predictions, file.path(output_dir, paste0(group_name, "_final_cv_predictions.csv")), row.names = FALSE)
  write.csv(final_importance, file.path(output_dir, paste0(group_name, "_final_variable_importance.csv")), row.names = FALSE)

  cv_wide <- reshape(
    final_cv_results[, c("group", "cv_method", "RMSE", "R2", "MAE")],
    idvar = "group",
    timevar = "cv_method",
    direction = "wide"
  )

  summary <- data.frame(
    group = group_name,
    source_file = basename(file),
    target = target_col,
    model_type = "Random Forest regression",
    final_fit_trees = cfg$final_fit_ntree,
    cv_trees = cfg$final_cv_ntree,
    best_mtry = best_mtry,
    candidate_numeric_predictors = paste(predictor_pool, collapse = "; "),
    excluded_predictor_columns = paste(c(lat_col, lon_col, excluded_non_numeric), collapse = "; "),
    selected_predictors = paste(selected_vars, collapse = "; "),
    n_selected_predictors = length(selected_vars),
    n_rows_used = nrow(df),
    selection_rule = "Choose candidate with best combined rank: high mean CV R2, low mean CV RMSE, and high worst-case CV R2 across repeated 10-fold, biome-stratified 10-fold, spatial kNNDM, and LOOCV refinement.",
    selection_cv_methods = paste0("repeated_10fold; biome_stratified_10fold; ", spatial_method_name, "; LOOCV refinement for top candidates"),
    spatial_cv_method = spatial_method_name,
    spatial_cv_note = spatial_cv_note,
    final_tuning_method = "mtry selected by repeated 10-fold CV x 5 repeats",
    selection_mean_RMSE = refined_score$mean_RMSE[1],
    selection_mean_R2 = refined_score$mean_R2[1],
    selection_min_R2 = refined_score$min_R2[1],
    final_fit_RMSE = final_fit_metrics$RMSE,
    final_fit_R2 = final_fit_metrics$R2,
    moran_residual_source = "full-data final random forest residuals",
    moran_weight = moran$moran_weight,
    moran_observed_I = moran$moran_observed_I,
    moran_expected_I = moran$moran_expected_I,
    moran_p_value = moran$moran_p_value,
    stringsAsFactors = FALSE
  )

  merge(summary, cv_wide, by = "group", all.x = TRUE)
}

all_results <- do.call(rbind, Map(run_group, files, groups))
write.csv(all_results, file.path(output_dir, "rf_group_model_summary_wide.csv"), row.names = FALSE)

all_final_cv <- do.call(rbind, lapply(groups, function(g) {
  read.csv(file.path(output_dir, paste0(g, "_final_cv_results.csv")), check.names = FALSE)
}))
write.csv(all_final_cv, file.path(output_dir, "rf_group_model_final_cv_long.csv"), row.names = FALSE)

all_selection <- do.call(rbind, lapply(groups, function(g) {
  read.csv(file.path(output_dir, paste0(g, "_candidate_combination_cv_results.csv")), check.names = FALSE)
}))
write.csv(all_selection, file.path(output_dir, "rf_group_candidate_combination_cv_long.csv"), row.names = FALSE)

compact_cols <- c(
  "group", "target", "model_type", "selected_predictors", "best_mtry",
  "selection_mean_R2", "selection_mean_RMSE",
  "final_fit_R2", "final_fit_RMSE",
  sort(grep("^(R2|RMSE)\\.", names(all_results), value = TRUE)),
  "moran_observed_I", "moran_expected_I", "moran_p_value"
)
compact_cols <- compact_cols[compact_cols %in% names(all_results)]
compact <- all_results[, compact_cols, drop = FALSE]

md <- c(
  "# Random Forest group model summary",
  "",
  paste0("Run date: ", Sys.Date()),
  "",
  "Lat/Lon were excluded from model predictors and used only for spatial CV and Moran's I.",
  "The non-numeric biome column was excluded from RF predictors and used for biome-stratified CV when present.",
  "If CAST::knndm returned a random CV assignment, coordinate-cluster spatial folds were used and marked as spatial_cluster_10fold.",
  "Predictive R2 is computed as 1 - SSE/SST on cross-validated predictions.",
  "",
  "## Compact Summary",
  "",
  paste(capture.output(print(compact, row.names = FALSE)), collapse = "\n"),
  ""
)
writeLines(md, file.path(output_dir, "rf_group_model_summary.md"))

message("Done. Summary written to: ", file.path(output_dir, "rf_group_model_summary_wide.csv"))
