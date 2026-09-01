#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(caret)
  library(ranger)
  library(spdep)
  library(CAST)
})

set.seed(20260616)

args <- commandArgs(trailingOnly = TRUE)
guild_number <- if (length(args) > 0) as.integer(gsub("[^0-9]", "", args[1])) else 3
if (is.na(guild_number) || !guild_number %in% 1:5) stop("Use red1, red2, red3, red4, or red5 as the optional argument.")

red_scaled_file <- sprintf("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F4/map_prediction/env.meta.red%d_scaled.csv", guild_number)
red_target <- paste0("Guild.", guild_number)
abundance_target <- paste0("Guild ", guild_number)
abundance_file <- "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F4/map_prediction/env.abundance_total.csv"
environment_file <- "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/paper/写作参考/nature/submit/Supplementary Table/Supplementary Table 10.csv"
output_dir <- sprintf("model_outputs_red%d_pa", guild_number)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cfg <- list(
  seed = 20260616,
  trees = 1000,
  screen_trees = 400,
  k_folds = 10,
  missing_threshold = 0.20,
  correlation_cutoff = 0.90,
  mem_count = 5,
  moran_k = 8
)

mode_value <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

auc_roc <- function(obs, prob) {
  ok <- is.finite(obs) & is.finite(prob)
  obs <- as.integer(obs[ok])
  prob <- prob[ok]
  n1 <- sum(obs == 1)
  n0 <- sum(obs == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  ranks <- rank(prob, ties.method = "average")
  (sum(ranks[obs == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

tss_best <- function(obs, prob) {
  ok <- is.finite(obs) & is.finite(prob)
  obs <- as.integer(obs[ok])
  prob <- prob[ok]
  if (length(unique(obs)) < 2) {
    return(data.frame(TSS = NA_real_, threshold = NA_real_, sensitivity = NA_real_, specificity = NA_real_))
  }
  thresholds <- sort(unique(prob))
  if (length(thresholds) > 300) {
    thresholds <- unique(as.numeric(quantile(prob, probs = seq(0, 1, length.out = 300), na.rm = TRUE)))
  }
  rows <- lapply(thresholds, function(th) {
    pred <- as.integer(prob >= th)
    tp <- sum(pred == 1 & obs == 1)
    tn <- sum(pred == 0 & obs == 0)
    fp <- sum(pred == 1 & obs == 0)
    fn <- sum(pred == 0 & obs == 1)
    sens <- ifelse((tp + fn) == 0, NA_real_, tp / (tp + fn))
    spec <- ifelse((tn + fp) == 0, NA_real_, tn / (tn + fp))
    data.frame(TSS = sens + spec - 1, threshold = th, sensitivity = sens, specificity = spec)
  })
  out <- do.call(rbind, rows)
  out[order(-out$TSS, out$threshold), ][1, , drop = FALSE]
}

metric <- function(obs, prob) {
  tss <- tss_best(obs, prob)
  data.frame(
    AUC = auc_roc(obs, prob),
    TSS = tss$TSS,
    threshold = tss$threshold,
    sensitivity = tss$sensitivity,
    specificity = tss$specificity,
    Brier = mean((obs - prob)^2, na.rm = TRUE),
    prevalence = mean(obs == 1, na.rm = TRUE),
    n = length(obs),
    positives = sum(obs == 1),
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
  train_list <- createMultiFolds(factor(y), k = k, times = repeats)
  lapply(names(train_list), function(nm) {
    train <- sort(unique(train_list[[nm]]))
    list(name = nm, train = train, test = setdiff(seq_along(y), train))
  })
}

stratified_group_folds <- function(y, group, k = 10) {
  if (!is.null(group) && length(unique(group)) > 1) {
    strata <- interaction(ifelse(y == 1, "present", "absent"), group, drop = TRUE)
  } else {
    strata <- factor(y)
  }
  test_list <- createFolds(strata, k = k, list = TRUE, returnTrain = FALSE)
  make_split_list(test_list, length(y), "Strat")
}

loocv_folds <- function(n) {
  lapply(seq_len(n), function(i) list(name = paste0("LOO", i), train = setdiff(seq_len(n), i), test = i))
}

spatial_kmeans <- function(coords, k) {
  set.seed(cfg$seed)
  km <- kmeans(scale(coords), centers = k, nstart = 100)
  make_split_list(split(seq_len(nrow(coords)), km$cluster), nrow(coords), paste0("K", k, "_"))
}

spatial_grid <- function(coords, nx, ny) {
  lat_bin <- cut(coords$Lat, breaks = nx, include.lowest = TRUE, labels = FALSE)
  lon_bin <- cut(coords$Lon, breaks = ny, include.lowest = TRUE, labels = FALSE)
  make_split_list(split(seq_len(nrow(coords)), paste(lat_bin, lon_bin, sep = "_")), nrow(coords), paste0("Grid", nx, "x", ny, "_"))
}

knndm_feature <- function(coords, k = 10) {
  tryCatch({
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
    if (random_assignment) stop("kNNDM returned random assignment")
    list(
      splits = lapply(seq_along(knn$indx_test), function(i) {
        list(name = paste0("kNNDM", i), train = sort(unique(knn$indx_train[[i]])), test = sort(unique(knn$indx_test[[i]])))
      }),
      label = "kNNDM_10fold",
      note = paste0("CAST::knndm in coordinate feature space; W=", signif(knn$W, 6))
    )
  }, error = function(e) {
    list(splits = spatial_kmeans(coords, k), label = "spatial_kmeans_10fold", note = paste0("kNNDM fallback: ", conditionMessage(e)))
  })
}

make_mem <- function(coords, k = 8, n_mem = 5) {
  coords <- as.matrix(coords)
  if (any(duplicated(as.data.frame(coords)))) {
    set.seed(cfg$seed)
    coords <- coords + matrix(rnorm(length(coords), 0, 1e-9), ncol = 2)
  }
  knn <- knearneigh(coords, k = k, longlat = FALSE)
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
  lw <- nb2listw(knn2nb(knn), style = "W", zero.policy = TRUE)
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

dedupe_site <- function(d) d[!duplicated(d$site), , drop = FALSE]

aggregate_full_env <- function() {
  ab <- dedupe_site(read.csv(abundance_file, check.names = FALSE))
  env <- dedupe_site(read.csv(environment_file, check.names = FALSE))
  d <- merge(ab[, c("site", abundance_target), drop = FALSE], env, by = "site", sort = FALSE)
  d$coord_id <- paste(round(d$Lat, 6), round(d$Lon, 6), sep = ",")
  numeric_cols <- names(d)[vapply(d, is.numeric, logical(1))]
  numeric_env <- setdiff(names(env)[vapply(env, is.numeric, logical(1))], c("", "Lat", "Lon"))
  pieces <- split(d, d$coord_id)
  rows <- lapply(names(pieces), function(id) {
    z <- pieces[[id]]
    out <- as.data.frame(lapply(z[, c(numeric_env, "Lat", "Lon"), drop = FALSE], function(x) mean(x, na.rm = TRUE)), check.names = FALSE)
    out$presence <- as.integer(any(z[[abundance_target]] > 0, na.rm = TRUE))
    out$biome <- if ("Ecosystem" %in% names(z)) mode_value(z$Ecosystem) else NA_character_
    out$coord_id <- id
    out$n_samples_at_coord <- nrow(z)
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  list(data = out, numeric_env = numeric_env, label = sprintf("red%d_full_env_aggregated", guild_number))
}

prepare_red_scaled <- function() {
  d <- read.csv(red_scaled_file, check.names = FALSE)
  d$presence <- as.integer(d[[red_target]] > min(d[[red_target]]) + 1e-10)
  numeric_env <- setdiff(names(d)[vapply(d, is.numeric, logical(1))], c(red_target, "Lat", "Lon", "presence"))
  list(data = d, numeric_env = numeric_env, label = sprintf("red%d_scaled_limited", guild_number))
}

prepare_dataset <- function(src) {
  d <- src$data
  numeric_env <- src$numeric_env
  miss <- sapply(d[, numeric_env, drop = FALSE], function(x) mean(is.na(x)))
  numeric_env <- names(miss)[miss <= cfg$missing_threshold]
  x <- median_impute(d[, numeric_env, drop = FALSE])
  nzv <- nearZeroVar(x, saveMetrics = TRUE)
  x <- x[, rownames(nzv)[!nzv$nzv], drop = FALSE]
  if (ncol(x) > 1) {
    drop_cor <- findCorrelation(cor(x, use = "pairwise.complete.obs"), cutoff = cfg$correlation_cutoff, names = TRUE, exact = TRUE)
    x <- x[, setdiff(names(x), drop_cor), drop = FALSE]
  }
  safe_names <- make.names(names(x), unique = TRUE)
  name_map <- setNames(names(x), safe_names)
  names(x) <- safe_names
  list(
    label = src$label,
    data = d,
    x_env = x,
    name_map = name_map,
    y = as.integer(d$presence),
    coords = d[, c("Lat", "Lon"), drop = FALSE],
    biome = if ("biome" %in% names(d)) d$biome else rep("all", nrow(d))
  )
}

fit_rf <- function(train_x, train_y, test_x, model_kind) {
  train_y <- factor(ifelse(train_y == 1, "present", "absent"), levels = c("absent", "present"))
  if (length(unique(train_y)) < 2) {
    p <- mean(train_y == "present")
    return(rep(p, nrow(test_x)))
  }
  class_weights <- if (model_kind == "RF_balanced") {
    tab <- table(train_y)
    c(absent = 1, present = as.numeric(tab["absent"] / tab["present"]))
  } else {
    NULL
  }
  fit <- ranger(
    .y ~ .,
    data = data.frame(.y = train_y, train_x, check.names = FALSE),
    num.trees = cfg$trees,
    probability = TRUE,
    mtry = max(1, floor(sqrt(ncol(train_x)))),
    min.node.size = 3,
    class.weights = class_weights,
    seed = cfg$seed
  )
  as.numeric(predict(fit, data = test_x)$predictions[, "present"])
}

cv_predict <- function(x, y, splits, model_kind) {
  rows <- vector("list", length(splits))
  for (i in seq_along(splits)) {
    sp <- splits[[i]]
    prob <- fit_rf(x[sp$train, , drop = FALSE], y[sp$train], x[sp$test, , drop = FALSE], model_kind)
    rows[[i]] <- data.frame(row_id = sp$test, obs = y[sp$test], prob = prob)
  }
  do.call(rbind, rows)
}

full_prob <- function(x, y, model_kind) fit_rf(x, y, x, model_kind)

rank_vars <- function(x, y) {
  yy <- factor(ifelse(y == 1, "present", "absent"), levels = c("absent", "present"))
  tab <- table(yy)
  weights <- c(absent = 1, present = as.numeric(tab["absent"] / tab["present"]))
  fit <- ranger(
    .y ~ .,
    data = data.frame(.y = yy, x, check.names = FALSE),
    num.trees = cfg$trees,
    probability = TRUE,
    importance = "permutation",
    mtry = max(1, floor(sqrt(ncol(x)))),
    class.weights = weights,
    seed = cfg$seed
  )
  sort(fit$variable.importance, decreasing = TRUE)
}

build_x <- function(ds, vars, feature_mode) {
  x <- ds$x_env[, vars, drop = FALSE]
  if (grepl("biome", feature_mode)) {
    mm <- as.data.frame(model.matrix(~ factor(ds$biome) - 1))
    names(mm) <- make.names(gsub("^factor\\(ds\\$biome\\)", "biome_", names(mm)), unique = TRUE)
    x <- cbind(x, mm)
  }
  if (grepl("coord", feature_mode)) {
    coord <- data.frame(coord_Lat = ds$coords$Lat, coord_Lon = ds$coords$Lon, coord_Lat2 = ds$coords$Lat^2, coord_Lon2 = ds$coords$Lon^2, coord_LatLon = ds$coords$Lat * ds$coords$Lon)
    x <- cbind(x, coord)
  }
  if (grepl("MEM", feature_mode)) x <- cbind(x, make_mem(ds$coords, cfg$moran_k, cfg$mem_count))
  as.data.frame(x, check.names = FALSE)
}

score_candidate <- function(d) {
  d$score <- d$mean_AUC + 0.75 * d$best_spatial_AUC + 0.5 * d$mean_TSS + 0.5 * d$best_spatial_TSS
  d[order(-d$score, -d$best_spatial_AUC, -d$mean_AUC), ]
}

sources <- list(prepare_red_scaled(), aggregate_full_env())
all_best <- list()
all_final_cv <- list()

for (src in sources) {
  ds <- prepare_dataset(src)
  message("Dataset: ", ds$label, " n=", length(ds$y), " positives=", sum(ds$y))
  imp <- rank_vars(ds$x_env, ds$y)
  imp_df <- data.frame(variable_safe = names(imp), variable = ds$name_map[names(imp)], importance = as.numeric(imp), stringsAsFactors = FALSE)
  write.csv(imp_df, file.path(output_dir, paste0(ds$label, "_importance.csv")), row.names = FALSE)
  ranked <- names(imp)
  sizes <- unique(pmin(length(ranked), c(2, 3, 5, 8, 12, 16, 20, 25, 30)))
  sizes <- sizes[sizes >= 1]
  feature_modes <- c("env_only", "env_plus_biome", "env_plus_MEM", "env_plus_biome_MEM", "env_plus_coord", "env_plus_biome_coord_MEM")
  cv_methods <- list(
    random_10fold = random_folds(ds$y, cfg$k_folds, repeats = 1),
    biome_stratified_10fold = stratified_group_folds(ds$y, ds$biome, cfg$k_folds),
    spatial_kmeans_10fold = spatial_kmeans(ds$coords, cfg$k_folds),
    spatial_grid_3x3 = spatial_grid(ds$coords, 3, 3)
  )
  knn <- knndm_feature(ds$coords, cfg$k_folds)
  cv_methods[[knn$label]] <- knn$splits

  candidate_rows <- list()
  candidate_cv <- list()
  cid <- 1
  for (model_kind in c("RF_balanced", "RF_unweighted")) {
    for (feature_mode in feature_modes) {
      for (sz in sizes) {
        vars <- ranked[seq_len(sz)]
        x <- build_x(ds, vars, feature_mode)
        cv_rows <- list()
        for (method_name in names(cv_methods)) {
          pred <- cv_predict(x, ds$y, cv_methods[[method_name]], model_kind)
          met <- metric(pred$obs, pred$prob)
          cv_rows[[method_name]] <- data.frame(dataset = ds$label, candidate_key = paste0("cand", cid), model_kind = model_kind, feature_mode = feature_mode, n_env_factors = length(vars), environmental_factors = paste(ds$name_map[vars], collapse = "; "), cv_method = method_name, AUC = met$AUC, TSS = met$TSS, threshold = met$threshold, sensitivity = met$sensitivity, specificity = met$specificity, Brier = met$Brier, positives = met$positives, stringsAsFactors = FALSE)
        }
        cv_eval <- do.call(rbind, cv_rows)
        prob_full <- full_prob(x, ds$y, model_kind)
        moran <- calc_moran(ds$y - prob_full, ds$coords, cfg$moran_k)
        spatial <- cv_eval[cv_eval$cv_method %in% c("kNNDM_10fold", "spatial_kmeans_10fold", "spatial_grid_3x3"), ]
        best_spatial <- spatial[order(-spatial$AUC, -spatial$TSS), ][1, ]
        candidate_rows[[cid]] <- data.frame(
          dataset = ds$label,
          candidate_key = paste0("cand", cid),
          model_kind = model_kind,
          feature_mode = feature_mode,
          n_env_factors = length(vars),
          environmental_factors = paste(ds$name_map[vars], collapse = "; "),
          spatial_terms = paste(c(
            if (grepl("coord", feature_mode)) c("coord_Lat", "coord_Lon", "coord_Lat2", "coord_Lon2", "coord_LatLon") else character(0),
            if (grepl("MEM", feature_mode)) paste0("MEM", seq_len(cfg$mem_count)) else character(0)
          ), collapse = "; "),
          mean_AUC = mean(cv_eval$AUC, na.rm = TRUE),
          mean_TSS = mean(cv_eval$TSS, na.rm = TRUE),
          best_spatial_method = best_spatial$cv_method,
          best_spatial_AUC = best_spatial$AUC,
          best_spatial_TSS = best_spatial$TSS,
          best_spatial_Brier = best_spatial$Brier,
          moran_observed_I = moran$moran_observed_I,
          moran_expected_I = moran$moran_expected_I,
          moran_p_value = moran$moran_p_value,
          moran_pass = moran$moran_p_value > 0.05,
          stringsAsFactors = FALSE
        )
        candidate_cv[[cid]] <- cv_eval
        cid <- cid + 1
      }
    }
  }
  candidates <- score_candidate(do.call(rbind, candidate_rows))
  write.csv(candidates, file.path(output_dir, paste0(ds$label, "_candidate_summary.csv")), row.names = FALSE)
  write.csv(do.call(rbind, candidate_cv), file.path(output_dir, paste0(ds$label, "_candidate_cv_long.csv")), row.names = FALSE)

  pass <- candidates[candidates$moran_pass & candidates$best_spatial_AUC >= 0.6 & candidates$best_spatial_TSS > 0, ]
  if (nrow(pass) == 0) pass <- candidates[candidates$moran_pass, ]
  if (nrow(pass) == 0) pass <- candidates
  chosen <- pass[order(-pass$score, -pass$best_spatial_AUC, pass$n_env_factors, -pass$mean_AUC), ][1, ]
  vars_chosen <- names(ds$name_map)[match(trimws(strsplit(chosen$environmental_factors, ";", fixed = TRUE)[[1]]), ds$name_map)]
  x_chosen <- build_x(ds, vars_chosen, chosen$feature_mode)
  final_methods <- list(
    LOOCV = loocv_folds(length(ds$y)),
    random_10fold_x5 = random_folds(ds$y, cfg$k_folds, repeats = 5),
    biome_stratified_10fold = stratified_group_folds(ds$y, ds$biome, cfg$k_folds),
    spatial_kmeans_10fold = spatial_kmeans(ds$coords, cfg$k_folds),
    spatial_grid_3x3 = spatial_grid(ds$coords, 3, 3)
  )
  final_knn <- knndm_feature(ds$coords, cfg$k_folds)
  final_methods[[final_knn$label]] <- final_knn$splits
  final_cv_rows <- list()
  for (method_name in names(final_methods)) {
    pred <- cv_predict(x_chosen, ds$y, final_methods[[method_name]], chosen$model_kind)
    met <- metric(pred$obs, pred$prob)
    final_cv_rows[[method_name]] <- data.frame(dataset = ds$label, model_kind = chosen$model_kind, feature_mode = chosen$feature_mode, cv_method = method_name, AUC = met$AUC, TSS = met$TSS, threshold = met$threshold, sensitivity = met$sensitivity, specificity = met$specificity, Brier = met$Brier, positives = met$positives, stringsAsFactors = FALSE)
  }
  final_cv <- do.call(rbind, final_cv_rows)
  prob_full <- full_prob(x_chosen, ds$y, chosen$model_kind)
  moran <- calc_moran(ds$y - prob_full, ds$coords, cfg$moran_k)
  spatial_final <- final_cv[final_cv$cv_method %in% c("kNNDM_10fold", "spatial_kmeans_10fold", "spatial_grid_3x3"), ]
  best_spatial_final <- spatial_final[order(-spatial_final$AUC, -spatial_final$TSS), ][1, ]
  best <- data.frame(
    dataset = ds$label,
    n = length(ds$y),
    positives = sum(ds$y),
    prevalence = mean(ds$y),
    model_kind = chosen$model_kind,
    feature_mode = chosen$feature_mode,
    environmental_factors = chosen$environmental_factors,
    spatial_terms = chosen$spatial_terms,
    best_spatial_method = best_spatial_final$cv_method,
    best_spatial_AUC = best_spatial_final$AUC,
    best_spatial_TSS = best_spatial_final$TSS,
    best_spatial_sensitivity = best_spatial_final$sensitivity,
    best_spatial_specificity = best_spatial_final$specificity,
    LOOCV_AUC = final_cv$AUC[final_cv$cv_method == "LOOCV"][1],
    LOOCV_TSS = final_cv$TSS[final_cv$cv_method == "LOOCV"][1],
    random_10fold_AUC = final_cv$AUC[final_cv$cv_method == "random_10fold_x5"][1],
    random_10fold_TSS = final_cv$TSS[final_cv$cv_method == "random_10fold_x5"][1],
    biome_10fold_AUC = final_cv$AUC[final_cv$cv_method == "biome_stratified_10fold"][1],
    biome_10fold_TSS = final_cv$TSS[final_cv$cv_method == "biome_stratified_10fold"][1],
    kNNDM_AUC = final_cv$AUC[final_cv$cv_method == "kNNDM_10fold"][1],
    kNNDM_TSS = final_cv$TSS[final_cv$cv_method == "kNNDM_10fold"][1],
    moran_observed_I = moran$moran_observed_I,
    moran_expected_I = moran$moran_expected_I,
    moran_p_value = moran$moran_p_value,
    moran_pass = moran$moran_p_value > 0.05,
    stringsAsFactors = FALSE
  )
  write.csv(final_cv, file.path(output_dir, paste0(ds$label, "_final_cv.csv")), row.names = FALSE)
  all_best[[ds$label]] <- best
  all_final_cv[[ds$label]] <- final_cv
}

best_summary <- do.call(rbind, all_best)
write.csv(best_summary, file.path(output_dir, sprintf("red%d_pa_best_model_summary.csv", guild_number)), row.names = FALSE)
rounded <- best_summary
num <- vapply(rounded, is.numeric, logical(1))
rounded[num] <- lapply(rounded[num], function(x) round(x, 6))
write.csv(rounded, file.path(output_dir, sprintf("red%d_pa_best_model_summary_rounded.csv", guild_number)), row.names = FALSE)
print(rounded, row.names = FALSE)
