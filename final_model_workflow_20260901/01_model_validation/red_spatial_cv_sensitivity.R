#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ranger)
  library(spdep)
})

set.seed(20260616)

base_dir <- "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F4/map_prediction"
summary_file <- "model_outputs_red_optimized/red_optimized_best_model_summary.csv"
output_dir <- "model_outputs_red_optimized"

cfg <- list(seed = 20260616, trees = 500, mem_count = 5, moran_k = 8)

rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
pred_r2 <- function(obs, pred) {
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  denom <- sum((obs - mean(obs))^2)
  if (denom <= 0) return(NA_real_)
  1 - sum((obs - pred)^2) / denom
}
metric <- function(obs, pred) data.frame(R2 = pred_r2(obs, pred), RMSE = rmse(obs, pred), MAE = mean(abs(obs - pred), na.rm = TRUE))

make_split_list <- function(test_list, n, prefix) {
  lapply(seq_along(test_list), function(i) {
    test <- sort(unique(test_list[[i]]))
    list(name = paste0(prefix, i), train = setdiff(seq_len(n), test), test = test)
  })
}

spatial_kmeans <- function(coords, k) {
  set.seed(cfg$seed)
  km <- kmeans(scale(coords), centers = k, nstart = 100)
  make_split_list(split(seq_len(nrow(coords)), km$cluster), nrow(coords), paste0("K", k, "_"))
}

spatial_grid <- function(coords, nx, ny) {
  lat_bin <- cut(coords$Lat, breaks = nx, include.lowest = TRUE, labels = FALSE)
  lon_bin <- cut(coords$Lon, breaks = ny, include.lowest = TRUE, labels = FALSE)
  block <- paste(lat_bin, lon_bin, sep = "_")
  make_split_list(split(seq_len(nrow(coords)), block), nrow(coords), paste0("Grid", nx, "x", ny, "_"))
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

fit_predict_once <- function(train_x, train_y, test_x, model_kind, mtry) {
  mtry <- max(1, min(ncol(train_x), mtry))
  if (model_kind == "RF_raw") {
    fit <- ranger(.y ~ ., data = data.frame(.y = train_y, train_x, check.names = FALSE), num.trees = cfg$trees, mtry = mtry, min.node.size = 5, seed = cfg$seed)
    return(as.numeric(predict(fit, data = test_x)$predictions))
  }
  baseline <- min(train_y, na.rm = TRUE)
  shifted <- train_y - baseline
  if (model_kind == "RF_shift_log1p") {
    fit <- ranger(.y ~ ., data = data.frame(.y = log1p(shifted), train_x, check.names = FALSE), num.trees = cfg$trees, mtry = mtry, min.node.size = 5, seed = cfg$seed)
    return(as.numeric(expm1(predict(fit, data = test_x)$predictions) + baseline))
  }
  if (model_kind == "RF_hurdle_shift_log1p") {
    tol <- sqrt(.Machine$double.eps)
    present <- factor(ifelse(train_y > baseline + tol, "present", "absent"), levels = c("absent", "present"))
    class_fit <- ranger(.y ~ ., data = data.frame(.y = present, train_x, check.names = FALSE), num.trees = cfg$trees, mtry = mtry, min.node.size = 5, probability = TRUE, seed = cfg$seed)
    prob <- predict(class_fit, data = test_x)$predictions[, "present"]
    pos <- train_y > baseline + tol
    if (sum(pos) < 8) {
      excess <- rep(mean(shifted[pos], na.rm = TRUE), nrow(test_x))
    } else {
      reg_fit <- ranger(.y ~ ., data = data.frame(.y = log1p(shifted[pos]), train_x[pos, , drop = FALSE], check.names = FALSE), num.trees = cfg$trees, mtry = mtry, min.node.size = 3, seed = cfg$seed + 1)
      excess <- expm1(predict(reg_fit, data = test_x)$predictions)
    }
    return(as.numeric(baseline + prob * excess))
  }
  stop("Unknown model kind")
}

cv_predict <- function(x, y, splits, model_kind, mtry) {
  rows <- vector("list", length(splits))
  for (i in seq_along(splits)) {
    sp <- splits[[i]]
    pred <- fit_predict_once(x[sp$train, , drop = FALSE], y[sp$train], x[sp$test, , drop = FALSE], model_kind, mtry)
    rows[[i]] <- data.frame(row_id = sp$test, obs = y[sp$test], pred = pred)
  }
  do.call(rbind, rows)
}

build_x <- function(dat, environmental_factors, feature_mode) {
  vars <- trimws(strsplit(environmental_factors, ";", fixed = TRUE)[[1]])
  x <- dat[, vars, drop = FALSE]
  names(x) <- make.names(names(x), unique = TRUE)
  if (grepl("biome", feature_mode)) {
    mm <- as.data.frame(model.matrix(~ factor(dat$biome) - 1))
    names(mm) <- make.names(gsub("^factor\\(dat\\$biome\\)", "biome_", names(mm)), unique = TRUE)
    x <- cbind(x, mm)
  }
  if (grepl("coord", feature_mode)) {
    coord <- data.frame(coord_Lat = dat$Lat, coord_Lon = dat$Lon, coord_Lat2 = dat$Lat^2, coord_Lon2 = dat$Lon^2, coord_LatLon = dat$Lat * dat$Lon)
    x <- cbind(x, coord)
  }
  if (grepl("MEM", feature_mode)) {
    x <- cbind(x, make_mem(dat[, c("Lat", "Lon")], k = cfg$moran_k, n_mem = cfg$mem_count))
  }
  as.data.frame(x, check.names = FALSE)
}

summary <- read.csv(summary_file, check.names = FALSE)
rows <- list()
rid <- 1
for (i in seq_len(nrow(summary))) {
  group <- summary$group[i]
  dat <- read.csv(file.path(base_dir, paste0("env.meta.", group, "_scaled.csv")), check.names = FALSE)
  target <- summary$target[i]
  y <- dat[[target]]
  x <- build_x(dat, summary$environmental_factors[i], summary$feature_mode[i])
  coords <- dat[, c("Lat", "Lon")]
  methods <- list(
    spatial_kmeans_5fold = spatial_kmeans(coords, 5),
    spatial_kmeans_8fold = spatial_kmeans(coords, 8),
    spatial_kmeans_10fold = spatial_kmeans(coords, 10),
    spatial_grid_2x2 = spatial_grid(coords, 2, 2),
    spatial_grid_3x3 = spatial_grid(coords, 3, 3)
  )
  for (method_name in names(methods)) {
    pred <- cv_predict(x, y, methods[[method_name]], summary$model_type[i], summary$mtry[i])
    met <- metric(pred$obs, pred$pred)
    rows[[rid]] <- data.frame(group = group, target = target, model_type = summary$model_type[i], feature_mode = summary$feature_mode[i], spatial_cv_method = method_name, R2 = met$R2, RMSE = met$RMSE, MAE = met$MAE, stringsAsFactors = FALSE)
    rid <- rid + 1
  }
}

out <- do.call(rbind, rows)
write.csv(out, file.path(output_dir, "red_final_models_spatial_cv_alternatives.csv"), row.names = FALSE)
rounded <- out
num <- vapply(rounded, is.numeric, logical(1))
rounded[num] <- lapply(rounded[num], function(x) round(x, 6))
write.csv(rounded, file.path(output_dir, "red_final_models_spatial_cv_alternatives_rounded.csv"), row.names = FALSE)
print(rounded[order(rounded$group, -rounded$R2), ], row.names = FALSE)
