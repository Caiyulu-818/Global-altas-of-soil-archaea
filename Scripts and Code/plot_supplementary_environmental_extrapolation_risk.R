#!/usr/bin/env Rscript

# Supplementary figure: continuous environmental extrapolation risk for G1-G5
# Input: 0.01-degree environmental_extrapolation_risk_ratio.tif layers
# Display: 0.1-degree mean aggregation for efficient, readable visualization

suppressPackageStartupMessages({
  library(terra)
  library(tidyterra)
  library(ggplot2)
  library(sf)
  library(dplyr)
  library(rnaturalearth)
  library(patchwork)
  library(scales)
})

# -------------------------------------------------------------------------
# 1. Input and output paths
# -------------------------------------------------------------------------
risk_files <- c(
  "Guild 1" = "/Volumes/DEEP/extra/Guild/risk/g1chi90environmental_extrapolation_risk_ratio_abundance.tif",
  "Guild 2" = "/Volumes/DEEP/extra/Guild/risk/g2environmental_extrapolation_risk_ratio.tif",
  "Guild 3" = "/Volumes/DEEP/extra/Guild/risk/g3environmental_extrapolation_risk_ratio_abundance.tif",
  "Guild 4" = "/Volumes/DEEP/extra/Guild/risk/g4environmental_extrapolation_risk_ratio.tif",
  "Guild 5" = "/Volumes/DEEP/extra/Guild/risk/g5environmental_extrapolation_risk_ratio.tif"
)

out_dir <- "/Volumes/DEEP/extra/Guild/supplementary_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

missing_files <- risk_files[!file.exists(risk_files)]
if (length(missing_files) > 0) {
  stop("Missing risk-ratio TIFF(s):\n", paste(missing_files, collapse = "\n"))
}

# -------------------------------------------------------------------------
# 2. Shared map settings
# -------------------------------------------------------------------------
# Ratio = ridge-Mahalanobis distance / chi-square 0.90 threshold.
# A ratio of 1 marks the training-environment boundary.
display_limit <- 3
aggregation_factor <- 10
robinson_crs <- "+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs"

world_robinson <- ne_countries(scale = "medium", returnclass = "sf") |>
  filter(continent != "Antarctica") |>
  st_transform(robinson_crs)

clean_risk_ratio <- function(r) {
  # Retain valid zero values while converting common numeric NoData values to NA.
  ifel(r <= -9990, NA, r)
}

make_risk_map <- function(file_path, guild_name) {
  message("Preparing ", guild_name)

  r <- rast(file_path)
  r <- clean_risk_ratio(r)
  names(r) <- "risk_ratio"

  # The full-resolution TIFF remains the data product; aggregation is only for
  # visualization and prevents a multi-gigabyte PDF.
  r_display <- aggregate(r, fact = aggregation_factor, fun = "mean", na.rm = TRUE)
  r_display <- crop(r_display, ext(-180, 180, -60, 85))
  r_robinson <- project(r_display, robinson_crs, method = "bilinear")
  names(r_robinson) <- "risk_ratio"

  ggplot() +
    geom_sf(
      data = world_robinson,
      fill = "grey92",
      colour = "grey85",
      linewidth = 0.08
    ) +
    geom_spatraster(
      data = r_robinson,
      aes(fill = risk_ratio),
      maxcell = Inf
    ) +
    viridis::scale_fill_viridis(
      option = "inferno",
      begin = 0.02,
      end = 0.98,
      limits = c(0, display_limit),
      oob = scales::squish,
      breaks = c(0, 0.5, 1, 2, 3),
      labels = c("0", "0.5", "1 (chi90)", "2", ">=3"),
      name = "Environmental\nextrapolation risk ratio",
      na.value = "transparent"
    ) +
    coord_sf(crs = robinson_crs, datum = NA, expand = FALSE) +
    labs(title = guild_name) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.position = "bottom",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      legend.key.width = grid::unit(1.8, "cm"),
      legend.key.height = grid::unit(0.35, "cm"),
      plot.margin = margin(2, 2, 2, 2, unit = "pt")
    )
}

# -------------------------------------------------------------------------
# 3. Create a common-scale five-panel figure
# -------------------------------------------------------------------------
plots <- Map(make_risk_map, risk_files, names(risk_files))

figure <- wrap_plots(plots, ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")

pdf_file <- file.path(out_dir, "Supplementary_Fig_environmental_extrapolation_risk_G1_G5chi90.pdf")
png_file <- file.path(out_dir, "Supplementary_Fig_environmental_extrapolation_risk_G1_G5chi90.png")
 
ggsave(pdf_file, figure, width = 14, height = 11, units = "in", dpi = 450)
ggsave(png_file, figure, width = 14, height = 11, units = "in", dpi = 450, bg = "white")

message("Finished. Outputs:")
message("  ", pdf_file)
message("  ", png_file)
