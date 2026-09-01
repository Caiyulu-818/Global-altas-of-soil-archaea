#!/usr/bin/env Rscript

# NFR distribution across soil ecosystems
# Reproduces the visual logic of panel d using the user's fr data.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

# -------------------------------------------------------------------------
# 1. Read and merge data
# -------------------------------------------------------------------------
sample_latest <- read_xlsx(
  "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/sampleinfo/sample.paddy.soil_filter1.xlsx"
)

fr <- read.csv(
  "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F4/FR/fr.csv",
  check.names = FALSE
)

fr <- fr %>%
  left_join(sample_latest, by = c("sample" = "Sample"))

if (!"Ecosystem" %in% names(fr)) {
  stop("The merged data do not contain an 'Ecosystem' column.")
}
if (!"nfr" %in% names(fr)) {
  stop("The fr table does not contain an 'nfr' column.")
}

# Keep the original ecosystem categories. Do not use Ecosystem1 here because
# that column collapses all non-agricultural ecosystems into one category.
ecosystem_order <- c(
  "Tundra",
  "Forest",
  "Grassland",
  "Agricultural Land",
  "Artificial Surfaces",
  "Wetland",
  "Paddy Soil"
)

plot_df <- fr %>%
  transmute(
    Ecosystem = as.character(Ecosystem),
    nfr = as.numeric(nfr)
  ) %>%
  filter(!is.na(Ecosystem), is.finite(nfr)) %>%
  mutate(
    Ecosystem = recode(
      Ecosystem,
      "Paddy soil" = "Paddy Soil",
      "Agricultural land" = "Agricultural Land",
      "Artificial surface" = "Artificial Surfaces",
      .default = Ecosystem
    ),
    Ecosystem = factor(Ecosystem, levels = ecosystem_order)
  ) %>%
  filter(!is.na(Ecosystem))

if (nrow(plot_df) == 0) {
  stop("No complete nfr observations remain after filtering.")
}

# -------------------------------------------------------------------------
# 2. Summary statistics and global Kruskal-Wallis test
# -------------------------------------------------------------------------
group_summary <- plot_df %>%
  group_by(Ecosystem) %>%
  summarise(
    n = n(),
    mean = mean(nfr),
    sd = sd(nfr),
    se = sd / sqrt(n),
    lower = mean - qt(0.975, df = max(n - 1, 1)) * se,
    upper = mean + qt(0.975, df = max(n - 1, 1)) * se,
    .groups = "drop"
  )

kw <- kruskal.test(nfr ~ Ecosystem, data = plot_df)
kw_label <- paste0(
  "Kruskal–Wallis, ",
  format.pval(kw$p.value, digits = 3, eps = 2.2e-16)
)

# Add sample sizes to x-axis labels.
axis_labels <- group_summary %>%
  mutate(label = paste0(as.character(Ecosystem), "\n(n = ", n, ")")) %>%
  { setNames(.$label, as.character(.$Ecosystem)) }

# Ecosystem-specific colours, reused consistently in later figures.
ecosystem_colours <- c(
  "Tundra" = "#8dd3c7",
  "Forest" = "#6ba76b",
  "Grassland" = "#1b9e57",
  "Agricultural Land" = "#d94f5c",
  "Artificial Surfaces" = "#4c78a8",
  "Wetland" = "#5f665f",
  "Paddy Soil" = "#e8b08f"
)

# -------------------------------------------------------------------------
# 3. Plot
# -------------------------------------------------------------------------
y_range <- range(plot_df$nfr, na.rm = TRUE)
y_span <- diff(y_range)
if (y_span == 0) y_span <- max(abs(y_range[1]), 1)

p_nfr <- ggplot(plot_df, aes(x = Ecosystem, y = nfr, fill = Ecosystem, colour = Ecosystem)) +
  # Distribution shape
  geom_violin(
    width = 0.85,
    alpha = 0.28,
    linewidth = 0.35,
    trim = FALSE,
    colour = NA
  ) +
  # Individual observations
  geom_jitter(
    width = 0.10,
    height = 0,
    size = 1.25,
    alpha = 0.25,
    stroke = 0
  ) +
  # Compact boxplot over the distribution
  geom_boxplot(
    width = 0.14,
    outlier.shape = NA,
    linewidth = 0.45,
    fill = "white",
    colour = "grey20",
    alpha = 0.85
  ) +
  # Mean and 95% confidence interval
  geom_errorbar(
    data = group_summary,
    aes(x = Ecosystem, ymin = lower, ymax = upper),
    inherit.aes = FALSE,
    width = 0.06,
    linewidth = 0.45,
    colour = "black"
  ) +
  geom_point(
    data = group_summary,
    aes(x = Ecosystem, y = mean),
    inherit.aes = FALSE,
    shape = 23,
    size = 2.5,
    fill = "white",
    colour = "black",
    stroke = 0.55
  ) +
  annotate(
    "text",
    x = 1,
    y = max(plot_df$nfr, na.rm = TRUE) + 0.09 * y_span,
    label = kw_label,
    hjust = 0,
    vjust = 0,
    size = 4.1,
    fontface = "italic"
  ) +
  scale_fill_manual(values = ecosystem_colours, drop = FALSE) +
  scale_colour_manual(values = ecosystem_colours, drop = FALSE) +
  scale_x_discrete(labels = axis_labels, drop = FALSE) +
  scale_y_continuous(
    name = "Normalized functional redundancy (NFR)",
    expand = expansion(mult = c(0.04, 0.16)),
    labels = label_number(accuracy = 0.1)
  ) +
  labs(x = NULL) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 28,
      hjust = 1,
      vjust = 1,
      colour = "grey15",
      size = 10.5,
      lineheight = 0.9
    ),
    axis.text.y = element_text(colour = "grey15"),
    axis.title.y = element_text(face = "bold", margin = margin(r = 9)),
    axis.line = element_line(linewidth = 0.45, colour = "grey20"),
    axis.ticks = element_line(linewidth = 0.4, colour = "grey20"),
    plot.margin = margin(8, 10, 18, 8),
    panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.major.x = element_blank()
  )

print(p_nfr)

ggsave(
  filename = "nfr_by_ecosystem.pdf",
  plot = p_nfr,
  width = 8.2,
  height = 5.8,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = "nfr_by_ecosystem.png",
  plot = p_nfr,
  width = 8.2,
  height = 5.8,
  units = "in",
  dpi = 600,
  bg = "white"
)
