setwd("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/chi90")
# ==========================================
# 模块八：基于外推安全网格的不确定性 (SD & CV) 掩膜与可视化
# ==========================================
cat("🚀 正在启动不确定性 (SD & CV) 的掩膜处理与可视化...\n")

library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rnaturalearth)
library(stringr)
library(patchwork)

# ==========================================
# 1. 基础配置与统一掩膜提取函数
# ==========================================
world_map <- ne_countries(scale = "medium", returnclass = "sf") %>% filter(admin != "Antarctica")

# 你的 TIF 文件夹根目录
base_mask_path <- "/Volumes/DEEP/extra/"
strictness <- "chi90"

# 统一读取、降采样并掩膜的核心函数
process_uncertainty <- function(uncert_path, var_type, scenario_folder, scenario_name) {
  cat("⏳ 正在处理:", var_type, "-", scenario_name, "...\n")
  
  # 1. 读取原始不确定性图层 (SD或CV) 并降采样
  r_uncert <- rast(uncert_path)
  r_uncert_agg <- aggregate(r_uncert, fact = 10, fun = "mean", na.rm = TRUE)
  
  # 2. 读取对应的 strict 预测结果作为 Mask (掩膜)
  # 注意：这里我们利用之前的 prediction 栅格，只要 prediction 是 NA，就屏蔽不确定性
  mask_path <- file.path(base_mask_path, var_type, scenario_folder, 
                         paste0("prediction_ridge_reliable_", strictness, ".tif"))
  r_mask <- rast(mask_path)
  r_mask_agg <- aggregate(r_mask, fact = 10, fun = "mean", na.rm = TRUE)
  
  # 3. 核心步骤：应用掩膜 (Masking)
  r_masked <- mask(r_uncert_agg, r_mask_agg)
  
  # 4. 转为 DataFrame
  df <- as.data.frame(r_masked, xy = TRUE, na.rm = TRUE)
  colnames(df)[3] <- "Value"
  df$Scenario <- scenario_name
  df$Variable <- var_type
  
  return(df)
}


# ==========================================
# 2. 批量加载并掩膜所有数据 (TD & GFM, SD & CV)
# ==========================================
cat("📂 正在构建全维度不确定性数据流...\n")

# --- A. 路径字典配置 ---
paths_td_sd <- c(
  "2020_SSP245" = "/Volumes/DEEP/tiff/uncertainty/outlier/TD_sd_2020SSP245.tif",
  "2100_SSP126" = "/Volumes/DEEP/tiff/uncertainty/outlier/TD_sd_2100SSP126.tif",
  "2100_SSP245" = "/Volumes/DEEP/tiff/uncertainty/outlier/1/TD_sd_2100SSP245.tif",
  "2100_SSP370" = "/Volumes/DEEP/tiff/uncertainty/outlier/1/TD_sd_2100ssp370.tif"
)


paths_td_cv <- c(
  "2020_SSP245" = "/Volumes/DEEP/tiff/uncertainty/cv/TD_uncertainty_cv_2020SSP245.tif",
  "2100_SSP126" = "/Volumes/DEEP/tiff/uncertainty/cv/TD_uncertainty2100SSP126.tif", # 注意你这个文件命名少了个_cv
  "2100_SSP245" = "/Volumes/DEEP/tiff/uncertainty/cv/TD_uncertainty_cv_2100SSP245.tif",
  "2100_SSP370" = "/Volumes/DEEP/tiff/uncertainty/cv/TD_uncertainty_cv_2100ssp370.tif"
)


# 强制排序的因子 Level
scen_levels <- c("2020 SSP2-4.5", "2100 SSP1-2.6", "2100 SSP2-4.5", "2100 SSP3-7.0")

# --- B. 批量提取与掩膜执行 ---
cat("⏳ 提取 TD SD...\n")
df_td_sd <- bind_rows(
  process_uncertainty(paths_td_sd["2020_SSP245"], "TD", "cropland2020_SSP245", "2020 SSP2-4.5"),
  process_uncertainty(paths_td_sd["2100_SSP126"], "TD", "cropland2100_SSP126", "2100 SSP1-2.6"),
  process_uncertainty(paths_td_sd["2100_SSP245"], "TD", "cropland2100_SSP245", "2100 SSP2-4.5"),
  process_uncertainty(paths_td_sd["2100_SSP370"], "TD", "cropland2100_SSP370", "2100 SSP3-7.0")
) %>% mutate(Scenario = factor(Scenario, levels = scen_levels))


cat("⏳ 提取 TD CV...\n")
df_td_cv <- bind_rows(
  process_uncertainty(paths_td_cv["2020_SSP245"], "TD", "cropland2020_SSP245", "2020 SSP2-4.5"),
  process_uncertainty(paths_td_cv["2100_SSP126"], "TD", "cropland2100_SSP126", "2100 SSP1-2.6"),
  process_uncertainty(paths_td_cv["2100_SSP245"], "TD", "cropland2100_SSP245", "2100 SSP2-4.5"),
  process_uncertainty(paths_td_cv["2100_SSP370"], "TD", "cropland2100_SSP370", "2100 SSP3-7.0")
) %>% mutate(Scenario = factor(Scenario, levels = scen_levels))



# ==========================================
# 3. 高级出图引擎 1：离散型 SD 面板图 (Discrete Binned Maps)
# ==========================================
# 自定义极其护眼的绿色色阶
green_palette <- c("#E3EBC5", "#B9D48F", "#9CCC8C", "#69B37B", "#239B56", "#006D2C")

plot_sd_atlas <- function(data_df, title_str) {
  
  # 自动计算当前指标的全局最大值，向上取整 (例如 0.087 -> 0.09)
  sd_max <- ceiling(max(data_df$Value, na.rm = TRUE) * 100) / 100
  breaks_val <- seq(0, sd_max, length.out = 7)
  
  p <- ggplot(data_df) +
    geom_sf(data = world_map, fill = "grey90", color = "white", linewidth = 0.1) +
    
    # ✨ 物理填缝魔法，拒绝百叶窗白纹
    geom_tile(aes(x = x, y = y, fill = Value), width = 0.12, height = 0.12) +
    
    # 使用连续映射并强制离散化 (等价于你之前的 classify 切片)
    scale_fill_stepsn(
      colors = green_palette,
      breaks = breaks_val,
      limits = c(0, sd_max),
      na.value = "transparent",
      name = "Standard Deviation (SD)"
    ) +
    
    # ✨ 自动 2x2 情景分面
    facet_wrap(~ Scenario, ncol = 2) +
    
    labs(title = title_str) +
    theme_void() +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 15, t = 10)),
      strip.text = element_text(size = 14, face = "bold", margin = ggplot2::margin(b = 8, t = 8)),
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold", vjust = 1),
      legend.text = element_text(size = 10, face = "bold")
    ) +
    coord_sf(ylim = c(-60, 90), expand = FALSE) +
    # 横向拉长图例
    guides(fill = guide_colorsteps(barwidth = 25, barheight = 0.8, title.position = "top", title.hjust = 0.5))
  
  return(p)
}

# ==========================================
# 4. 升级版连续型 CV 面板图 (使用顶刊蓝-黄-红渐变色卡)
# ==========================================

# 精准复刻参考图中的 Blue-Yellow-Red 连续色盘
cv_palette <- c(
  "#0F5B9E", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", # 蓝调区 (低 CV)
  "#FFFFBF",                                             # 浅黄过渡 (中等 CV)
  "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"  # 红橘区 (高 CV)
)

plot_cv_atlas <- function(data_df, title_str) {
  
  # 自动计算 CV 最大值，向上取整保留两位小数 (保证色阶统一)
  cv_max <- ceiling(max(data_df$Value, na.rm = TRUE) * 100) / 100
  
  p <- ggplot(data_df) +
    geom_sf(data = world_map, fill = "grey90", color = "white", linewidth = 0.1) +
    
    # ✨ 物理填缝机制，杜绝地图白边撕裂
    geom_tile(aes(x = x, y = y, fill = Value), width = 0.12, height = 0.12) +
    
    # 采用全新的蓝-黄-红连续色带
    scale_fill_gradientn(
      colors = cv_palette,
      limits = c(0, cv_max),
      na.value = "transparent",
      name = "Coefficient of Variation (CV)",
      breaks = seq(0, cv_max, length.out = 6),
      labels = scales::number_format(accuracy = 0.01)
    ) +
    
    facet_wrap(~ Scenario, ncol = 2) +
    
    labs(title = title_str) +
    theme_void() +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 15, t = 10)),
      strip.text = element_text(size = 14, face = "bold", margin = ggplot2::margin(b = 8, t = 8)),
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold", vjust = 1),
      legend.text = element_text(size = 10, face = "bold")
    ) +
    coord_sf(ylim = c(-60, 90), expand = FALSE) +
    # 拉长底部图例，使其与你的参考图一致
    guides(fill = guide_colorbar(barwidth = 25, barheight = 0.8, title.position = "top", title.hjust = 0.5, ticks.linewidth = 1))
  
  return(p)
}



# 1. TD 的 SD 图谱
p_td_sd_atlas <- plot_sd_atlas(df_td_sd, "Taxonomic Diversity (TD) Prediction Uncertainty (SD)")
ggsave("Uncertainty_TD_SD_Masked_Atlas.pdf", p_td_sd_atlas, width = 12, height = 9, dpi = 450, bg = "transparent")

# 3. TD 的 CV 图谱
p_td_cv_atlas <- plot_cv_atlas(df_td_cv, "Taxonomic Diversity (TD) Coefficient of Variation (CV)")
ggsave("Uncertainty_TD_CV_Masked_Atlas.png", p_td_cv_atlas, width = 12, height = 9, dpi = 450, bg = "transparent")








 #GFM####

# ==========================================
# 5. GFM (GS) 专属：基于物理坐标系的强力掩膜与可视化
# ==========================================
cat("\n🚀 正在启动 GFM 专属的强力掩膜处理...\n")

# ------------------------------------------
# A. GFM 强力物理掩膜提取函数
# ------------------------------------------
process_gfm_uncertainty_robust <- function(uncert_path, scenario_folder, scenario_name) {
  cat("⏳ 正在进行精准物理掩膜: GFM -", scenario_name, "...\n")
  
  # 1. 读取原始不确定性图层 (SD 或 CV)，降采样，并立即转为 DataFrame (锁死 2 位小数)
  r_uncert <- rast(uncert_path)
  r_uncert_agg <- aggregate(r_uncert, fact = 10, fun = "mean", na.rm = TRUE)
  
  df_uncert <- as.data.frame(r_uncert_agg, xy = TRUE, na.rm = TRUE) %>%
    mutate(x = round(x, 2), y = round(y, 2))
  colnames(df_uncert)[3] <- "Value"
  
  # 2. 读取 GFM 专属的外推预测 strict 图层作为 Mask
  mask_path <- file.path(base_mask_path, "GFM", scenario_folder, 
                         paste0("prediction_ridge_reliable_", strictness, ".tif"))
  if(!file.exists(mask_path)) stop(paste("❌ 找不到 GFM 掩膜文件:", mask_path))
  
  r_mask <- rast(mask_path)
  r_mask_agg <- aggregate(r_mask, fact = 10, fun = "mean", na.rm = TRUE)
  
  # 转为 DataFrame，同样锁死 2 位小数
  df_mask <- as.data.frame(r_mask_agg, xy = TRUE, na.rm = TRUE) %>%
    mutate(x = round(x, 2), y = round(y, 2))
  
  # 3. 💥 核心魔法：抛弃 terra::mask，直接用 inner_join 求物理坐标交集！
  # 只要这个点在掩膜 (df_mask) 里存在，它就被保留；否则直接剔除，实现完美抠图！
  df_final <- inner_join(df_uncert, dplyr::select(df_mask, x, y), by = c("x", "y")) %>%
    mutate(Scenario = scenario_name, Variable = "GFM")
  
  return(df_final)
}

# ------------------------------------------
# B. 定义 GFM 的原始 TIF 路径字典
# ------------------------------------------
# 1. GFM (GS) 的 SD 路径
paths_gfm_sd <- c(
  "2020_SSP245" = "/Volumes/DEEP/tiff/uncertainty/outlier/GS_sd_2020SSP245.tif",
  "2100_SSP126" = "/Volumes/DEEP/tiff/uncertainty/outlier/GS_sd_2100SSP126.tif",
  "2100_SSP245" = "/Volumes/DEEP/tiff/uncertainty/outlier/GS_sd_2100SSP245.tif",
  "2100_SSP370" = "/Volumes/DEEP/tiff/uncertainty/outlier/GS_sd_2100SSP370.tif"
)

# 2. GFM (GS) 的 CV 路径 (根据你刚提供的路径)
paths_gfm_cv <- c(
  "2020_SSP245" = "/Volumes/DEEP/tiff/uncertainty/cv/GS_uncertainty_cv_2020SSP245.tif",
  "2100_SSP126" = "/Volumes/DEEP/tiff/uncertainty/cv/GS_uncertainty_cv_2100SSP126.tif",
  "2100_SSP245" = "/Volumes/DEEP/tiff/uncertainty/cv/GS_uncertainty_cv_2100SSP245.tif",
  "2100_SSP370" = "/Volumes/DEEP/tiff/uncertainty/cv/GS_uncertainty_cv_2100SSP370.tif"
)

# ------------------------------------------
# C. 批量提取并强力掩膜执行
# ------------------------------------------
scen_levels <- c("2020 SSP2-4.5", "2100 SSP1-2.6", "2100 SSP2-4.5", "2100 SSP3-7.0")

cat("\n📥 提取 GFM SD...\n")
df_gfm_sd <- bind_rows(
  process_gfm_uncertainty_robust(paths_gfm_sd["2020_SSP245"], "cropland2020_SSP245", "2020 SSP2-4.5"),
  process_gfm_uncertainty_robust(paths_gfm_sd["2100_SSP126"], "cropland2100_SSP126", "2100 SSP1-2.6"),
  process_gfm_uncertainty_robust(paths_gfm_sd["2100_SSP245"], "cropland2100_SSP245", "2100 SSP2-4.5"),
  process_gfm_uncertainty_robust(paths_gfm_sd["2100_SSP370"], "cropland2100_SSP370", "2100 SSP3-7.0")
) %>% mutate(Scenario = factor(Scenario, levels = scen_levels))

cat("\n📥 提取 GFM CV...\n")
df_gfm_cv <- bind_rows(
  process_gfm_uncertainty_robust(paths_gfm_cv["2020_SSP245"], "cropland2020_SSP245", "2020 SSP2-4.5"),
  process_gfm_uncertainty_robust(paths_gfm_cv["2100_SSP126"], "cropland2100_SSP126", "2100 SSP1-2.6"),
  process_gfm_uncertainty_robust(paths_gfm_cv["2100_SSP245"], "cropland2100_SSP245", "2100 SSP2-4.5"),
  process_gfm_uncertainty_robust(paths_gfm_cv["2100_SSP370"], "cropland2100_SSP370", "2100 SSP3-7.0")
) %>% mutate(Scenario = factor(Scenario, levels = scen_levels))

# ------------------------------------------
# D. 直接复用之前的绘图函数，极速出图
# ------------------------------------------
cat("\n🎨 正在渲染绝对完美掩膜的 GFM 图集...\n")

# 1. GFM 的 SD 图谱 (绿色阶梯图)
p_gfm_sd_atlas <- plot_sd_atlas(df_gfm_sd, "Guild Functional Menhinick (GFM) Prediction Uncertainty (SD)")
ggsave("Uncertainty_GFM_SD_Masked_Atlas.png", p_gfm_sd_atlas, width = 12, height = 9, dpi = 450, bg = "transparent")
print(p_gfm_sd_atlas)

# 2. GFM 的 CV 图谱 (蓝黄红连续图)
p_gfm_cv_atlas <- plot_cv_atlas(df_gfm_cv, "Guild Functional Menhinick (GFM) Coefficient of Variation (CV)")
ggsave("Uncertainty_GFM_CV_Masked_Atlas.png", p_gfm_cv_atlas, width = 12, height = 9, dpi = 450, bg = "transparent")
print(p_gfm_cv_atlas)

cat("🎉 恭喜！GFM (GS) 专属的物理强力掩膜已完美生效！海岸线与边界已精准抠除。\n")



# ==========================================
# 6. 导出 TD 与 GFM 的不确定性 (SD & CV) 汇总统计表格
# ==========================================
cat("\n📊 正在汇总 TD 与 GFM 的掩膜后不确定性统计信息...\n")

# 1. 提取并合并所有数据，同时打上对应的 Metric (SD 或 CV) 标签
df_summary_all <- bind_rows(
  df_td_sd %>% mutate(Metric = "SD (Uncertainty)"),
  df_td_cv %>% mutate(Metric = "CV (Coefficient of Variation)"),
  df_gfm_sd %>% mutate(Metric = "SD (Uncertainty)"),
  df_gfm_cv %>% mutate(Metric = "CV (Coefficient of Variation)")
)

# 2. 分组计算全局统计指标
uncertainty_summary_table <- df_summary_all %>%
  group_by(Variable, Metric, Scenario) %>%
  summarise(
    Min        = round(min(Value, na.rm = TRUE), 4),
    Max        = round(max(Value, na.rm = TRUE), 4),
    Mean       = round(mean(Value, na.rm = TRUE), 4),
    Spatial_SD = round(sd(Value, na.rm = TRUE), 4), # 这里的 SD 指的是不确定性在空间分布上的离散程度
    .groups    = "drop"
  ) %>%
  # 排序：先按变量 (TD/GFM)，再按指标，最后按时间情景排序，保证表格逻辑清晰
  arrange(Variable, desc(Metric), Scenario)

# 3. 打印预览到控制台
cat("\n✨ 统计结果预览:\n")
print(uncertainty_summary_table)

# 4. 导出为精美的 CSV 表格
output_csv_name <- "Supplementary_Table_TD_GFM_Uncertainty_Summary.csv"
write.csv(uncertainty_summary_table, output_csv_name, row.names = FALSE, quote = FALSE)

cat(sprintf("\n✅ TD 与 GFM 的统计表格已成功导出至工作目录:\n   %s\n", file.path(getwd(), output_csv_name)))
