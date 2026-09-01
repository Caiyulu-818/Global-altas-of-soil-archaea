#fd&TD
library(terra)
library(ggplot2)
library(rnaturalearth) # 用于获取高清世界底图
library(patchwork)# 用于最后拼图
library(dplyr)# 用于数据处理

# ==========================================
# 1. 读取和聚合栅格数据
# ==========================================
FD2100370 <- rast("/Volumes/DEEP/tiff/GS2100SS370.rf.tif")
FD2020370 <- rast("/Volumes/DEEP/tiff/GS2020ssp370.rf.tif")
FD2100126 <- rast("/Volumes/DEEP/tiff/GS2100SS126.rf.tif")
FD2020126 <- rast("/Volumes/DEEP/tiff/GS2020SS126.rf.tif")
FD2100245<-rast("/Volumes/DEEP/tiff/GS2100SS245.rf.tif")
FD2020245<-rast("/Volumes/DEEP/tiff/GS2020ssp245.rf.tif")

TD2100370 <- rast("/Volumes/DEEP/tiff/TD2100SS370.rf.tif")
TD2020370 <- rast("/Volumes/DEEP/tiff/TD2020SS370.rf.tif")
TD2100126 <- rast("/Volumes/DEEP/tiff/TD2100SS126.rf.tif")
TD2020126 <- rast("/Volumes/DEEP/tiff/TD2020SS126.rf.tif")
TD2100245<-rast("/Volumes/DEEP/tiff/TD2100SS245.rf.tif")
TD2020245<-rast("/Volumes/DEEP/tiff/TD2020SS245.rf.tif")



# 聚合降采样 (0.1度)
FD2100370_01deg <- aggregate(FD2100370, fact = 10, fun = mean, na.rm = TRUE)
FD2020370_01deg <- aggregate(FD2020370, fact = 10, fun = mean, na.rm = TRUE)
FD2100126_01deg <- aggregate(FD2100126, fact = 10, fun = mean, na.rm = TRUE)
FD2020126_01deg <- aggregate(FD2020126, fact = 10, fun = mean, na.rm = TRUE)
FD2100245_01deg <- aggregate(FD2100245, fact = 10, fun = mean, na.rm = TRUE)
FD2020245_01deg <- aggregate(FD2020245, fact = 10, fun = mean, na.rm = TRUE)

# 聚合降采样 (0.1度)
TD2100370_01deg <- aggregate(TD2100370, fact = 10, fun = mean, na.rm = TRUE)
TD2020370_01deg <- aggregate(TD2020370, fact = 10, fun = mean, na.rm = TRUE)
TD2100126_01deg <- aggregate(TD2100126, fact = 10, fun = mean, na.rm = TRUE)
TD2020126_01deg <- aggregate(TD2020126, fact = 10, fun = mean, na.rm = TRUE)
TD2100245_01deg <- aggregate(TD2100245, fact = 10, fun = mean, na.rm = TRUE)
TD2020245_01deg <- aggregate(TD2020245, fact = 10, fun = mean, na.rm = TRUE)


library(terra)
library(ggplot2)
library(rnaturalearth)
library(patchwork)
library(dplyr)
library(data.table) # fread 需要

# === 前置准备：你的读取和聚合代码保持不变 ===
# (假设你已经跑完了 FD2020245_01deg, TD2020245_01deg 等 aggregate 步骤)
position1.info1 <- fread("/Volumes/DEEP/position/results1231.csv")
landclass<-read.csv("/Volumes/DEEP/position/landclass.csv")
#global2020245####
# 1. 定义数据框提取函数
get_df <- function(r_fd, r_td, scenario_year) {
  df_fd <- as.data.frame(r_fd, xy = TRUE, na.rm = TRUE)
  df_td <- as.data.frame(r_td, xy = TRUE, na.rm = TRUE)
  colnames(df_fd)[3] <- "FD"
  colnames(df_td)[3] <- "TD"
  
  df_merge <- inner_join(df_fd, df_td, by = c("x", "y")) %>%
   mutate(x = round(x, 2), y = round(y, 2), Scenario = scenario_year)
  
  
  df_merge <- df_merge %>%
    mutate(landcover = position1.info1$value) %>%
    mutate(landcover = ifelse(is.na(landcover), 80, landcover)) %>%
    mutate(Habitat = case_when(
      landcover == 10 ~ "Forest",
      landcover == 20 ~ "Shrubland",
      landcover == 30 ~ "Grassland",
      landcover == 40 ~ "Agricultural Land",
      landcover %in% c(50, 60) ~ "Wetland",
      landcover == 70 ~ "Tundra",
      landcover == 80 ~ "Bare Land",
      landcover == 90 ~ "Artificial Surfaces",
      landcover == 100 ~ "Permanent water bodies",
      landcover == 110 ~ "Snow and ice",
      landcover == 254 ~ "Unclassified",
      landcover == 255 ~ "Nodata",
      TRUE~ as.character(landcover)
    ))
  
  return(df_merge)
}

# 2. 提取数据 (特别加入 2020_SSP245 作为 Present)
df_2020_245 <- get_df(FD2020245_01deg, TD2020245_01deg, "Present_2020_SSP245")
df_2100_245 <- get_df(FD2100245_01deg, TD2100245_01deg, "Present_2100_SSP245")
df_2100_126 <- get_df(FD2100126_01deg, TD2100126_01deg, "Present_2100_SSP126")
df_2100_370 <- get_df(FD2100370_01deg, TD2100370_01deg, "Present_2100_SSP370")

# ==========================================
# 1. 数据合并与 5 级分位数计算
# ==========================================
# ✨ 修复 Bug：正确合并两年的数据，确保全局断点一致
all_data_abs <- bind_rows(df_2020_245, df_2100_245,df_2100_126,df_2100_370)

# 3. 计算断点 (25%, 50%, 75%) -> 完美四等分
fd_quantiles <- quantile(all_data_abs$FD, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)
td_quantiles <- quantile(all_data_abs$TD, probs = c(0.25, 0.50, 0.75), na.rm = TRUE)

hist(df_2020_245$FD)
hist(df_2100_245$FD)
hist(df_2100_126$FD)
hist(df_2100_370$FD)
hist(df_2020_245$TD)
hist(df_2100_245$TD)
hist(df_2100_126$TD)
hist(df_2100_370$TD)


# 4. 执行 4 级分类
all_data_abs <- all_data_abs %>%
  mutate(
    FD_class = case_when(
      FD < fd_quantiles[1] ~ 1,
      FD < fd_quantiles[2] ~ 2,
      FD < fd_quantiles[3] ~ 3,
      TRUE ~ 4
    ),
    TD_class = case_when(
      TD < td_quantiles[1] ~ 1,
      TD < td_quantiles[2] ~ 2,
      TD < td_quantiles[3] ~ 3,
      TRUE ~ 4
    ),
    bi_class = paste0(TD_class, "-", FD_class)
  )

# ==========================================
# 终极定制：4x4 完美象限色盘 (16色)
# 逻辑：完美的十字切割，四大象限互不干扰
# ==========================================
biv_colors_abs_16 <- c(
  # 🟩 左上象限：绿色渐变 (功能专精区 -> GFD高，TD低)
  "1-4" = "#4d7e54", # 极值：深森林绿
  "2-4" = "#669877", # 中绿
  "1-3" = "#a4cbb7", # 中绿
  "2-3" = "#cfeadf", # 原点：极浅微光绿 (靠近中心)
  
  # 🟨 右上象限：金黄渐变 (双重繁荣区 -> GFD高，TD高)
  "4-4" = "#E6A800", # 极值：深金黄
  "3-4" = "#FEC44F", # 中黄
  "4-3" = "#FEE391", # 中黄
  "3-3" = "#FFF7BC", # 原点：极浅微光黄 (靠近中心)
  
  # 🟦 左下象限：蓝色渐变 (双重枯竭区 -> GFD低，TD低)
  "1-1" = "#2C7BB6", # 极值：深海蓝
  "2-1" = "#75b5dc", # 中蓝
  "1-2" = "#ABD9E9", # 中蓝
  "2-2" = "#E0F3F8", # 原点：极浅微光蓝 (靠近中心)
  
  # 🚨 右下象限：红色渐变 (掩盖效应区 -> GFD低，TD高) —— 你的核心重灾区！
  "4-1" = "#DE2D26", # 极值：极致血红 (触目惊心！)
  "3-1" = "#d16d5b", # 中红
  "4-2" = "#FC9272", # 中红
  "3-2" = "#FEE5D9" # 原点：极浅微光红 (靠近中心)
)

# ==========================================
# 匹配的高级 4x4 图例
# ==========================================
p_legend_abs <- ggplot(expand.grid(TD_class = 1:4, FD_class = 1:4) %>% mutate(bi_class = paste0(TD_class, "-", FD_class)),
                       aes(x = TD_class, y = FD_class, fill = bi_class)) +
  geom_tile(color = "white", size = 1) + # 粗白线，强化十字切割感
  scale_fill_manual(values = biv_colors_abs_16) +
  
  # 修改刻度为 4 级
  scale_x_continuous(breaks = c(1, 2.5, 4), labels = c("Low", "Med", "High")) +
  scale_y_continuous(breaks = c(1, 2.5, 4), labels = c("Low", "Med", "High")) +
  
  labs(x = "Taxonomic Diversity (TD) →", y = "Guild Functional Menhinick (GFM) →") +
  theme_void() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 12, face = "bold", color = "black"),
    axis.text = element_text(size = 10, face = "bold", color = "grey30"),
    
    # 拉开距离防重叠
    axis.title.y = element_text(angle = 90, hjust = 0.5, margin = ggplot2::margin(r = 15)),
    axis.title.x = element_text(hjust = 0.5, margin = ggplot2::margin(t = 15)),
    
    axis.text.y = element_text(hjust = 1, margin = ggplot2::margin(r = 5)),
    axis.text.x = element_text(vjust = 1, margin = ggplot2::margin(t = 5))
  ) +
  coord_fixed()

p_legend_abs


# ==========================================
#====================
# 5. 核心制图与拼接函数
# ==========================================
world_map <- ne_countries(scale = "medium", returnclass = "sf") %>% filter(admin != "Antarctica")

plot_bivariate_map <- function(plot_data, title_text) {
  p_map <- ggplot() +
    geom_sf(data = world_map, fill = "grey90", color = "white", size = 0.1) +
    geom_tile(data = plot_data, aes(x = x, y = y, fill = bi_class)) +
    scale_fill_manual(values = biv_colors_abs_16) +
    labs(title = title_text) + theme_void() +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) +
    coord_sf(ylim = c(-60, 90), expand = FALSE)
  
  # ✨ 修复 Bug：此处必须使用 p_legend_abs，不能写 p_delta_legend
  final_plot <- p_map + patchwork::inset_element(p_legend_abs, left = 0.02, bottom = 0.05, right = 0.25, top = 0.35, align_to = 'plot')
  return(final_plot)
}
write_rds(all_data_abs,"./arc/01.arcinfo/F3/FD/FDdriver/all.data.abs.rds")
# ==========================================
# 6. 生成并保存 2020 和 2100 的全局绝对值地图
# ==========================================
p_present_2020 <- plot_bivariate_map(filter(all_data_abs, Scenario == "Present_2020_SSP245"), "Global Present Baseline (2020)")
p_present_2100 <- plot_bivariate_map(filter(all_data_abs, Scenario == "Present_2100_SSP245"), "2100 SSP2-4.5")
p_present_2020
p_present_2100
p_2100_126 <- plot_bivariate_map(filter(all_data_abs, Scenario == "Present_2100_SSP126"), "2100 SSP1-2.6")
p_2100_370 <- plot_bivariate_map(filter(all_data_abs, Scenario == "Present_2100_SSP370"), "2100 SSP3-7.0")
p_2100_126
p_2100_370

# 拼合上下两张图一起输出，方便对比
combined_abs_plot1 <- (p_present_2020)/ (p_2100_126)/(p_2100_370)
combined_abs_plot1
ggsave("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/ok/newssp1262100Absolute_5x5.pdf", p_2100_126, width = 14, height = 8, dpi = 300,bg="transparent")
ggsave("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/ok/newssp3702100Absolute_5x5.pdf", p_2100_370, width = 14, height = 8, dpi = 300,bg="transparent")

ggsave("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/ok/newPresent_2020Absolute_5x5.pdf", p_present_2020, width = 14, height = 8, dpi = 300,bg="transparent")
ggsave("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/ok/newPresent_2100_Absolute_5x5.pdf",p_present_2100, width = 14, height = 8, dpi = 300,bg="transparent")
ggsave("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/ok/Present_2020_and_2100_Absolute_5x5.pdf", combined_abs_plot1, width = 14, height = 16, dpi = 300,bg="transparent")


#hotspot####Van Nuland et al. 2025nature
library(dplyr)

# ==========================================
# 1. 严格提取 2020 Present 基准数据
# ==========================================
# 顶刊准则：阈值必须由基准年决定，作为衡量未来的统一标尺
df_baseline <- all_data_abs %>% 
  filter(Scenario == "Present_2020_SSP245")

# ==========================================
# 2. 计算上 95% 分位数 (Top 5%) 作为 Cut-off
# ==========================================
cutoff_gfd <- quantile(df_baseline$FD, probs = 0.95, na.rm = TRUE)
cutoff_td  <- quantile(df_baseline$TD, probs = 0.95, na.rm = TRUE)


# 在控制台打印出这两个决定命运的数值！
cat("=== Global Hotspot Cut-off Values (95th Percentile) ===\n")
cat("GFM Cut-off: ", round(cutoff_gfd, 4), "\n")#0.6499 
cat("TD Cut-off:  ", round(cutoff_td, 4), "\n")# 2.6301
# ==========================================
# 3. 鉴定每个情景下的 Hotspots 
# ==========================================
# 我们在原数据框里新增两列标签，标记该网格在当前情景下是否属于热点
all_data_abs <- all_data_abs %>%
  mutate(
    is_GFD_hotspot = ifelse(FD >= cutoff_gfd, "Hotspot", "Non-hotspot"),
    is_TD_hotspot  = ifelse(TD >= cutoff_td,  "Hotspot", "Non-hotspot"),
    
    # 你还可以定义“双重热点区 (Coupled Hotspots)” - 即既是物种热点又是功能热点！
    is_Coupled_hotspot = ifelse(FD >= cutoff_gfd & TD >= cutoff_td, "Coupled Hotspot", "Non-hotspot")
  )

# 看看在不同情景下，GFD 热点网格的数量怎么变化？
table(all_data_abs$Scenario, all_data_abs$is_TD_hotspot)

write_rds(all_data_abs,"hotspot.all.abs.rds")

# 安装 rgeoda
#install.packages("rgeoda")
library(sf)
library(rgeoda)
library(dplyr)

# 1. 假设你已经把全球栅格转成了 sf 点对象 (仅保留陆地有效值，剔除 NA)
# 这里用你之前代码生成的 df_baseline，将其转为 sf
sf_baseline <- st_as_sf(df_baseline, coords = c("x", "y"), crs = 4326)

# 2. 建立空间权重矩阵 (Spatial Weights) 
# 因为是全球数据，使用 K-Nearest Neighbors (KNN，比如 K=8 找周围一圈) 是最稳的，防止孤岛报错
knn_weights <- knn_weights(sf_baseline, k = 8)

# 3. 计算 Local Getis-Ord Gi* # rgeoda 的 local_g 函数极其强大，速度极快
g_star <- local_g(knn_weights, sf_baseline["FD"])

# 4. 提取 P 值并手动进行 FDR 校正！
# 顶刊核心操作就在这一行：
fdr_p_values <- p.adjust(g_star$p_vals, method = "fdr")

# 5. 鉴定最终的统计学热点 (95% 置信度)
# 条件：Z 分数为正 (高值聚集)，且 FDR 校正后的 P 值 < 0.05
# ==========================================
# 5. 鉴定最终的统计学热点并进行严格的置信度分级 (ArcGIS 标准)
# ==========================================
sf_baseline <- sf_baseline %>%
  mutate(
    Z_score = g_star$lisa_vals,
    P_value_FDR = fdr_p_values,
    
    # 🚨 顶刊核心：利用 Z 值的正负和 FDR 校正后的 P 值进行 7 级阶梯分级
    # 注意 case_when 的执行顺序，必须先过滤最严格的 0.01，再过滤 0.05
    Hotspot_Status = case_when(
      # --- 高值聚集区 (Hotspots) ---
      Z_score > 0 & P_value_FDR < 0.01 ~ "Hotspot (99% Confidence)",
      Z_score > 0 & P_value_FDR < 0.05 ~ "Hotspot (95% Confidence)",
      Z_score > 0 & P_value_FDR < 0.10 ~ "Hotspot (90% Confidence)",
      
      # --- 低值聚集区 (Coldspots) ---
      Z_score < 0 & P_value_FDR < 0.01 ~ "Coldspot (99% Confidence)",
      Z_score < 0 & P_value_FDR < 0.05 ~ "Coldspot (95% Confidence)",
      Z_score < 0 & P_value_FDR < 0.10 ~ "Coldspot (90% Confidence)",
      
      # --- 统计学不显著区 (随机分布的噪音) ---
      TRUE                             ~ "Not Significant"
    ),
    
    # 为了方便后续画图和统计，把状态变成有序因子 (Ordered Factor)
    Hotspot_Status = factor(Hotspot_Status, levels = c(
      "Hotspot (99% Confidence)", "Hotspot (95% Confidence)", "Hotspot (90% Confidence)",
      "Not Significant",
      "Coldspot (90% Confidence)", "Coldspot (95% Confidence)", "Coldspot (99% Confidence)"
    ))
  )
write_rds(sf_baseline,"hotspot.fdr.level.rds")
# 查看分级后的网格数量统计！
table(sf_baseline$Hotspot_Status)

#deltaglobal####

# ==========================================
# 1. 以 2020_SSP245 为统一基准计算差值 (Delta)
# ==========================================
# ✨ 统一减去 Present 的数据！
Delta_GFD_126_vs_Present <- FD2100126_01deg - FD2020245_01deg
Delta_TD_126_vs_Present <- TD2100126_01deg - TD2020245_01deg

Delta_GFD_370_vs_Present <- FD2100370_01deg - FD2020245_01deg
Delta_TD_370_vs_Present <- TD2100370_01deg - TD2020245_01deg

Delta_GFD_245_vs_Present <- FD2100245_01deg - FD2020245_01deg
Delta_TD_245_vs_Present <- TD2100245_01deg - TD2020245_01deg

# 定义提取差值数据的函数
get_delta_df <- function(r_delta_fd, r_delta_td, scenario_name) {
  df_fd <- as.data.frame(r_delta_fd, xy = TRUE, na.rm = TRUE)
  df_td <- as.data.frame(r_delta_td, xy = TRUE, na.rm = TRUE)
  colnames(df_fd)[3] <- "Delta_GFD"
  colnames(df_td)[3] <- "Delta_TD"
  
  df_merge <- inner_join(df_fd, df_td, by = c("x", "y")) %>%
    mutate(x = round(x, 2), y = round(y, 2), Scenario = scenario_name)
  
  df_merge <- df_merge %>%
    mutate(landcover = position1.info1$value) %>%
    mutate(landcover = ifelse(is.na(landcover), 80, landcover)) %>%
    mutate(Habitat = case_when(
      landcover == 10 ~ "Forest",
      landcover == 20 ~ "Shrubland",
      landcover == 30 ~ "Grassland",
      landcover == 40 ~ "Agricultural Land",
      landcover %in% c(50, 60) ~ "Wetland",
      landcover == 70 ~ "Tundra",
      landcover == 80 ~ "Bare Land",
      landcover == 90 ~ "Artificial Surfaces",
      landcover == 100 ~ "Permanent water bodies",
      landcover == 110 ~ "Snow and ice",
      landcover == 254 ~ "Unclassified",
      landcover == 255 ~ "Nodata",
      TRUE~ as.character(landcover)
    ))
  
  return(df_merge)
}


df_delta_126 <- get_delta_df(Delta_GFD_126_vs_Present, Delta_TD_126_vs_Present, "Delta_SSP126_vs_Present")
df_delta_370 <- get_delta_df(Delta_GFD_370_vs_Present, Delta_TD_370_vs_Present, "Delta_SSP370_vs_Present")
df_delta_245 <- get_delta_df(Delta_GFD_245_vs_Present, Delta_TD_245_vs_Present, "Delta_SSP245_vs_Present")

all_delta_data <- bind_rows(df_delta_126, df_delta_370,df_delta_245)
write_rds(all_delta_data,"./arc/01.arcinfo/F3/FD/FDdriver/ok/all_delta_data.rds")
all_delta_data<-readRDS("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/ok/all_delta_data.rds")
# ==========================================
# 2. 基于绝对零点的 3x3 分类 (中间代表 0)
# ==========================================
#th_buffer <- 0.04 # 缓冲阈值：-0.02 到 +0.02 之间算作“没有变化(0)”

all_delta_data <- all_delta_data %>%
  mutate(
    FD_class = case_when(
      Delta_GFD < -0.02~ 1, # 1: Loss
      Delta_GFD >0.02  ~ 3, # 2: Zero / Stable
      TRUE ~ 2# 3: Gain
    ),
    TD_class = case_when(
      Delta_TD < -0.02~ 1,
      Delta_TD >0.02 ~ 3,
      TRUE ~ 2
    ),
    bi_class = paste0(TD_class, "-", FD_class)
  )

# ✨ 第一步：构建带有“区间文本”的绘图数据框
legend_data <- expand.grid(TD_class = 1:3, FD_class = 1:3) %>%
  mutate(
    bi_class = paste0(TD_class, "-", FD_class),
    
    # 提取 X 轴 (TD) 的具体范围
    td_range = case_when(
      TD_class == 1 ~ "Δ TD < -0.02",
      TD_class == 2 ~ "-0.02 < Δ TD < 0.02 ",
      TD_class == 3 ~ "Δ TD > 0.02"
    ),
    
    # 提取 Y 轴 (GFD) 的具体范围
    fd_range = case_when(
      FD_class == 1 ~ "Δ GFM < -0.02",
      FD_class == 2 ~ "-0.02 < Δ GFM < 0.02",
      FD_class == 3 ~ "Δ GFM > 0.02"
    ),
    
    # 将 X 和 Y 的范围拼接成两行文字 (\n 换行)
    label_text = paste0(td_range, "\n", fd_range),
    
    # 🎨 智能配色：1-1(深绿) 和 3-1(深红) 色块太深，用白色字体，其他用黑色
    text_color = ifelse(bi_class %in% c("1-1", "3-1"), "black", "black")
  )
# ==========================================
# 3. 差值专用的 3x3 高对比度图例盘 (重点！)
# ==========================================
# 3. 差值专用的 3x3 原始经典图例盘 (Green-Red-Yellow)
# ==========================================
# 替换回最原始的 9 宫格颜色
library(dplyr)
library(ggplot2)



# 1. 定义经典的 9 宫格颜色
biv_colors_delta <- c(
  "1-3" = "#44DE44", "2-3" = "#9BE241", "3-3" = "#F2E62E", # 顶层
  "1-2" = "#5C9432", "2-2" = "#A4A836", "3-2" = "#E19C29", # 中层
  "1-1" = "#435222", "2-1" = "#886B26", "3-1" = "#CE3327" # 底层
)
biv_colors_delta <- c(
  "1-3" = "#cae0a2", "2-3" = "#c2c8d8", "3-3" = "#ffd19d", # 顶层
  "1-2" = "#8CAD8D", "2-2" = "#B1B3D9", "3-2" = "#f39b7e", # 中层
  "1-1" = "#455D1C", "2-1" = "#43324a", "3-1" = "#dc0000" # 底层
)
#c2c8d8"#dc0000
#9b9ac1
#cabad7
# 2. 赋予每个格子简洁的“生态学意义” (而不是冷冰冰的数学公式)
legend_data <- expand.grid(TD_class = 1:3, FD_class = 1:3) %>%
  mutate(
    bi_class = paste0(TD_class, "-", FD_class),
    
    # 用简短的定性描述代替数学公式
    meaning = case_when(
      bi_class == "3-3" ~ "Coupled\nGain",# 协同增长 (双赢)
      bi_class == "1-1" ~ "Coupled\nLoss",# 协同下降 (双输)
      bi_class == "3-1" ~ "Δ Decoupling",# 严重解耦 (分类涨，功能跌，最危险)
      bi_class == "1-3" ~ "Functional\nBuffer",# 功能缓冲 (分类跌，功能涨)
      TRUE ~ ""# 中间过渡地带留白，保持清爽！
    ),
    
    # 智能反色：在深色块上用白色字体
    text_color = ifelse(bi_class %in% c("1-1", "3-1"), "white", "black")
  )

# 3. 绘制高级降噪版图例
p_delta_legend <- ggplot(legend_data, aes(x = TD_class, y = FD_class, fill = bi_class)) +
  geom_tile(color = "white", linewidth = 0.8) +
  
  # 贴上极简的生态学标签
  geom_text(aes(label = meaning, color = text_color),
            size = 6, fontface = "bold", lineheight = 0.9) +
  scale_color_identity() +
  scale_fill_manual(values = biv_colors_delta) +
  
  # 🚨 核心更新：X 轴和 Y 轴分别展示其专属的独立阈值
  scale_x_continuous(breaks = 1:3, 
                     labels = c("Loss\n(Δ < -0.02)", "Stable\n(± 0.02)", "Gain\n(Δ > 0.02)")) +
  scale_y_continuous(breaks = 1:3, 
                     labels = c("Loss\n(Δ < -0.02)", "Stable\n(± 0.02)", "Gain\n(Δ > 0.02)")) +
  
  labs(x = expression(bold(Delta~"Taxonomic Diversity (TD)")),
       y = expression(bold(Delta~"Guild Functional Menhinick (GFM)"))) +
  
  theme_void() +
  theme(
    legend.position = "none",
    
    # 轴标题优化
    axis.title.x = element_text(hjust = 0.5, margin = ggplot2::margin(t = 10), size = 18),
    axis.title.y = element_text(angle = 90, hjust = 0.5, margin = ggplot2::margin(r = 10), size = 12),
    
    # 轴刻度文本优化 (包含具体阈值)
    axis.text.x = element_text(vjust = 1, hjust = 0.5, margin = ggplot2::margin(t = 5), size = 18, color = "black", lineheight = 0.8),
    axis.text.y = element_text(hjust = 1, vjust = 0.5, margin = ggplot2::margin(r = 5), size = 18, color = "black", lineheight = 0.8)
  ) +
  coord_fixed()

print(p_delta_legend)

ggsave("p_delta_legend.pdf",height=8,width = 8,bg="transparent")
# ==========================================
# 4. 绘制差值地图并拼接
# ==========================================
plot_delta_map <- function(plot_data, title_text) {
  p_map <- ggplot() +
    geom_sf(data = world_map, fill = "grey90", color = "white", size = 0.1) +
    geom_tile(data = plot_data, aes(x = x, y = y, fill = bi_class
                                    )) +
    scale_fill_manual(values = biv_colors_delta) +
    labs(title = title_text) + theme_void() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) +
    coord_sf(ylim = c(-60, 90), expand = FALSE)
  
  #return(p_map + patchwork::inset_element(p_delta_legend, left = 0.02, bottom = 0.04, right = 0.36, top = 0.36, align_to ="plot"))
}

p_delta_126_final <- plot_delta_map(filter(all_delta_data, Scenario == "Delta_SSP126_vs_Present"), "Shift to 2100 SSP1-2.6 (baseline: 2020 SSP2-4.5)")
p_delta_370_final <- plot_delta_map(filter(all_delta_data, Scenario == "Delta_SSP370_vs_Present"), "Shift to 2100 SSP3-7.0 (baseline: 2020 SSP2-4.5)")
p_delta_245_final <- plot_delta_map(filter(all_delta_data, Scenario == "Delta_SSP245_vs_Present"), "Shift to 2100 SSP2-4.5 (baseline: 2020 SSP2-4.5)")

p_delta_126_final
p_delta_370_final
p_delta_245_final

library(dplyr)
# 1. Agricultural Land
data_crop_delta <- filter(all_delta_data, Habitat == "Agricultural Land")
p_crop_delta_126 <- plot_delta_map(filter(data_crop_delta, Scenario == "Delta_SSP126_vs_Present"), "Agricultural Land Shift (SSP1-2.6)")
p_crop_delta_370 <- plot_delta_map(filter(data_crop_delta, Scenario == "Delta_SSP370_vs_Present"), "Agricultural Land Shift (SSP3-7.0)")
p_crop_delta_245 <- plot_delta_map(filter(data_crop_delta, Scenario == "Delta_SSP245_vs_Present"), "Agricultural Land Shift (SSP2-4.5)")

# 2. Bare Land
data_bare_delta <- filter(all_delta_data, Habitat == "Bare Land")
p_bare_delta_126 <- plot_delta_map(filter(data_bare_delta, Scenario == "Delta_SSP126_vs_Present"), "Bare Land Shift (SSP1-2.6)")
p_bare_delta_370 <- plot_delta_map(filter(data_bare_delta, Scenario == "Delta_SSP370_vs_Present"), "Bare Land Shift (SSP3-7.0)")
p_bare_delta_245 <- plot_delta_map(filter(data_bare_delta, Scenario == "Delta_SSP245_vs_Present"), "Bare Land Shift (SSP2-4.5)")

# 3. Snow and ice
data_snow_delta <- filter(all_delta_data, Habitat == "Snow and ice")
p_snow_delta_126 <- plot_delta_map(filter(data_snow_delta, Scenario == "Delta_SSP126_vs_Present"), "Snow and ice Shift (SSP1-2.6)")
p_snow_delta_370 <- plot_delta_map(filter(data_snow_delta, Scenario == "Delta_SSP370_vs_Present"), "Snow and ice Shift (SSP3-7.0)")
p_snow_delta_245 <- plot_delta_map(filter(data_snow_delta, Scenario == "Delta_SSP245_vs_Present"), "Snow and ice Shift (SSP2-4.5)")

# 4. Permanent water bodies
data_water_delta <- filter(all_delta_data, Habitat == "Permanent water bodies")
p_water_delta_126 <- plot_delta_map(filter(data_water_delta, Scenario == "Delta_SSP126_vs_Present"), "Permanent water bodies Shift (SSP1-2.6)")
p_water_delta_370 <- plot_delta_map(filter(data_water_delta, Scenario == "Delta_SSP370_vs_Present"), "Permanent water bodies Shift (SSP3-7.0)")
p_water_delta_245 <- plot_delta_map(filter(data_water_delta, Scenario == "Delta_SSP245_vs_Present"), "Permanent water bodies Shift (SSP2-4.5)")

# 5. Tundra
data_tundra_delta <- filter(all_delta_data, Habitat == "Tundra")
p_tundra_delta_126 <- plot_delta_map(filter(data_tundra_delta, Scenario == "Delta_SSP126_vs_Present"), "Tundra Shift (SSP1-2.6)")
p_tundra_delta_370 <- plot_delta_map(filter(data_tundra_delta, Scenario == "Delta_SSP370_vs_Present"), "Tundra Shift (SSP3-7.0)")
p_tundra_delta_245 <- plot_delta_map(filter(data_tundra_delta, Scenario == "Delta_SSP245_vs_Present"), "Tundra Shift (SSP2-4.5)")

# 6. Grassland
data_grass_delta <- filter(all_delta_data, Habitat == "Grassland")
p_grass_delta_126 <- plot_delta_map(filter(data_grass_delta, Scenario == "Delta_SSP126_vs_Present"), "Grassland Shift (SSP1-2.6)")
p_grass_delta_370 <- plot_delta_map(filter(data_grass_delta, Scenario == "Delta_SSP370_vs_Present"), "Grassland Shift (SSP3-7.0)")
p_grass_delta_245 <- plot_delta_map(filter(data_grass_delta, Scenario == "Delta_SSP245_vs_Present"), "Grassland Shift (SSP2-4.5)")

# 7. Wetland
data_wetland_delta <- filter(all_delta_data, Habitat == "Wetland")
p_wetland_delta_126 <- plot_delta_map(filter(data_wetland_delta, Scenario == "Delta_SSP126_vs_Present"), "Wetland Shift (SSP1-2.6)")
p_wetland_delta_370 <- plot_delta_map(filter(data_wetland_delta, Scenario == "Delta_SSP370_vs_Present"), "Wetland Shift (SSP3-7.0)")
p_wetland_delta_245 <- plot_delta_map(filter(data_wetland_delta, Scenario == "Delta_SSP245_vs_Present"), "Wetland Shift (SSP2-4.5)")

# 8. Forest
data_forest_delta <- filter(all_delta_data, Habitat == "Forest")
p_forest_delta_126 <- plot_delta_map(filter(data_forest_delta, Scenario == "Delta_SSP126_vs_Present"), "Forest Shift (SSP1-2.6)")
p_forest_delta_370 <- plot_delta_map(filter(data_forest_delta, Scenario == "Delta_SSP370_vs_Present"), "Forest Shift (SSP3-7.0)")
p_forest_delta_245 <- plot_delta_map(filter(data_forest_delta, Scenario == "Delta_SSP245_vs_Present"), "Forest Shift (SSP2-4.5)")

# 9. Artificial Surfaces
data_art_delta <- filter(all_delta_data, Habitat == "Artificial Surfaces")
p_art_delta_126 <- plot_delta_map(filter(data_art_delta, Scenario == "Delta_SSP126_vs_Present"), "Artificial Surfaces Shift (SSP1-2.6)")
p_art_delta_370 <- plot_delta_map(filter(data_art_delta, Scenario == "Delta_SSP370_vs_Present"), "Artificial Surfaces Shift (SSP3-7.0)")
p_art_delta_245 <- plot_delta_map(filter(data_art_delta, Scenario == "Delta_SSP245_vs_Present"), "Artificial Surfaces Shift (SSP2-4.5)")

# 10. Shrubland
data_shrub_delta <- filter(all_delta_data, Habitat == "Shrubland")
p_shrub_delta_126 <- plot_delta_map(filter(data_shrub_delta, Scenario == "Delta_SSP126_vs_Present"), "Shrubland Shift (SSP1-2.6)")
p_shrub_delta_370 <- plot_delta_map(filter(data_shrub_delta, Scenario == "Delta_SSP370_vs_Present"), "Shrubland Shift (SSP3-7.0)")
p_shrub_delta_245 <- plot_delta_map(filter(data_shrub_delta, Scenario == "Delta_SSP245_vs_Present"), "Shrubland Shift (SSP2-4.5)")

# 11paddy
data_paddy_delta <- filter(paddy_delta_data , Habitat == "Paddy soil")
p_paddy_delta_126 <- plot_delta_map(filter(paddy_delta_data, Scenario == "Delta_SSP126_vs_Present"), "Paddy soil Shift (SSP1-2.6)")
p_paddy_delta_370 <- plot_delta_map(filter(paddy_delta_data, Scenario == "Delta_SSP370_vs_Present"), "Paddy soil Shift (SSP3-7.0)")
p_paddy_delta_245 <- plot_delta_map(filter(paddy_delta_data, Scenario == "Delta_SSP245_vs_Present"), "Paddy soil Shift (SSP2-4.5)")


# （假设你还有整体的图表）
ggsave("Delta_ssp126_from_Present2020_3x34.pdf", p_delta_126_final, width = 16, height = 8, dpi = 300,bg="transparent")
ggsave("Delta_ssp370_from_Present2020_3x34.pdf", p_delta_370_final, width = 16, height = 8, dpi = 300,bg="transparent")
# 如果你有 245 的总图，可以把下一行取消注释：
ggsave("Delta_ssp245_from_Present2020_3x34.pdf", p_delta_245_final, width = 16, height = 8, dpi = 300,bg="transparent")

# 1. Agricultural Land
ggsave("Delta_Agricultural_Land_SSP126.pdf", p_crop_delta_126, width = 16, height = 8, dpi = 300,bg="transparent")
ggsave("Delta_Agricultural_Land_SSP245.pdf", p_crop_delta_245, width = 16, height = 8, dpi = 300,bg="transparent")
ggsave("Delta_Agricultural_Land_SSP370.pdf", p_crop_delta_370, width = 16, height = 8, dpi = 300,bg="transparent")

# 2. Bare Land
ggsave("Delta_Bare_Land_SSP126.pdf", p_bare_delta_126, width = 16, height = 8, dpi = 300,bg="transparent")
ggsave("Delta_Bare_Land_SSP245.pdf", p_bare_delta_245, width = 16, height = 8, dpi = 300,bg="transparent")
ggsave("Delta_Bare_Land_SSP370.pdf", p_bare_delta_370, width = 16, height = 8, dpi = 300,bg="transparent")

# 3. Snow and ice
ggsave("Delta_Snow_and_ice_SSP126.pdf", p_snow_delta_126, width = 16, height = 8, dpi = 300,bg="transparent")
ggsave("Delta_Snow_and_ice_SSP245.pdf", p_snow_delta_245, width = 16, height = 8, dpi = 300,bg="transparent")
ggsave("Delta_Snow_and_ice_SSP370.pdf", p_snow_delta_370, width = 16, height = 8, dpi = 300,bg="transparent")

# 4. Permanent water bodies
ggsave("Delta_Water_Bodies_SSP126.pdf", p_water_delta_126, width = 16, height = 8, dpi = 300,bg="transparent")
ggsave("Delta_Water_Bodies_SSP245.pdf", p_water_delta_245, width = 16, height = 8, dpi = 300,bg="transparent")
ggsave("Delta_Water_Bodies_SSP370.pdf", p_water_delta_370, width = 16, height = 8, dpi = 300,bg="transparent")

# 5. Tundra
ggsave("Delta_Tundra_SSP126.pdf", p_tundra_delta_126, width = 16, height = 8, dpi = 300,bg="transparent")
ggsave("Delta_Tundra_SSP245.pdf", p_tundra_delta_245, width = 16, height = 8, dpi = 300,bg="transparent")
ggsave("Delta_Tundra_SSP370.pdf", p_tundra_delta_370, width = 16, height = 8, dpi = 300,bg="transparent")

# 6. Grassland
ggsave("Delta_Grassland_SSP126.pdf", p_grass_delta_126, width = 16, height = 8, dpi = 300,bg="transparent")
ggsave("Delta_Grassland_SSP245.pdf", p_grass_delta_245, width = 16, height = 8, dpi = 300,bg="transparent")
ggsave("Delta_Grassland_SSP370.pdf", p_grass_delta_370, width = 16, height = 8, dpi = 300,bg="transparent")

# 7. Wetland
ggsave("Delta_Wetland_SSP126.pdf", p_wetland_delta_126, width = 16, height = 8, dpi = 300,bg="transparent")
ggsave("Delta_Wetland_SSP245.pdf", p_wetland_delta_245, width = 16, height = 8, dpi = 300,bg="transparent")
ggsave("Delta_Wetland_SSP370.pdf", p_wetland_delta_370, width = 16, height = 8, dpi = 300,bg="transparent")

# 8. Forest
ggsave("Delta_Forest_SSP126.pdf", p_forest_delta_126, width = 16, height = 8, dpi = 300,bg="transparent")
ggsave("Delta_Forest_SSP245.pdf", p_forest_delta_245, width = 16, height = 8, dpi = 300,bg="transparent")
ggsave("Delta_Forest_SSP370.pdf", p_forest_delta_370, width = 16, height = 8, dpi = 300,bg="transparent")

# 9. Artificial Surfaces
ggsave("Delta_Artificial_Surfaces_SSP126.pdf", p_art_delta_126, width = 16, height = 10, dpi = 300,bg="transparent")
ggsave("Delta_Artificial_Surfaces_SSP245.pdf", p_art_delta_245, width = 16, height = 8, dpi = 300,bg="transparent")
ggsave("Delta_Artificial_Surfaces_SSP370.pdf", p_art_delta_370, width = 16, height = 8, dpi = 300,bg="transparent")

# 10. Shrubland
ggsave("Delta_Shrubland_SSP126.pdf", p_shrub_delta_126, width = 16, height = 8, dpi = 300,bg="transparent")
ggsave("Delta_Shrubland_SSP245.pdf", p_shrub_delta_245, width = 16, height = 8, dpi = 300,bg="transparent")
ggsave("Delta_Shrubland_SSP370.pdf", p_shrub_delta_370, width = 16, height = 8, dpi = 300,bg="transparent")

# 11 paddy
ggsave("Delta_paddy_SSP126.pdf", p_paddy_delta_126, width = 16, height = 8, dpi = 300,bg="transparent")
ggsave("Delta_paddy_SSP245.pdf", p_paddy_delta_245, width = 16, height = 8, dpi = 300,bg="transparent")
ggsave("Delta_paddy_SSP370.pdf", p_paddy_delta_370, width = 16, height = 8, dpi = 300,bg="transparent")

#GFD /TD global altas####
library(terra)
library(ggplot2)
library(rnaturalearth)
library(patchwork)
library(dplyr)

# 1. 计算差值 (Delta)
# 注意：用 2100 减去 2020 Present
Delta_TD_126 <- TD2100126_01deg - TD2020245_01deg
Delta_TD_370 <- TD2100370_01deg - TD2020245_01deg
Delta_TD_245 <- TD2100245_01deg - TD2020245_01deg

Delta_GFD_126 <- FD2100126_01deg - FD2020245_01deg
Delta_GFD_370 <- FD2100370_01deg - FD2020245_01deg
Delta_GFD_245 <- FD2100245_01deg - FD2020245_01deg


# 2. 统一的数据框提取函数 (方便绘图)
rast_to_df <- function(r, var_name, scenario_name) {
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  colnames(df)[3] <- "Value"
  df$Variable <- var_name
  df$Scenario <- scenario_name
  return(df)
}

# --- 提取差值数据 ---
df_delta_td_126 <- rast_to_df(Delta_TD_126, "Delta_TD", "2100 SSP126 - 2020")
df_delta_td_370 <- rast_to_df(Delta_TD_370, "Delta_TD", "2100 SSP370 - 2020")
df_delta_td_245 <- rast_to_df(Delta_TD_245, "Delta_TD", "2100 SSP245 - 2020")

df_delta_gfd_126 <- rast_to_df(Delta_GFD_126, "Delta_GFD", "2100 SSP126 - 2020")
df_delta_gfd_370 <- rast_to_df(Delta_GFD_370, "Delta_GFD", "2100 SSP370 - 2020")
df_delta_gfd_245 <- rast_to_df(Delta_GFD_245, "Delta_GFD", "2100 SSP245 - 2020")

df_delta_GFD <- bind_rows( df_delta_gfd_245,df_delta_gfd_126, df_delta_gfd_370)
df_delta_TD <- bind_rows( df_delta_td_245,df_delta_td_126, df_delta_td_370 )

# --- 提取绝对值数据 ---
df_abs_td_2020 <- rast_to_df(TD2020245_01deg, "Absolute_TD", "Present")
df_abs_td_126 <- rast_to_df(TD2100126_01deg, "Absolute_TD", "2100 SSP126")
df_abs_td_370 <- rast_to_df(TD2100370_01deg, "Absolute_TD", "2100 SSP370")
df_abs_td_245 <- rast_to_df(TD2100245_01deg, "Absolute_TD", "2100 SSP245")


df_abs_gfd_2020 <- rast_to_df(FD2020245_01deg, "Absolute_GFD", "Present")
df_abs_gfd_126 <- rast_to_df(FD2100126_01deg, "Absolute_GFD", "2100 SSP126")
df_abs_gfd_370 <- rast_to_df(FD2100370_01deg, "Absolute_GFD", "2100 SSP370")
df_abs_gfd_245 <- rast_to_df(FD2100245_01deg, "Absolute_GFD", "2100 SSP245")

df_abs_td <- bind_rows(df_abs_td_2020,df_abs_td_245, df_abs_td_126, df_abs_td_370)
df_abs_gfd<-bind_rows(df_abs_gfd_2020, df_abs_gfd_245,df_abs_gfd_126, df_abs_gfd_370)

world_map <- ne_countries(scale = "medium", returnclass = "sf") %>% filter(admin != "Antarctica")

# 1. 定义 TD 的差值色盘 (蓝 - 白 - 红)
# 定义 TD 的差值断点与色盘
# 白色区间将完美锁定在 [-0.005, 0.005]
breaks_td <- c(-0.1, -0.02, -0.005, 0.005, 0.05, 0.2)
colors_td <- c(
  "#226b68", # 极值：深青绿 (< -0.03)
  "#059f86", # 中值：中青绿 (-0.03 ~ -0.01)
  "#91d1c2", # 浅值：浅青绿 (-0.01 ~ -0.002)
  "#f7f7f7", # 零值缓冲带：纯白 (-0.002 ~ 0.002) - 完美包含中位数！
  "#dfc27d", # 浅值：浅黄棕 (0.002 ~ 0.01)
  "#bf812d", # 中值：中等泥棕 (0.01 ~ 0.03)
  "#8c510a" # 极值：极深焦棕 (> 0.03)
)

# 2. 定义 GFD 的差值色盘 (青绿 - 白 - 焦棕)
# 完美贴合数据的 6 个内部断点
breaks_gfd <- c(-0.10, -0.02, -0.005, 0.005, 0.02, 0.10)
colors_gfd <- c(
  "#08519c", # 极值：深海蓝 (< -0.1)
  "#3182bd", # 中值：中蓝 (-0.1 ~ -0.02)
  "#9ecae1", # 浅值：浅蓝 (-0.02 ~ -0.005)
  "#f7f7f7", # 零值缓冲带：纯白 (-0.005 ~ 0.005) - 过滤掉微小噪音！
  "#fcae91", # 浅值：浅红 (0.005 ~ 0.05)
  "#fb6a4a", # 中值：中红 (0.05 ~ 0.2)
  "#cb181d"# 极值：深红 (> 0.2 甚至到 0.8 的极端值都在这里)
)
# 2. 修复绘图函数：加上 ggplot2:: 强制调用排版边距
plot_delta <- function(data, var_filter, breaks_val, colors_val, title) {
  p <- ggplot(filter(data, Variable == var_filter)) +
    geom_sf(data = world_map, fill = "grey95", color = "white", size = 0.1) +
    geom_tile(aes(x = x, y = y, fill = Value)) +
    # 使用 stepsn 函数，将连续数据严格映射到你定义的离散区间
    scale_fill_stepsn(colors = colors_val, breaks = breaks_val,
                      limits = c(min(breaks_val)-0.1, max(breaks_val)+0.1),
                      na.value = "transparent", show.limits = TRUE) +
    facet_wrap(~ Scenario, ncol = 3) + # 一键双情景分屏
    labs(title = title, fill = "Change") +
    theme_void() +
    theme(
      # ✨ 核心修复：把所有的 margin() 改成 ggplot2::margin()
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 10)),
      strip.text = element_text(size = 14, face = "bold", margin = ggplot2::margin(b = 5)),
      legend.position = "bottom",
      legend.key.width = unit(2.5, "cm"), # 拉长图例，显示清晰的色阶
      legend.key.height = unit(0.5, "cm")
    ) +
    coord_sf(ylim = c(-60, 90), expand = FALSE)
  return(p)
}

# 现在可以安全地绘制并输出了！
p_delta_td <- plot_delta(df_delta_TD, "Delta_TD", breaks_td, colors_td, "Taxonomic Diversity (TD) Shift: 2100 vs present")
p_delta_gfd <- plot_delta(df_delta_GFD, "Delta_GFD", breaks_gfd, colors_gfd, "Guild Functional Menhinick (GFM) Shift: 2100 vs present")

p_delta_td
p_delta_gfd
# 拼图展示
combined_delta_maps <- p_delta_td / p_delta_gfd
ggsave("Global_DeltaTD_Distribution.pdf", p_delta_td, width = 12, height = 4, dpi = 300,bg="transparent")
ggsave("Global_DeltaGFD_Distribution.pdf", p_delta_gfd, width = 12, height = 4, dpi = 300,bg="transparent")


#loss gain area####

library(dplyr)
library(tidyr)
library(readr)
all_delta_data2<-all_delta_data1[,-c(6,11:14)]

all_delta_data3<-rbind(all_delta_data2,paddy_delta_data)
# ==========================================
# 1. 顶刊级面积计算：加入纬度余弦校正 (Cosine Correction)
# ==========================================
# 地球赤道上 1度 ≈ 111.32 km。
# 一个 0.1 x 0.1 度的网格，面积约为：(11.132 km) * (11.132 km * cos(纬度))
all_delta_data3 <- all_delta_data3 %>%#filter(Habitat!="Paddy soil") %>% 
  mutate(
    # 计算每个独立网格的真实物理面积 (km2)
    cell_area_km2 = (11.132) * (11.132 * cos(y * pi / 180)),
    
    # 严格根据你色盘的“白色缓冲带”定义真正的变化 (过滤掉极微小噪音)
    GFD_status = case_when(
      Delta_GFD < -0.02 ~ "Loss",
      Delta_GFD > 0.02  ~ "Gain",
      TRUE               ~ "Stable"
    ),
    TD_status = case_when(
      Delta_TD < -0.02 ~ "Loss",
      Delta_TD > 0.02 ~ "Gain",
      TRUE              ~ "Stable"
    )
  )

# ==========================================
# 2. Global Level (全球整体统计)
# ==========================================
global_stats <- all_delta_data %>%filter(Habitat!="Paddy soil") %>% 
  group_by(Scenario) %>%
  summarise(
    Total_Area_km2 = sum(cell_area_km2, na.rm = TRUE),
    
    # GFD (功能冗余) 统计
    GFD_Loss_Area = sum(cell_area_km2[GFD_status == "Loss"], na.rm = TRUE),
    GFD_Gain_Area = sum(cell_area_km2[GFD_status == "Gain"], na.rm = TRUE),
    GFD_Loss_Pct  = (GFD_Loss_Area / Total_Area_km2) * 100,
    GFD_Gain_Pct  = (GFD_Gain_Area / Total_Area_km2) * 100,
    GFD_Net_Pct   = GFD_Gain_Pct - GFD_Loss_Pct, # 净变化比例
    
    # TD (分类多样性) 统计
    TD_Loss_Area = sum(cell_area_km2[TD_status == "Loss"], na.rm = TRUE),
    TD_Gain_Area = sum(cell_area_km2[TD_status == "Gain"], na.rm = TRUE),
    TD_Loss_Pct  = (TD_Loss_Area / Total_Area_km2) * 100,
    TD_Gain_Pct  = (TD_Gain_Area / Total_Area_km2) * 100,
    TD_Net_Pct   = TD_Gain_Pct - TD_Loss_Pct
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 2))) # 保留两位小数，极其舒适

write_csv(global_stats, "Global_Area_Shift_Stats.csv")
print("=== Global Statistics ===")
print(global_stats)

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(patchwork)

# ==========================================
# 1. 数据重组 (宽表转长表，为 dodged 柱状图准备)
# ==========================================
plot_data.change <- global_stats %>%
  dplyr::select(Scenario, matches("Pct")) %>%
  # 将列名拆分为 Variable (GFD/TD) 和 Metric (Loss/Gain/Net)
  pivot_longer(
    cols = -Scenario,
    names_to = c("Variable", "Metric"),
    names_pattern = "^(.*)_(.*)_Pct$"
  ) %>%
  mutate(
    # 将 GFD 重命名为 GFM
    Variable = ifelse(Variable == "GFD", "GFM", Variable),
    
    # ✨ 核心：将 Loss 的数值转为负数，确保它画在 X 轴左侧
    Plot_Value = ifelse(Metric == "Loss", -value, value),
    
    # 清洗情景名称，去掉多余的字符串
    Scenario_Clean = str_replace(Scenario, "Delta_", ""),
    Scenario_Clean = str_replace(Scenario_Clean, "_vs_Present", ""),
    
    # 设定因子的固定顺序：TD 在上，GFM 在下；指标顺序：Loss, Net, Gain
    Variable = factor(Variable, levels = c("TD", "GFM")),
    Metric = factor(Metric, levels = c("Loss", "Net", "Gain"))
  )

# ==========================================
# 2. 定义顶刊级配色方案 (按指标着色)
# ==========================================
# Loss: 警告红 | Net: 沉稳灰黑 | Gain: 生态绿
metric_colors <- c(
  "Loss" = "#5E6C82", 
  "Net"  = "#F5E4C8", 
  "Gain" = "#81B3A9"
)
# ==========================================
# 3. 封装单张图绘制函数 (包含动态数值标签)
# ==========================================
plot_scenario <- function(data, scenario_name, show_x_label = FALSE) {
  
  x_title <- if(show_x_label) "Proportion of Global Area (%) \n ← Loss Area                                  Gain Area →" else NULL
  
  df_sub <- filter(data, Scenario_Clean == scenario_name)
  
  p <- ggplot(df_sub, aes(y = Variable, x = Plot_Value, fill = Metric)) +
    
    # 柱形图主体
    geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 1, linewidth = 0.3) +
    
    # ✨ 核心新增：在柱子外侧打上数值标签 (带 % 符号)
    geom_text(
      aes(
        # 保留 1 位小数，并加上 %
        label = paste0(round(Plot_Value, 1), "%"), 
        # 魔法对齐：正数往右偏 (hjust=-0.1)，负数往左偏 (hjust=1.1)
        hjust = ifelse(Plot_Value >= 0, -0.1, 1.1)
      ),
      position = position_dodge(width = 0.8), # 必须和柱子的 dodge 宽度完全一致！
      size = 3.5,
      color = "black",
      fontface = "bold"
    ) +
    
    # 0 基准线 (置于顶层)
    geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
    
    # ✨ 核心修改：恢复负数显示，并扩大 X 轴界限防止文字被切掉
    scale_x_continuous(
      limits = c(-40, 40),  # 必须适度拉宽 limits，给两边的文字标签留出物理空间
      breaks = seq(-40, 40, by = 10),
      labels = function(x) paste0(x, "%") # 直接显示原本的 x，保留负号
    ) +
    
    scale_fill_manual(values = metric_colors) +
    
    labs(
     # title = paste("Scenario:", scenario_name),
      x = x_title,
      y = NULL
    ) +
    
    theme_classic(base_size = 14) +
    theme(aspect.ratio = 0.6,
      plot.title = element_text(face = "bold", size = 12, hjust = 0),
      axis.title.x = element_text(face = "bold", size = 12, margin = ggplot2::margin(t = 10)),
      axis.text.y = element_text(face = "bold", size = 12, color = "black"),
      axis.text.x = element_text(color = "black", size = 11),
      
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank()
      
     # panel.grid.major.x = element_line(color = "grey85", linetype = "dashed")
    )
  
  return(p)
}

# ==========================================
# 4. 生成单独的情景图并用 patchwork 纵向拼接
# ==========================================
p_ssp126 <- plot_scenario(plot_data.change, "SSP126", show_x_label = FALSE)
p_ssp245 <- plot_scenario(plot_data.change, "SSP245", show_x_label = FALSE)
p_ssp370 <- plot_scenario(plot_data.change, "SSP370", show_x_label = FALSE)
p_ssp126
final_plot <- (p_ssp126 / p_ssp245 / p_ssp370) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "top", 
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold"))
ggsave("Global_Shift_Bars_p_ssp126.pdf",p_ssp126 , width = 6, height = 7, dpi = 300,bg="transparent")
ggsave("Global_Shift_Bars_p_ssp245.pdf",p_ssp245 , width = 6, height = 7, dpi = 300,bg="transparent")
ggsave("Global_Shift_Bars_p_ssp370.pdf",p_ssp370 , width = 6, height = 7, dpi = 300,bg="transparent")

#气泡图###
plot_scenario_lollipop <- function(data, scenario_name, show_x_label = FALSE) {
  
  x_title <- if(show_x_label) "Proportion of Global Area (%)" else NULL
  df_sub <- filter(data, Scenario_Clean == scenario_name)
  
  # 🌟 关键：统一的 dodge 参数
  pd <- position_dodge(width = 0.7) 
  
  p <- ggplot(df_sub, aes(y = Variable, x = Plot_Value, color = Metric, group = Metric)) +
    
    # 🌟 1. 用 geom_linerange 替代 geom_segment
    # xmin=0 确保从中心线开始，xmax=Plot_Value 连到球
    geom_linerange(
      aes(xmin = 0, xmax = Plot_Value),
      position = pd,
      linewidth = 1,
      alpha = 0.6
    ) +
    
    # 🌟 2. 气泡图层 (必须与 pd 保持一致)
    geom_point(
      aes(size=abs(Plot_Value),fill = Metric), 
      position = pd,
      # size = 10, 
      shape = 21, 
      stroke = 1.2, 
      color = "white"
    ) +
    
    # 🌟 3. 文字标签 (水平位移用 hjust，垂直位移用 pd)
    geom_text(
      aes(
        label = paste0(round(abs(Plot_Value), 1), "%"),
        hjust = ifelse(Plot_Value >= 0, -0.6, 1.6)
      ),
      position = pd,
      size = 3.5,
      color = "black",
      fontface = "bold"
    ) +scale_size_continuous(range = c(4, 10)) +
    
    # 辅助线：0 基准轴
    geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
    
    # X 轴美化
    scale_x_continuous(
      limits = c(-45, 45), 
      breaks = seq(-40, 40, by = 10),
      labels = function(x) paste0(abs(x), "%")
    ) +
    
    # 配色
    scale_color_manual(values = metric_colors) +
    scale_fill_manual(values = metric_colors) +
    
    labs(x = x_title, y = NULL) +
    
    # 绘图风格
    theme_classic(base_size = 14) +
    theme(
      aspect.ratio = 0.5,
      legend.position = "NULL",
      axis.text.y = element_text(face = "bold", size = 16, color = "black"),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.x = element_blank()
    )
  
  return(p)
}
#==========================================
# 4. 生成单独的情景图并用 patchwork 纵向拼接
# ==========================================
p_ssp126 <- plot_scenario_lollipop(plot_data.change, "SSP126", show_x_label = FALSE)
p_ssp245 <- plot_scenario_lollipop(plot_data.change, "SSP245", show_x_label = FALSE)
p_ssp370 <- plot_scenario_lollipop(plot_data.change, "SSP370", show_x_label = FALSE)
p_ssp126
final_plot <- (p_ssp126 / p_ssp245 / p_ssp370) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "top", 
        legend.title = element_blank(),
        legend.text = element_text(size = 12, face = "bold"))
ggsave("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/ok/5/Global_Shift_Bars_p_ssp126.pop.pdf",p_ssp126 , width = 6, height = 10, dpi = 300,bg="transparent")
ggsave("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/ok/5/Global_Shift_Bars_p_ssp245.pop.pdf",p_ssp245 , width = 6, height = 10, dpi = 300,bg="transparent")
ggsave("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/ok/5/Global_Shift_Bars_p_ssp370.pop.pdf",p_ssp370 , width = 6, height = 10, dpi = 300,bg="transparent")

# 打印最终效果
print(final_plot)

# 保存高分辨率 PDF (如果文字太挤，可以适当增加 width 宽度)
# ggsave("Global_Shift_Bars_Stacked_Labeled.pdf", final_plot, width = 9, height = 7, dpi = 300)


# 保存高分辨率 PDF
# ggsave("Global_Shift_Bars_Stacked.pdf", final_plot, width = 8, height = 7, dpi = 300)
# 保存高分辨率 PDF
# ggsave("Global_Loss_Gain_Net_Diverging_Bars.pdf", p_bars, width = 8, height = 4, dpi = 300)
# ==========================================
# 3. Habitat Level (不同生态系统级别统计)
# ==========================================
habitat_stats1<- all_delta_data3 %>%
  # 过滤掉无分类的杂讯数据
  filter(!Habitat %in% c("Nodata", "Unclassified", "80", "Snow and ice")) %>%
  group_by(Scenario, Habitat) %>%
  summarise(
    Total_Area_km2 = sum(cell_area_km2, na.rm = TRUE),
    
    GFD_Loss_Pct  = (sum(cell_area_km2[GFD_status == "Loss"], na.rm = TRUE) / Total_Area_km2) * 100,
    GFD_Gain_Pct  = (sum(cell_area_km2[GFD_status == "Gain"], na.rm = TRUE) / Total_Area_km2) * 100,
    GFD_Net_Pct   = GFD_Gain_Pct - GFD_Loss_Pct,
    
    TD_Loss_Pct  = (sum(cell_area_km2[TD_status == "Loss"], na.rm = TRUE) / Total_Area_km2) * 100,
    TD_Gain_Pct  = (sum(cell_area_km2[TD_status == "Gain"], na.rm = TRUE) / Total_Area_km2) * 100,
    TD_Net_Pct   = TD_Gain_Pct - TD_Loss_Pct,
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ round(., 2))) %>%
  arrange(Scenario, desc(Total_Area_km2))

write_csv(habitat_stats1, "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/submit/data/Habitat_Area_Shift_Stats_paddyall.csv")
print("=== Habitat Statistics ===")
print(head(habitat_stats, 10))

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(grDevices) # 用于输出多页 PDF

# ==========================================
# 1. 数据重组 (宽表转长表，加入 Habitat 维度)
# ==========================================
# 注意：这里使用你刚生成的 habitat_stats1
plot_data_hab <- habitat_stats1 %>%
  dplyr::select(Scenario, Habitat, matches("Pct")) %>%
  pivot_longer(
    cols = matches("Pct"),
    names_to = c("Variable", "Metric"),
    names_pattern = "^(.*)_(.*)_Pct$"
  ) %>%
  mutate(
    Variable = ifelse(Variable == "GFD", "GFM", Variable),
    
    # Loss 转负数，用于 X 轴左侧绘制
    Plot_Value = ifelse(Metric == "Loss", -value, value),
    
    # 清洗情景名称
    Scenario_Clean = str_replace(Scenario, "Delta_", ""),
    Scenario_Clean = str_replace(Scenario_Clean, "_vs_Present", ""),
    Scenario_Clean = factor(Scenario_Clean, levels = c("SSP126", "SSP245", "SSP370")), # 让分面从上到下按情景排序
    
    Variable = factor(Variable, levels = c("TD", "GFM")),
    Metric = factor(Metric, levels = c("Loss", "Net", "Gain"))
  )

# ==========================================
# 2. 定义配色方案
# ==========================================

# ==========================================
# 3. 封装单张 Habitat 绘图函数 (采用 facet_wrap 自动堆叠 3 个情景)
# ==========================================
plot_single_habitat <- function(hab_name, data) {
  
  df_sub <- data %>% filter(Habitat == hab_name)
  
  p <- ggplot(df_sub, aes(y = Variable, x = Plot_Value, fill = Metric)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 1, linewidth = 0.3) +
    
    # 数值标签 (动态偏移)
    geom_text(
      aes(
        label = paste0(round(Plot_Value, 1), "%"), 
        hjust = ifelse(Plot_Value >= 0, -0.1, 1.1)
      ),
      position = position_dodge(width = 0.8), 
      size = 3.2, color = "black", fontface = "bold"
    ) +
    
    # 0 基准线
    geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
    
    # X 轴刻度：所有生态系统统一定死 (-100 到 100)，这样翻看多页 PDF 时，不同系统间的变化幅度才有视觉可比性！
    scale_x_continuous(
      limits = c(-80, 80), 
      breaks = seq(-80, 80, by = 20),
      labels = function(x) paste0(abs(x), "%")
    ) +
    
    scale_fill_manual(values = metric_colors) +
    
    # ✨ 核心魔法：直接使用 facet_wrap 把 3 个情景分成 3 行，极其整洁
    facet_wrap(~ Scenario_Clean, ncol = 1) +
    
    labs(
      title = paste("Ecosystem:", hab_name),
      x = "Proportion of Habitat Area (%) \n ← Loss Area                                  Gain Area →",
      y = NULL
    ) +
    
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = ggplot2::margin(b = 15)),
      strip.text = element_text(face = "bold", size = 12, color = "white"),
      strip.background = element_rect(fill = "grey30", color = "black"), # 分面标题带深灰底色，高级感拉满
      
      axis.title.x = element_text(face = "bold", size = 12, margin = ggplot2::margin(t = 10)),
      axis.text.y = element_text(face = "bold", size = 12, color = "black"),
      axis.text.x = element_text(color = "black", size = 11),
      
      panel.grid.major= element_blank(),
      panel.grid.minor = element_blank(),
      
      legend.position = "top", 
      legend.title = element_blank(),
      legend.text = element_text(size = 12, face = "bold")
    )
  
  return(p)
}

# ==========================================
# 4. 自动化批量输出多页 PDF (Supplementary 材料神器)
# ==========================================
# 提取所有不重复的生态系统名称
habitat_list <- unique(plot_data_hab$Habitat)

# 打开一个 PDF 设备
pdf("./change/Supplementary_Ecosystem_Shift_Profiles.pdf", width = 6, height = 8)

# 循环遍历每个生态系统，画图并写入 PDF 的新一页
for (hab in habitat_list) {
  p_hab <- plot_single_habitat(hab, plot_data_hab)
  print(p_hab)
  cat("已成功生成图像：", hab, "\n")
}

# 关闭并保存 PDF
dev.off()

 cat("\n✅ 所有生态系统图像已批量打包至 'Supplementary_Ecosystem_Shift_Profiles.pdf'，请在工作目录查收！\n")

#density lat纬度####
library(dplyr)
library(tidyr)
library(ggplot2)

# ==========================================
# 1. 数据预处理：长宽转换与按生态系统聚合
# ==========================================
# 将 all_delta_data (宽表) 转换为绘图所需的 (长表)
df_delta_long <- all_delta_data %>%
  pivot_longer(
    cols = c("Delta_TD", "Delta_GFD"),
    names_to = "Variable",
    values_to = "Value"
  )

# 按纬度带（每1度）、情景、生态系统和变量计算 Mean 和 SD
lat_summary_habitat <- df_delta_long %>%
  mutate(lat_band = round(y)) %>% # 划定 1 度的纬度带
  group_by(Habitat, Scenario, Variable, lat_band) %>%
  summarise(
    mean_val = mean(Value, na.rm = TRUE),
    sd_val = sd(Value, na.rm = TRUE),
    .groups = "drop"
  )

# 提取所有有效的生态系统名称 (剔除无效分类)
valid_habitats <- unique(lat_summary_habitat$Habitat)
valid_habitats <- valid_habitats[!valid_habitats %in% c("Nodata", "Unclassified", "80")]

# ==========================================
# 2. 开启批量生成循环：为每个生态系统绘图
# ==========================================
for (hab in valid_habitats) {
  
  # 提取当前生态系统的数据
  df_hab <- lat_summary_habitat %>% filter(Habitat == hab)
  
  # 如果该生态系统数据为空，跳过
  if(nrow(df_hab) == 0) next
  
  # ✨ 核心魔法 1：动态计算专属比例系数 (Scaling Coefficient)
  max_td <- max(abs(df_hab$mean_val[df_hab$Variable == "Delta_TD"]), na.rm = TRUE)
  max_gfd <- max(abs(df_hab$mean_val[df_hab$Variable == "Delta_GFD"]), na.rm = TRUE)
  
  # 防护机制：如果全是 NA 或 0，系数设为 1
  if (is.infinite(max_td) | is.infinite(max_gfd) | is.na(max_td) | is.na(max_gfd) | max_gfd == 0) {
    coeff <- 1 
  } else {
    coeff <- max_td / max_gfd
  }
  
  # ✨ 核心魔法 2：根据专属比例缩放 GFD 数据
  df_hab_scaled <- df_hab %>%
    mutate(
      plot_mean = ifelse(Variable == "Delta_GFD", mean_val * coeff, mean_val),
      plot_sd   = ifelse(Variable == "Delta_GFD", sd_val * coeff, sd_val)
    )
  
  # ==========================================
  # 3. 绘制专属双 X 轴图表
  # ==========================================
  p_hab <- ggplot(df_hab_scaled, aes(y = lat_band, x = plot_mean, color = Variable, fill = Variable)) +
    # 0 基准虚线
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 0.6) +
    
    # 阴影与核心曲线
    geom_ribbon(aes(xmin = plot_mean - plot_sd, xmax = plot_mean + plot_sd), alpha = 0.3, color = NA) +
    geom_path(linewidth = 0.8) +
    
    # 按情景分面
    facet_wrap(~ Scenario, ncol = 3) + 
    scale_y_continuous(limits = c(-60, 90), breaks = seq(-60, 90, 30), expand = c(0,0)) +
    
    # ✨ 核心魔法 3：应用专属的双 X 轴
    scale_x_continuous(
      name = "Mean Δ TD", 
      sec.axis = sec_axis(~ . / coeff, name = "Mean Δ GFM") 
    ) +
    ylab("Latitude (°N)") +
    labs(title = paste0("Latitudinal Shift in ", hab)) + # 加上生态系统大标题
    
    theme_bw() +
    theme(
      aspect.ratio = 1.2, # 稍微拉长一点以适应纬度分布
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 15)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 11, face = "bold", margin = ggplot2::margin(b = 5, t = 5)),
      axis.title.y = element_text(size = 12, face = "bold", color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      
      # 坐标轴颜色匹配
      axis.title.x.bottom = element_text(size = 12, face = "bold", color = "#059f86", margin = ggplot2::margin(t = 5)),
      axis.text.x.bottom  = element_text(size = 10, color = "#059f86", face = "bold"),
      
      axis.title.x.top    = element_text(size = 12, face = "bold", color = "#08519c", margin = ggplot2::margin(b = 5)),
      axis.text.x.top     = element_text(size = 10, color = "#08519c", face = "bold"),
      
      # 图例悬浮在右上角
      legend.position = c(0.98, 0.98),        
      legend.justification = c(1, 1),         
      legend.background = element_rect(fill = alpha("white", 0.8), color = "grey80", linewidth = 0.5), 
      legend.title = element_blank(),        
      legend.text = element_text(size = 9, face = "bold"),
      legend.key = element_blank()
    ) +
    
    scale_color_manual(
      values = c("Delta_TD" = "#059f86", "Delta_GFD" = "#08519c"),
      labels = c("Delta_TD" = "Δ TD", "Delta_GFD" = "Δ GFD")
    ) +
    scale_fill_manual(
      values = c("Delta_TD" = "#91d1c2", "Delta_GFD" = "#9ecae1"),
      labels = c("Delta_TD" = "Δ TD", "Delta_GFD" = "Δ GFD")
    )
  
  # ==========================================
  # 4. 自动保存每张图
  # ==========================================
  # 生成安全的文件名 (去掉可能引起路径错误的空格或特殊字符)
  safe_hab_name <- gsub(" & | ", "_", hab)
  file_name <- paste0("Lat_Curve_", safe_hab_name, ".pdf")
  
  ggsave(file_name, p_hab, width = 12, height = 6, dpi = 300,bg="transparent")
  
  # 在控制台打印进度
  cat("Saved:", file_name, "\n")
}
#不分系统，按照情景####
lat_summary_all <- df_delta_long %>%
  mutate(lat_band = round(y)) %>% # 划定 1 度的纬度带
  group_by(Scenario, Variable, lat_band) %>%
  summarise(
    mean_val = mean(Value, na.rm = TRUE),
    sd_val = sd(Value, na.rm = TRUE),
    .groups = "drop"
  )
max_td <- max(abs(lat_summary_all$mean_val[lat_summary_all$Variable == "Delta_TD"]), na.rm = TRUE)

max_gfd <- max(abs(lat_summary_all$mean_val[lat_summary_all$Variable == "Delta_GFD"]), na.rm = TRUE)

coeff <- max_td / max_gfd  # 自动计算出的放大倍数 (大概在 10~15 倍左右)



# ==========================================

# 2. 缩放 GFD 的绘图数据
# ==========================================
lat_summary_scaled <- lat_summary_all %>%
  mutate(
    # 如果是 GFD，就把它的 mean 和 sd 都乘以 coeff 放大；如果是 TD 就保持不变
    plot_mean = ifelse(Variable == "Delta_GFD", mean_val * coeff, mean_val),
    plot_sd   = ifelse(Variable == "Delta_GFD", sd_val * coeff, sd_val)
  )
# ==========================================

# 1. 计算真实的全局 Mean 和 SD (绝不含缩放)
# ==========================================
anno_df <- df_delta_long %>%
  group_by(Scenario, Variable) %>%
  summarise(
    real_mean = mean(Value, na.rm = TRUE),
    real_sd   = sd(Value, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  mutate(
    # 拼接极其学术的高级文本
    label_text = sprintf("%s: %.3f ± %.3f", 
                         ifelse(Variable == "Delta_TD", "Δ TD", "Δ GFM"), 
                         real_mean, real_sd),
    
    # 摆放位置：Y轴 (纬度) 放在南半球底部
    y_pos = ifelse(Variable == "Delta_TD", -42, -55)
  )
# 在控制台查看真实的数值
print(anno_df)
# ==========================================

# 3. 绘制双 X 轴高级图表

# ==========================================
p_lat_dual <- ggplot(lat_summary_scaled, aes(y = lat_band, x = plot_mean, color = Variable, fill = Variable)) +
  # 0 基准虚线
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 0.6) +
  # 阴影与核心曲线
  geom_ribbon(aes(xmin = plot_mean - plot_sd, xmax = plot_mean + plot_sd), alpha = 0.3, color = NA) +
  geom_path(linewidth = 0.8) +
  # ✨ 核心魔法图层：将 Mean ± SD 文本印在图上
  geom_text(
    data = anno_df,
    aes(y = y_pos, label = label_text, color = Variable), 
    x = Inf,                 # X轴钉死在最右侧
    hjust = 1.05,            # 略微向左推一点点，防止贴着边框
    inherit.aes = FALSE,     # 不继承主图的 X/Y 映射，避免报错
    size = 3.5, 
    fontface = "bold",
    show.legend = FALSE      # 文本图层不需要额外生成图例
  ) +
  facet_wrap(~ Scenario, ncol = 3) +
  scale_y_continuous(limits = c(-60, 90), breaks = seq(-60, 90, 30), expand = c(0,0)) +
  # ✨ 核心魔法：双 X 轴设置！
  scale_x_continuous(
    name = "Mean Δ TD", # 底部的主 X 轴 (给 TD 用)
    # 顶部的副 X 轴 (给 GFD 用)，刻度自动除以 coeff 还原为真实数值！
    sec.axis = sec_axis(~ . / coeff, name = "Mean Δ GFM") 
  ) +ylab("Latitude (°N)")+
  theme_bw() +
  theme(aspect.ratio = 0.6,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 12, face = "bold", margin = ggplot2::margin(b = 5, t = 5)),
        axis.title.y = element_text(size = 12, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        # ✨ 色彩匹配魔法：让坐标轴的颜色和曲线的颜色完全一致！
        axis.title.x.bottom = element_text(size = 12, face = "bold", color = "#059f86", margin = ggplot2::margin(t = 5)),
        axis.text.x.bottom  = element_text(size = 10, color = "#059f86", face = "bold"),
     axis.title.x.top    = element_text(size = 12, face = "bold", color = "#08519c", margin = ggplot2::margin(b = 5)),
        axis.text.x.top     = element_text(size = 10, color = "#08519c", face = "bold"),
        # 图例设置 (悬浮在右上角)
        legend.position = c(0.98, 0.98),       
        
        legend.justification = c(1, 1),        
        
        legend.background = element_rect(fill = alpha("white", 0.8), color = "grey80", size = 0.5), 
        
        legend.title = element_blank(),        
        
        legend.text = element_text(size = 9, face = "bold"),
        
        legend.key = element_blank()
        
  ) +
  scale_color_manual(
    values = c("Delta_TD" = "#059f86", "Delta_GFD" = "#08519c"),
    labels = c("Delta_TD" = "Δ TD", "Delta_GFD" = "Δ GFM") # 优化图例显示的文字
  ) +
  scale_fill_manual(
    values = c("Delta_TD" = "#91d1c2", "Delta_GFD" = "#9ecae1"),
    labels = c("Delta_TD" = "Δ TD", "Delta_GFD" = "Δ GFM")
  )
# 预览并保存
p_lat_dual
ggsave("lat_curve_Globalssp_Delta_Distribution.pdf",p_lat_dual, width = 12, height = 4, dpi = 300,bg="transparent")


#lat 不分面####
library(dplyr)
library(tidyr)
library(ggplot2)

# ==========================================
# 1. 数据预处理：长宽转换与聚合 (保持不变)
# ==========================================
df_delta_long <- all_delta_data %>%
  pivot_longer(
    cols = c("Delta_TD", "Delta_GFD"),
    names_to = "Variable",
    values_to = "Value"
  )

lat_summary_habitat <- df_delta_long %>%
  mutate(lat_band = round(y)) %>% 
  group_by(Habitat, Scenario, Variable, lat_band) %>%
  summarise(
    mean_val = mean(Value, na.rm = TRUE),
    sd_val = sd(Value, na.rm = TRUE),
    .groups = "drop"
  )

valid_habitats <- unique(lat_summary_habitat$Habitat)
valid_habitats <- valid_habitats[!valid_habitats %in% c("Nodata", "Unclassified", "80")]

# 提取所有有效的情景 (Scenario)
valid_scenarios <- unique(lat_summary_habitat$Scenario)

# ==========================================
# 2. 开启双重批量生成循环 (Habitat x Scenario)
# ==========================================
for (hab in valid_habitats) {
  for (scen in valid_scenarios) {
    
    # 提取当前生态系统 + 当前情景的数据
    df_plot <- lat_summary_habitat %>% filter(Habitat == hab, Scenario == scen)
    
    # 如果没数据，直接跳过
    if(nrow(df_plot) == 0) next
    
    # ✨ 核心魔法 1：针对当前这张唯一的图，计算最完美的比例系数
    max_td <- max(abs(df_plot$mean_val[df_plot$Variable == "Delta_TD"]), na.rm = TRUE)
    max_gfd <- max(abs(df_plot$mean_val[df_plot$Variable == "Delta_GFD"]), na.rm = TRUE)
    
    if (is.infinite(max_td) | is.infinite(max_gfd) | is.na(max_td) | is.na(max_gfd) | max_gfd == 0) {
      coeff <- 1 
    } else {
      coeff <- max_td / max_gfd
    }
    
    # 根据比例缩放 GFD 数据
    df_plot_scaled <- df_plot %>%
      mutate(
        plot_mean = ifelse(Variable == "Delta_GFD", mean_val * coeff, mean_val),
        plot_sd   = ifelse(Variable == "Delta_GFD", sd_val * coeff, sd_val)
      )
    
    # ==========================================
    # 3. 绘制单张独立图表 (移除了 facet_wrap)
    # ==========================================
    p_single <- ggplot(df_plot_scaled, aes(y = lat_band, x = plot_mean, color = Variable, fill = Variable)) +
      geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 0.6) +
      geom_ribbon(aes(xmin = plot_mean - plot_sd, xmax = plot_mean + plot_sd), alpha = 0.3, color = NA) +
      geom_path(linewidth = 0.8) +
      
      scale_y_continuous(limits = c(-60, 90), breaks = seq(-60, 90, 30), expand = c(0,0)) +
      
      scale_x_continuous(
        name = "Mean Δ TD", 
        sec.axis = sec_axis(~ . / coeff, name = "Mean Δ GFM") 
      ) +
      ylab("Latitude (°N)") +
      
      # ✨ 核心修改：标题现在包含 生态系统 + 情景
    #  labs(title = paste0("Latitudinal Shift in ", hab, "\n(", scen, ")")) + 
      
      theme_bw() +
      theme(
        aspect.ratio = 1.5, # 单图没必要拉得太长，0.8 看起来更饱满
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 15), lineheight = 1.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        axis.title.y = element_text(size = 12, face = "bold", color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        
        axis.title.x.bottom = element_text(size = 12, face = "bold", color = "#059f86", margin = ggplot2::margin(t = 5)),
        axis.text.x.bottom  = element_text(size = 10, color = "#059f86", face = "bold"),
        
        axis.title.x.top    = element_text(size = 12, face = "bold", color = "#08519c", margin = ggplot2::margin(b = 5)),
        axis.text.x.top     = element_text(size = 10, color = "#08519c", face = "bold"),
        
        legend.position = c(0.98, 0.98),        
        legend.justification = c(1, 1),         
        legend.background = element_rect(fill = alpha("white", 0.8), color = "grey80", linewidth = 0.5), 
        legend.title = element_blank(),        
        legend.text = element_text(size = 9, face = "bold"),
        legend.key = element_blank()
      ) +
      scale_color_manual(
        values = c("Delta_TD" = "#059f86", "Delta_GFD" = "#08519c"),
        labels = c("Delta_TD" = "Δ TD", "Delta_GFD" = "Δ GFD")
      ) +
      scale_fill_manual(
        values = c("Delta_TD" = "#91d1c2", "Delta_GFD" = "#9ecae1"),
        labels = c("Delta_TD" = "Δ TD", "Delta_GFD" = "Δ GFD")
      )
    
    # ==========================================
    # 4. 自动生成双重命名文件并保存
    # ==========================================
    # 清洗掉可能导致报错的特殊符号 (如空格, &, 减号)
    safe_hab <- gsub(" & | ", "_", hab)
    safe_scen <- gsub(" |-|\\.", "_", scen) 
    
    file_name <- paste0("Lat_Curve_", safe_hab, "_", safe_scen, ".svg")
    
    # 由于是单图，宽度可以稍微窄一点
    ggsave(file_name, p_single, width = 4, height = 5, dpi = 300,bg="transparent")
    
    cat("Saved:", file_name, "\n")
  }
}


