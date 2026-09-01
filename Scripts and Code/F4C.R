library(dplyr)
library(purrr)
library(ggplot2)
library(minpack.lm)
library(ggrepel) # 🌟 确保引入智能标签包

# ==========================================
# 1. 数据准备
# ==========================================
p.guild1 <- shannon_guild4 %>% 
  dplyr::select("Menhinick_Coef", "Shannon.TD", "Ecosystem.y") %>% 
  na.omit() %>%
  dplyr::filter(Shannon.TD > 0)

p.guild1$Ecosystem.y <- gsub("Artificial surfaces", "Artificial Surfaces", p.guild1$Ecosystem.y)

# ==========================================
# 2. 数据装箱 (Binning) 
# ==========================================
bin_width <- 0.4
binned_data <- p.guild1 %>%
  dplyr::mutate(
    bin_x = floor(Shannon.TD / bin_width) * bin_width + (bin_width / 2)
  ) %>%
  dplyr::group_by(Ecosystem.y, bin_x) %>%
  dplyr::summarise(
    bin_y = mean(Menhinick_Coef, na.rm = TRUE),
    n_points = n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_points >= 3) 

# ==========================================
# 3. 重新拟合，并计算极值 Tipping Point
# ==========================================
eq_data_binned <- binned_data %>%
  dplyr::group_by(Ecosystem.y) %>% 
  dplyr::filter(!Ecosystem.y %in% c("Shrubland", "Bare Land")) %>% 
  dplyr::group_modify(~ {
    
    b_start <- 0.5
    a_start <- max(.x$bin_y, na.rm = TRUE) * b_start * exp(1)
    
    model <- tryCatch({
      nlsLM(
        bin_y ~ a * bin_x * exp(-b * bin_x), 
        data = .x,
        start = list(a = a_start, b = b_start),
        lower = c(a = 0, b = 0.0001),        
        upper = c(a = Inf, b = Inf), 
        control = nls.lm.control(maxiter = 1000)
      )
    }, error = function(e) return(NULL))
    
    if (is.null(model)) {
      return(tibble(eq_text = "Fit failed", a=NA, b=NA, pval=NA, r2=NA, Tipping_X=NA, Tipping_Y=NA, point_label=NA))
    }
    
    coefs <- coef(model)
    a_val <- coefs["a"]
    b_val <- coefs["b"]
    
    coef_summary <- summary(model)$coefficients
    col_idx <- if(ncol(coef_summary) >= 4) 4 else ncol(coef_summary)
    pval <- tryCatch({ coef_summary["b", col_idx] }, error = function(e) NA_real_)
    p_text <- if (is.na(pval)) "p = NA" else if (pval < 0.001) "p < 0.001" else sprintf("p = %.3f", pval)
    
    rss <- sum(residuals(model)^2)
    tss <- sum((.x$bin_y - mean(.x$bin_y))^2)
    pseudo_r2 <- ifelse(tss == 0, NA, 1 - (rss / tss))
    if (!is.na(pseudo_r2) && pseudo_r2 < 0) pseudo_r2 <- 0.01 
    
    r2_text <- if (is.na(pseudo_r2)) "R² = NA" else sprintf("R² = %.2f", pseudo_r2)
    
    eq_text <- sprintf(
      "%s: y = %.2f*x*e^(-%.2f*x), %s, %s",
      .y$Ecosystem.y, a_val, b_val, r2_text, p_text
    )
    
    # 🌟【核心新增】：微积分精确解求 Ricker 模型 Y 最大值时的 Tipping Point 坐标
    tipping_x_val <- 1 / b_val
    tipping_y_val <- (a_val / b_val) * exp(-1)
    
    # 格式化图中钻石点旁边的文字标签
    point_label_text <- sprintf("Max\n(%.2f, %.2f)", tipping_x_val, tipping_y_val)
    
    tibble(
      eq_text = eq_text, a = a_val, b = b_val, pval = pval, r2 = pseudo_r2,
      Tipping_X = tipping_x_val, Tipping_Y = tipping_y_val, point_label = point_label_text
    )
  }) %>%
  ungroup()

# ==========================================
# 4. 绘图准备
# ==========================================
p.guild2 <- p.guild1 %>%
  dplyr::filter(!Ecosystem.y %in% c("Shrubland", "Bare Land")) 

eco_levels <- c("Agricultural Land", "Grassland", "Artificial Surfaces","Tundra", "Forest", "Wetland", "Paddy soil")
p.guild2$Ecosystem.y <- factor(p.guild2$Ecosystem.y, levels = eco_levels)
binned_data$Ecosystem.y <- factor(binned_data$Ecosystem.y, levels = eco_levels)
eq_data_binned$Ecosystem.y <- factor(eq_data_binned$Ecosystem.y, levels = eco_levels)

global_b_start <- 0.5
global_a_start <- max(binned_data$bin_y, na.rm = TRUE) * global_b_start * exp(1)

# ==========================================
# 5. 绘图展现
# ==========================================
p.guild_micha_binned <- ggplot() +
  # 图层 1: 原始蓝色系透明散点打底
  geom_point(
    data = p.guild2,
    aes(x = Shannon.TD, y = Menhinick_Coef, color = Ecosystem.y),
    alpha = 0.15, size = 1.0
  ) +
  
  # 图层 2: 装箱均值黑圈点
  geom_point(
    data = binned_data %>% filter(!is.na(Ecosystem.y)),
    aes(x = bin_x, y = bin_y, fill = Ecosystem.y),
    shape = 21, color = "black", size = 1, stroke = 0, alpha = 0.6
  ) +
  
  # 图层 3: geom_smooth 拟合 Ricker 曲线
  geom_smooth(
    data = binned_data %>% filter(!is.na(Ecosystem.y)),
    aes(x = bin_x, y = bin_y, color = Ecosystem.y),
    method = "nls",
    formula = y ~ a * x * exp(-b * x), 
    method.args = list(
      start = list(a = global_a_start, b = global_b_start),
      algorithm = "port",
      lower = c(a = 0, b = 0.0001),     
      upper = c(a = Inf, b = Inf), 
      control = nls.control(maxiter = 1000, warnOnly = TRUE) 
    ),
    se = FALSE,
    fullrange = TRUE,  
    linewidth = 0.3
  ) +
  
  # 🌟【新图层 A】：从波峰精准垂直下落到 X 轴的虚线
  geom_segment(
    data = eq_data_binned %>% filter(!is.na(Tipping_X)),
    aes(x = Tipping_X, xend = Tipping_X, y = 0, yend = Tipping_Y, color = Ecosystem.y),
    linetype = "dashed", alpha = 0.6, linewidth = 0.2
  ) +
  
  # 🌟【新图层 B】：死死钉在曲线最高峰顶部的钻石标记（Tipping Point）
  geom_point(
    data = eq_data_binned %>% filter(!is.na(Tipping_X)),
    aes(x = Tipping_X, y = Tipping_Y, fill = Ecosystem.y),
    shape = 23, color = "black", size = 2.5, stroke = 0, alpha = 0.9
  ) +
  
  # 🌟【新图层 C】：波峰坐标的智能避让文字标注
  ggrepel::geom_text_repel(
    data = eq_data_binned %>% filter(!is.na(Tipping_X)),
    aes(x = Tipping_X, y = Tipping_Y, label = point_label, color = Ecosystem.y),
    size = 2.2, fontface = "bold", lineheight = 0.9,
    box_padding = 0.4, point_padding = 0.3,
    force = 2, max.overlaps = Inf, segment.color = "gray60", segment.linewidth = 0.3
  ) +
  
  facet_wrap(~Ecosystem.y, ncol = 7, scales = "free_x") +
  
  geom_text(
    data = eq_data_binned,
    aes(x = Inf, y = Inf, label = eq_text),
    inherit.aes = FALSE,
    hjust = 1.05, vjust = 1.2, size = 2, color = "black"
  ) +
  
  scale_color_manual(values = eco_color_map) +
  scale_fill_manual(values = eco_color_map) +
  theme_bw() +
  labs(x = "Shannon index", y = "Guild Functional Menhinick") +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    aspect.ratio = 1.2,
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black")
  )

p.guild_micha_binned

# 保存全局 PDF 大图
ggsave(
  "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/ok/p.guild_new.pdf",
  plot = p.guild_micha_binned, height = 5, width = 6.5, bg = "transparent"
)


###uplanded####
library(dplyr)
library(purrr)
library(ggplot2)
library(minpack.lm)
library(ggrepel) 

# ==========================================
# 1. 数据准备
# ==========================================
p.guild1 <- shannon_guild4 %>% 
  dplyr::select("Menhinick_Coef", "Shannon.TD", "Ecosystem.y") %>% 
  na.omit() 
#%>%
 # dplyr::filter(Shannon.TD > 0)

p.guild1$Ecosystem.y <- gsub("Artificial surfaces", "Artificial Surfaces", p.guild1$Ecosystem.y)

# ==========================================
# 2. 数据装箱 (Binning) 并划分两大阵营
# ==========================================
bin_width <- 0.4
binned_data <- p.guild1 %>%
  dplyr::mutate(
    bin_x = floor(Shannon.TD / bin_width) * bin_width + (bin_width / 2),
    # 🌟【核心新增】：自动划定陆地与水淹阵营
    Eco_Type = ifelse(Ecosystem.y %in% c("Paddy soil", "Wetland"), "Flooded ecosystems", "Upland ecosystems")
  ) %>%
  dplyr::group_by(Ecosystem.y, Eco_Type, bin_x) %>% # 包含 Eco_Type 进行分组
  dplyr::summarise(
    bin_y = mean(Menhinick_Coef, na.rm = TRUE),
    n_points = n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_points >= 3) 

# ==========================================
# 3. 重新拟合，并计算极值 Tipping Point
# ==========================================
eq_data_binned <- binned_data %>%
  dplyr::group_by(Ecosystem.y, Eco_Type) %>% 
  dplyr::filter(!Ecosystem.y %in% c("Shrubland", "Bare Land")) %>% 
  dplyr::group_modify(~ {
    
    b_start <- 0.5
    a_start <- max(.x$bin_y, na.rm = TRUE) * b_start * exp(1)
    
    model <- tryCatch({
      nlsLM(
        bin_y ~ a * bin_x * exp(-b * bin_x), 
        data = .x,
        start = list(a = a_start, b = b_start),
        lower = c(a = 0, b = 0.0001),        
        upper = c(a = Inf, b = Inf), 
        control = nls.lm.control(maxiter = 1000)
      )
    }, error = function(e) return(NULL))
    
    if (is.null(model)) {
      return(tibble(eq_text = "Fit failed", a=NA, b=NA, pval=NA, r2=NA, Tipping_X=NA, Tipping_Y=NA, point_label=NA))
    }
    
    coefs <- coef(model)
    a_val <- coefs["a"]
    b_val <- coefs["b"]
    
    coef_summary <- summary(model)$coefficients
    col_idx <- if(ncol(coef_summary) >= 4) 4 else ncol(coef_summary)
    pval <- tryCatch({ coef_summary["b", col_idx] }, error = function(e) NA_real_)
    p_text <- if (is.na(pval)) "p = NA" else if (pval < 0.001) "p < 0.001" else sprintf("p = %.3f", pval)
    
    rss <- sum(residuals(model)^2)
    tss <- sum((.x$bin_y - mean(.x$bin_y))^2)
    pseudo_r2 <- ifelse(tss == 0, NA, 1 - (rss / tss))
    if (!is.na(pseudo_r2) && pseudo_r2 < 0) pseudo_r2 <- 0.01 
    
    r2_text <- if (is.na(pseudo_r2)) "R² = NA" else sprintf("R² = %.2f", pseudo_r2)
    
    # 紧凑单行公式排版
    eq_text <- sprintf(
      "%s: y = %.2f*x*e^(-%.2f*x), %s, %s",
      .y$Ecosystem.y, a_val, b_val, r2_text, p_text
    )
    
    # 微积分精确解求 Ricker 波峰
    tipping_x_val <- 1 / b_val
    tipping_y_val <- (a_val / b_val) * exp(-1)
    
    point_label_text <- sprintf("%s Tipping point: (%.2f, %.2f)", .y$Ecosystem.y, tipping_x_val, tipping_y_val)
    
    tibble(
      eq_text = eq_text, a = a_val, b = b_val, pval = pval, r2 = pseudo_r2,
      Tipping_X = tipping_x_val, Tipping_Y = tipping_y_val, point_label = point_label_text
    )
  }) %>%
  ungroup()

# ==========================================
# 3.5 核心新增：将同一个阵营的公式合并为多行文本块
# ==========================================
text_labels_binned <- eq_data_binned %>%
  dplyr::group_by(Eco_Type) %>%
  dplyr::summarise(
    combined_text = paste(eq_text, collapse = "\n"),
    .groups = "drop"
  )

# ==========================================
# ==========================================
# 4. 绘图准备（全数据框 Eco_Type 因子化同步解锁）
# ==========================================
p.guild2 <- p.guild1 %>%
  dplyr::filter(!Ecosystem.y %in% c("Shrubland", "Bare Land")) %>%
  dplyr::mutate(
    Eco_Type = ifelse(Ecosystem.y %in% c("Paddy soil", "Wetland"), "Flooded ecosystems", "Upland ecosystems")
  )

# 🌟【核心绝杀】：让所有数据框的分面标签共享完全一模一样的因子顺序！
target_facet_levels <- c("Upland ecosystems", "Flooded ecosystems")

p.guild2$Eco_Type           <- factor(p.guild2$Eco_Type, levels = target_facet_levels)
binned_data$Eco_Type         <- factor(binned_data$Eco_Type, levels = target_facet_levels)
eq_data_binned$Eco_Type      <- factor(eq_data_binned$Eco_Type, levels = target_facet_levels)
text_labels_binned$Eco_Type  <- factor(text_labels_binned$Eco_Type, levels = target_facet_levels) # 连文本框也不能放过

# 保持你原本的生境级别排序
eco_levels <- c("Agricultural Land", "Grassland", "Artificial Surfaces", "Tundra", "Forest", "Wetland", "Paddy soil")
p.guild2$Ecosystem.y <- factor(p.guild2$Ecosystem.y, levels = eco_levels)
binned_data$Ecosystem.y <- factor(binned_data$Ecosystem.y, levels = eco_levels)
eq_data_binned$Ecosystem.y <- factor(eq_data_binned$Ecosystem.y, levels = eco_levels)

global_b_start <- 0.5
global_a_start <- max(binned_data$bin_y, na.rm = TRUE) * global_b_start * exp(1)

# ==========================================
# 5. 绘图展现 (双面板对峙模式)
# ==========================================
p.guild_micha_binned <- ggplot() +
  # 图层 1: 原始散点打底
 # geom_point(
  #  data = p.guild2,
   # aes(x = Shannon.TD, y = Menhinick_Coef, color = Ecosystem.y),
    #alpha = 0.12, size = 0.8
  #) +
  
  # 图层 2: 装箱均值点
  #geom_point(
   # data = binned_data %>% filter(!is.na(Ecosystem.y)),
    #aes(x = bin_x, y = bin_y, fill = Ecosystem.y),
    #shape = 21, color = "black", size = 1.8, stroke = 0.4, alpha = 0.7
  #) +
  
  # 图层 3: geom_smooth 拟合曲线 (此时同一 Panel 内会自适应绘制多条线)
  geom_smooth(
    data = binned_data %>% filter(!is.na(Ecosystem.y)),
    aes(x = bin_x, y = bin_y, color = Ecosystem.y),
    method = "nls",
    formula = y ~ a * x * exp(-b * x), 
    method.args = list(
      start = list(a = global_a_start, b = global_b_start),
      algorithm = "port",
      lower = c(a = 0, b = 0.0001),     
      upper = c(a = Inf, b = Inf), 
      control = nls.control(maxiter = 1000, warnOnly = TRUE) 
    ),
    se = FALSE,
    fullrange = TRUE,  
    linewidth = 0.8
  ) +
  
  # 图层 4：各曲线波峰垂直下落虚线
  geom_segment(
    data = eq_data_binned %>% filter(!is.na(Tipping_X)),
    aes(x = Tipping_X, xend = Tipping_X, y = 0, yend = Tipping_Y, color = Ecosystem.y),
    linetype = "dashed", alpha = 0.6, linewidth = 0.4
  ) +
  
  # 图层 5：死死钉在波峰顶部的钻石标记
  geom_point(
    data = eq_data_binned %>% filter(!is.na(Tipping_X)),
    aes(x = Tipping_X, y = Tipping_Y, fill = Ecosystem.y),
    shape = 23, color = "black", size = 2.8, stroke = 0.5, alpha = 0.9
  ) +
  
  # 图层 6：智能避让标签 (完美错开密集的生境标签与坐标)
  ggrepel::geom_text_repel(
    data = eq_data_binned %>% filter(!is.na(Tipping_X)),
    aes(x = Tipping_X, y = Tipping_Y, label = point_label, color = Ecosystem.y),
    size = 2.4, fontface = "bold", lineheight = 0.9,
    box_padding = 0.6, point_padding = 0.4,
    force = 3, max.overlaps = Inf, segment.color = "gray60", segment.linewidth = 0.3
  ) +
  
  # 🌟【关键改动】：由 7 列分面改为按阵营分面，X轴刻度独立自适应
  facet_wrap(~Eco_Type, scales = "free_x") +
  
  # 🌟【关键改动】：渲染拼接好的多行公式框
  geom_text(
    data = text_labels_binned,
    aes(x = Inf, y = Inf, label = combined_text),
    inherit.aes = FALSE,
    hjust = 1.05, vjust = 1.1, size = 2.2, color = "black", lineheight = 1.1
  ) +
  
  scale_color_manual(values = eco_color_map) +
  scale_fill_manual(values = eco_color_map) +
  theme_bw() +
  labs(
    x = "Shannon index", 
    y = "Guild Functional Menhinick"
  ) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none", # 🌟 因为多线重叠，在底部重新开启 Legend 图例以供对照
    aspect.ratio = 0.4,         # 完美的正方形多线图排版
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 9, color = "black")
  )

p.guild_micha_binned

# 保存大图
ggsave(
  "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/ok/p.guild_new_panels.pdf",
  plot = p.guild_micha_binned, height = 5, width = 10, bg = "transparent"
)


#total####

library(dplyr)
library(purrr)
library(ggplot2)
library(minpack.lm)
library(ggrepel) 

# ==========================================
# 1. 数据准备
# ==========================================
p.guild1 <- shannon_guild4 %>% 
  dplyr::select("Menhinick_Coef", "Shannon.TD") %>% # 🌟 移除了 Ecosystem.y，完全合并数据
  na.omit()

# ==========================================
# 2. 全局数据装箱 (Global Binning) 
# ==========================================
bin_width <- 0.4
binned_data <- p.guild1 %>%
  dplyr::mutate(
    bin_x = floor(Shannon.TD / bin_width) * bin_width + (bin_width / 2)
  ) %>%
  dplyr::group_by(bin_x) %>% # 🌟 只按 X 轴区间分组，计算全局均值
  dplyr::summarise(
    bin_y = mean(Menhinick_Coef, na.rm = TRUE),
    n_points = n(),
    .groups = "drop"
  ) %>%
  dplyr::filter(n_points >= 3) 

# ==========================================
# 3. 全局非线性拟合与 Tipping Point 计算
# ==========================================
b_start <- 0.5
a_start <- max(binned_data$bin_y, na.rm = TRUE) * b_start * exp(1)

# 🌟 直接对整张表进行单次非线性拟合
global_model <- tryCatch({
  nlsLM(
    bin_y ~ a * bin_x * exp(-b * bin_x), 
    data = binned_data,
    start = list(a = a_start, b = b_start),
    lower = c(a = 0, b = 0.0001),        
    upper = c(a = Inf, b = Inf), 
    control = nls.lm.control(maxiter = 1000)
  )
}, error = function(e) return(NULL))

# 提取全局拟合参数与统计量
if (!is.null(global_model)) {
  coefs <- coef(global_model)
  a_val <- coefs["a"]
  b_val <- coefs["b"]
  
  coef_summary <- summary(global_model)$coefficients
  col_idx <- if(ncol(coef_summary) >= 4) 4 else ncol(coef_summary)
  pval <- tryCatch({ coef_summary["b", col_idx] }, error = function(e) NA_real_)
  p_text <- if (is.na(pval)) "p = NA" else if (pval < 0.001) "p < 0.001" else sprintf("p = %.3f", pval)
  
  rss <- sum(residuals(global_model)^2)
  tss <- sum((binned_data$bin_y - mean(binned_data$bin_y))^2)
  pseudo_r2 <- ifelse(tss == 0, NA, 1 - (rss / tss))
  if (!is.na(pseudo_r2) && pseudo_r2 < 0) pseudo_r2 <- 0.01 
  r2_text <- if (is.na(pseudo_r2)) "R² = NA" else sprintf("R² = %.2f", pseudo_r2)
  
  eq_text <- sprintf("Global Fit: y = %.2f*x*e^(-%.2f*x)\n%s, %s", a_val, b_val, r2_text, p_text)
  
  # 计算全局最高峰坐标
  tipping_x_val <- 1 / b_val
  tipping_y_val <- (a_val / b_val) * exp(-1)
  point_label_text <- sprintf("Global Peak\n(%.2f, %.2f)", tipping_x_val, tipping_y_val)
  
  # 构建单行全局标签矩阵
  eq_data_binned <- tibble(
    eq_text = eq_text, a = a_val, b = b_val, pval = pval, r2 = pseudo_r2,
    Tipping_X = tipping_x_val, Tipping_Y = tipping_y_val, point_label = point_label_text
  )
} else {
  eq_data_binned <- tibble(eq_text = "Fit failed", a=NA, b=NA, pval=NA, r2=NA, Tipping_X=NA, Tipping_Y=NA, point_label=NA)
}

# ==========================================
# 4. 绘图准备 (采用高端清爽的经典学术配色)
# ==========================================
global_b_start <- 0.5
global_a_start <- max(binned_data$bin_y, na.rm = TRUE) * global_b_start * exp(1)

# ==========================================
# 5. 绘图展现 (单面板高效呈现)
# ==========================================
p.global_ricker <- ggplot() +
  # 图层 1: 原始全局散点（使用优雅的浅灰色打底，突显主体）
  geom_point(
    data = p.guild1,
    aes(x = Shannon.TD, y = Menhinick_Coef),
    color = "gray80", alpha = 0.2, size = 0.8
  ) +
  
  # 图层 2: 全局装箱均值点（深碳色圆点）
  geom_point(
    data = binned_data,
    aes(x = bin_x, y = bin_y),
    shape = 21, fill = "gray30", color = "black", size = 1.8, stroke = 0.5, alpha = 0.8
  ) +
  
  # 图层 3: 全局唯一的 geom_smooth Ricker 拟合曲线（学术深蓝色）
  geom_smooth(
    data = binned_data,
    aes(x = bin_x, y = bin_y),
    method = "nls",
    formula = y ~ a * x * exp(-b * x), 
    method.args = list(
      start = list(a = global_a_start, b = global_b_start),
      algorithm = "port",
      lower = c(a = 0, b = 0.0001),     
      upper = c(a = Inf, b = Inf), 
      control = nls.control(maxiter = 1000, warnOnly = TRUE) 
    ),
    se = FALSE,
    fullrange = TRUE,  
    color = "#1f78b4", linewidth = 1.0
  ) +
  
  # 图层 4: 全局波峰垂直下落虚线
  geom_segment(
    data = eq_data_binned %>% filter(!is.na(Tipping_X)),
    aes(x = Tipping_X, xend = Tipping_X, y = 0, yend = Tipping_Y),
    color = "firebrick", linetype = "dashed", alpha = 0.7, linewidth = 0.4
  ) +
  
  # 图层 5: 钉在全局曲线最高峰顶部的红色钻石标记
  geom_point(
    data = eq_data_binned %>% filter(!is.na(Tipping_X)),
    aes(x = Tipping_X, y = Tipping_Y),
    shape = 23, fill = "firebrick", color = "black", size = 3.0, stroke = 0.5
  ) +
  
  # 图层 6: 全局波峰坐标标注
  ggrepel::geom_text_repel(
    data = eq_data_binned %>% filter(!is.na(Tipping_X)),
    aes(x = Tipping_X, y = Tipping_Y, label = point_label),
    color = "firebrick", size = 2.6, fontface = "bold", lineheight = 0.9,
    box_padding = 0.5, point_padding = 0.4, force = 2
  ) +
  
  # 右上角全局公式文本框
  geom_text(
    data = eq_data_binned,
    aes(x = Inf, y = Inf, label = eq_text),
    inherit.aes = FALSE,
    hjust = 1.05, vjust = 1.3, size = 2.6, color = "black", lineheight = 1.1
  ) +
  
  theme_bw() +
  labs(
    x = "Taxonomic Diversity (Shannon index)", 
    y = "Guild Functional Divergence (Menhinick index)"
  ) +
  theme(
    panel.grid = element_blank(),
    aspect.ratio = 0.6, # 🌟 顶级期刊标准正方形排版
    axis.text = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 10, face = "bold")
  )

p.global_ricker

# 保存全局大图
ggsave(
  "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/ok/p.guild_global.pdf",
  plot = p.global_ricker, height = 4, width = 4.3, bg = "transparent"
)

#projection####

all_data_abs<-readRDS("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/all.data.abs.rds")

library(dplyr)
library(purrr)
library(tidyr)
library(minpack.lm)
library(ggplot2)

# ==========================================
# 1. 高效数据装箱 (Joint Binning)
# ==========================================
bin_width <- 0.4

binned_future <- all_data_abs %>%
  # 滤除无效数据
  dplyr::filter(TD > 0, !is.na(FD), !is.na(TD), !is.na(Habitat), !is.na(Scenario)) %>%
  dplyr::mutate(
    bin_x = floor(TD / bin_width) * bin_width + (bin_width / 2)
  ) %>%
  # 联合分组：按情景、生境、X区间
  dplyr::group_by(Scenario, Habitat, bin_x) %>%
  dplyr::summarise(
    bin_y = mean(FD, na.rm = TRUE),
    n_points = n(),
    .groups = "drop"
  ) %>%
  # 保持统计稳健性，滤除孤立点
  dplyr::filter(n_points >= 3)

# ==========================================
# 2. 批量拟合模型并提取 Tipping Point
# ==========================================
tipping_results <- binned_future %>%
  dplyr::group_by(Scenario, Habitat) %>%
  # 确保该组内至少有4个有效装箱点，否则 nlsLM 会因自由度不足报错
  dplyr::filter(n() >= 4) %>%
  dplyr::group_modify(~ {
    
    b_start <- 0.5
    a_start <- max(.x$bin_y, na.rm = TRUE) * b_start * exp(1)
    
    model <- tryCatch({
      nlsLM(
        bin_y ~ a * bin_x * exp(-b * bin_x), 
        data = .x,
        start = list(a = a_start, b = b_start),
        lower = c(a = 0, b = 0.0001),        
        upper = c(a = Inf, b = Inf), 
        control = nls.lm.control(maxiter = 1000)
      )
    }, error = function(e) return(NULL))
    
    # 捕获拟合失败的组
    if (is.null(model)) {
      return(tibble(a = NA, b = NA, Tipping_X = NA, Tipping_Y = NA, r2 = NA))
    }
    
    coefs <- coef(model)
    a_val <- coefs["a"]
    b_val <- coefs["b"]
    
    # 计算拟合优度
    rss <- sum(residuals(model)^2)
    tss <- sum((.x$bin_y - mean(.x$bin_y))^2)
    pseudo_r2 <- ifelse(tss == 0, NA, 1 - (rss / tss))
    if (!is.na(pseudo_r2) && pseudo_r2 < 0) pseudo_r2 <- 0.01 
    
    # 微积分解析解提取生态转折点
    tipping_x_val <- 1 / b_val
    tipping_y_val <- (a_val / b_val) * exp(-1)
    
    tibble(
      a = a_val, 
      b = b_val, 
      Tipping_X = tipping_x_val, 
      Tipping_Y = tipping_y_val, 
      r2 = pseudo_r2
    )
  }) %>%
  ungroup()

# 查看计算结果
print(tipping_results)



#different level relationship####

library(tidyverse)
library(vegan)


# abundance矩阵
abundance1645 <- read.csv(
  '/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/paper/写作参考/nature/submit/Supplementary Table/Supplementary Table 6.csv'
)
abundance1645<-abundance1645[,-1201]

# 查看列名
colnames(abundance1645)


# Genome作为行名
abund <- abundance1645 %>%
  column_to_rownames("Genome")


# 如果第一列不是sample丰度，删除
# 根据你的表格修改
abund <- abund[,2:ncol(abund)]


# 转numeric
abund <- as.data.frame(
  lapply(
    abund,
    as.numeric
  ),
  row.names = rownames(abund)
)


# 转matrix
abund <- as.matrix(abund)


# 检查
dim(abund)
head(abund)


make_taxon_id <- function(taxonomy, level){
  
  tax_levels <- c(
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
  )
  
  
  use_levels <- tax_levels[
    1:match(level, tax_levels)
  ]
  
  
  print(level)
  print(use_levels)
  
  
  taxonomy$taxon_id <- apply(
    taxonomy[, use_levels, drop=FALSE],
    1,
    function(x){
      
      x[is.na(x)] <- "unknown"
      
      paste(
        x,
        collapse="|"
      )
      
    }
  )
  
  
  return(taxonomy)
}

calculate_shannon_taxa <- function(abund,
                                   taxonomy,
                                   level){
  
  
  # 创建唯一taxon ID
  
  taxonomy <- make_taxon_id(
    taxonomy,
    level
  )
  
  
  # genome-taxonomy对应关系
  
  tax <- taxonomy %>%
    select(
      genome,
      taxon_id
    )
  
  
  # abundance转换
  
  abund_df <- abund %>%
    as.data.frame() %>%
    rownames_to_column("genome")
  
  
  # 匹配taxonomy
  
  merge_df <- merge(
    abund_df,
    tax,
    by="genome"
  )
  
  
  # 聚合到分类水平
  agg <- merge_df %>%
    group_by(
      taxon_id
    ) %>%
    summarise(
      across(
        where(is.numeric),
        sum
      )
    )
  
  # 转matrix
  
  mat <- agg %>%
    column_to_rownames(
      "taxon_id"
    ) %>%
    as.matrix()
  
  
  # 强制numeric
  
  storage.mode(mat) <- "numeric"
  
  
  # Shannon
  
  shannon <- diversity(
    t(mat),
    index="shannon"
  )
  
  
  result <- data.frame(
    Sample=colnames(mat),
    Shannon=shannon,
    Level=level
  )
  
  
  return(result)
  
}

levels <- c(
  "phylum",
  "class",
  "order",
  "family",
  "genus",
  "species"
)


shannon_all <- lapply(
  
  levels,
  
  function(level){
    
    calculate_shannon_taxa(
      abund = abund,
      taxonomy = arc1645,
      level = level
    )
    
  }
  
)


shannon_all <- bind_rows(
  shannon_all
)


head(shannon_all)


library(dplyr)


gfm_data <- shannon_guild4 %>%
  select(
    Sample,
    Menhinick_Coef,
    Ecosystem.y
  )

shannon_gfm <- shannon_all %>%
  left_join(
    gfm_data,
    by="Sample"
  ) %>%
  na.omit()

p.guild1 <- shannon_gfm %>%
  select(
    Shannon,
    Menhinick_Coef,
    Level,
    Ecosystem.y
  ) %>%
  na.omit()



# ==========================================
# 1. 数据准备
# ==========================================

# dplyr::filter(Shannon.TD > 0)

p.guild1$Ecosystem.y <- gsub("Artificial surfaces", "Artificial Surfaces", p.guild1$Ecosystem.y)

eco_levels <- c("Agricultural Land", "Grassland", "Artificial Surfaces","Tundra", "Forest", "Wetland", "Paddy soil")
p.guild1$Ecosystem.y <- factor(p.guild1$Ecosystem.y, levels = eco_levels)
# ==========================================
# 2. 数据装箱 (Binning) 并划分两大阵营
# ==========================================
bin_width <- 0.4
binned_data <- p.guild1 %>%
  mutate(
    bin_x =
      floor(Shannon / bin_width) *
      bin_width +
      (bin_width/2),
    # 🌟【核心新增】：自动划定陆地与水淹阵营
    Eco_Type = ifelse(Ecosystem.y %in% c("Paddy soil", "Wetland"), "Flooded ecosystems", "Upland ecosystems")
  ) %>%
  group_by(
    Level,
    Ecosystem.y,
    Eco_Type,
    bin_x
  ) %>%
  summarise(
    bin_y=mean(Menhinick_Coef),
    n_points=n(),
    .groups="drop"
  )

# ==========================================
# 3. 重新拟合，并计算极值 Tipping Point
# ==========================================
eq_data_binned <- binned_data %>%
  dplyr::group_by(Ecosystem.y, Eco_Type,Level) %>% 
  dplyr::filter(!Ecosystem.y %in% c("Shrubland", "Bare Land")) %>% 
  dplyr::group_modify(~ {
    
    b_start <- 0.5
    a_start <- max(.x$bin_y, na.rm = TRUE) * b_start * exp(1)
    
    model <- tryCatch({
      nlsLM(
        bin_y ~ a * bin_x * exp(-b * bin_x), 
        data = .x,
        start = list(a = a_start, b = b_start),
        lower = c(a = 0, b = 0.0001),        
        upper = c(a = Inf, b = Inf), 
        control = nls.lm.control(maxiter = 1000)
      )
    }, error = function(e) return(NULL))
    
    if (is.null(model)) {
      return(tibble(eq_text = "Fit failed", a=NA, b=NA, pval=NA, r2=NA, Tipping_X=NA, Tipping_Y=NA, point_label=NA))
    }
    
    coefs <- coef(model)
    a_val <- coefs["a"]
    b_val <- coefs["b"]
    
    coef_summary <- summary(model)$coefficients
    col_idx <- if(ncol(coef_summary) >= 4) 4 else ncol(coef_summary)
    pval <- tryCatch({ coef_summary["b", col_idx] }, error = function(e) NA_real_)
    p_text <- if (is.na(pval)) "p = NA" else if (pval < 0.001) "p < 0.001" else sprintf("p = %.3f", pval)
    
    rss <- sum(residuals(model)^2)
    tss <- sum((.x$bin_y - mean(.x$bin_y))^2)
    pseudo_r2 <- ifelse(tss == 0, NA, 1 - (rss / tss))
    if (!is.na(pseudo_r2) && pseudo_r2 < 0) pseudo_r2 <- 0.01 
    
    r2_text <- if (is.na(pseudo_r2)) "R² = NA" else sprintf("R² = %.2f", pseudo_r2)
    
    # 紧凑单行公式排版
    eq_text <- sprintf(
      "%s: y = %.2f*x*e^(-%.2f*x), %s, %s",
      .y$Ecosystem.y, a_val, b_val, r2_text, p_text
    )
    
    # 微积分精确解求 Ricker 波峰
    tipping_x_val <- 1 / b_val
    tipping_y_val <- (a_val / b_val) * exp(-1)
    
    point_label_text <- sprintf("%s Tipping point: (%.2f, %.2f)", .y$Ecosystem.y, tipping_x_val, tipping_y_val)
    
    tibble(
      eq_text = eq_text, a = a_val, b = b_val, pval = pval, r2 = pseudo_r2,
      Tipping_X = tipping_x_val, Tipping_Y = tipping_y_val, point_label = point_label_text
    )
  }) %>%
  ungroup()

# ==========================================
# 3.5 核心新增：将同一个阵营的公式合并为多行文本块
# ==========================================
text_labels_binned <- eq_data_binned %>%
  dplyr::group_by(Eco_Type,Level) %>%
  dplyr::summarise(
    combined_text = paste(eq_text, collapse = "\n"),
    .groups = "drop"
  )

# ==========================================
# ==========================================
# 4. 绘图准备（全数据框 Eco_Type 因子化同步解锁）
# ==========================================
p.guild2 <- p.guild1 %>%
  dplyr::filter(!Ecosystem.y %in% c("Shrubland", "Bare Land")) %>%
  dplyr::mutate(
    Eco_Type = ifelse(Ecosystem.y %in% c("Paddy soil", "Wetland"), "Flooded ecosystems", "Upland ecosystems")
  )

# 🌟【核心绝杀】：让所有数据框的分面标签共享完全一模一样的因子顺序！
target_facet_levels <- c("Upland ecosystems", "Flooded ecosystems")

p.guild2$Eco_Type           <- factor(p.guild2$Eco_Type, levels = target_facet_levels)
binned_data$Eco_Type         <- factor(binned_data$Eco_Type, levels = target_facet_levels)
eq_data_binned$Eco_Type      <- factor(eq_data_binned$Eco_Type, levels = target_facet_levels)
text_labels_binned$Eco_Type  <- factor(text_labels_binned$Eco_Type, levels = target_facet_levels) # 连文本框也不能放过

# 保持你原本的生境级别排序
eco_levels <- c("Agricultural Land", "Grassland", "Artificial Surfaces", "Tundra", "Forest", "Wetland", "Paddy soil")
p.guild2$Ecosystem.y <- factor(p.guild2$Ecosystem.y, levels = eco_levels)
binned_data$Ecosystem.y <- factor(binned_data$Ecosystem.y, levels = eco_levels)
eq_data_binned$Ecosystem.y <- factor(eq_data_binned$Ecosystem.y, levels = eco_levels)

global_b_start <- 0.5
global_a_start <- max(binned_data$bin_y, na.rm = TRUE) * global_b_start * exp(1)

# ==========================================
# 5. 绘图展现 (双面板对峙模式)
# ==========================================
p.guild_micha_binned <- ggplot() +
  # 图层 1: 原始散点打底
  # geom_point(
  #  data = p.guild2,
  # aes(x = Shannon.TD, y = Menhinick_Coef, color = Ecosystem.y),
  #alpha = 0.12, size = 0.8
  #) +
  
  # 图层 2: 装箱均值点
  #geom_point(
  # data = binned_data %>% filter(!is.na(Ecosystem.y)),
  #aes(x = bin_x, y = bin_y, fill = Ecosystem.y),
#shape = 21, color = "black", size = 1.8, stroke = 0.4, alpha = 0.7
#) +

# 图层 3: geom_smooth 拟合曲线 (此时同一 Panel 内会自适应绘制多条线)
geom_smooth(
  data = na.omit(binned_data) %>% filter(!is.na(Ecosystem.y)),
  aes(x = bin_x, y = bin_y, color = Ecosystem.y),
  method = "nls",
  formula = y ~ a * x * exp(-b * x), 
  method.args = list(
    start = list(a = global_a_start, b = global_b_start),
    algorithm = "port",
    lower = c(a = 0, b = 0.0001),     
   # upper = c(a = Inf, b = Inf), 
    control = nls.control(maxiter = 1000, warnOnly = TRUE) 
  ),
  se = FALSE,
  fullrange = TRUE,  
  linewidth = 0.8
) +
  
  # 图层 4：各曲线波峰垂直下落虚线
  geom_segment(
    data = eq_data_binned %>% filter(!is.na(Tipping_X)),
    aes(x = Tipping_X, xend = Tipping_X, y = 0, yend = Tipping_Y, color = Ecosystem.y),
    linetype = "dashed", alpha = 0.6, linewidth = 0.4
  ) +
  
  # 图层 5：死死钉在波峰顶部的钻石标记
  geom_point(
    data = eq_data_binned %>% filter(!is.na(Tipping_X)),
    aes(x = Tipping_X, y = Tipping_Y, fill = Ecosystem.y),
    shape = 23, color = "black", size = 2.8, stroke = 0.5, alpha = 0.9
  ) +
  
  # 图层 6：智能避让标签 (完美错开密集的生境标签与坐标)
  ggrepel::geom_text_repel(
    data = eq_data_binned %>% filter(!is.na(Tipping_X)),
    aes(x = Tipping_X, y = Tipping_Y, label = point_label, color = Ecosystem.y),
    size = 2.4, fontface = "bold", lineheight = 0.9,
    box_padding = 0.6, point_padding = 0.4,
    force = 3, max.overlaps = Inf, segment.color = "gray60", segment.linewidth = 0.3
  ) +
  
  # 🌟【关键改动】：由 7 列分面改为按阵营分面，X轴刻度独立自适应
  facet_wrap(Level~Eco_Type, scales = "free_x") +
  
  # 🌟【关键改动】：渲染拼接好的多行公式框
  geom_text(
    data = text_labels_binned,
    aes(x = Inf, y = Inf, label = combined_text),
    inherit.aes = FALSE,
    hjust = 1.05, vjust = 1.1, size = 2.2, color = "black", lineheight = 1.1
  ) +
  
  scale_color_manual(values = eco_color_map) +
  scale_fill_manual(values = eco_color_map) +
  theme_bw() +
  labs(
    x = "Shannon index", 
    y = "Guild Functional Menhinick"
  ) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none", # 🌟 因为多线重叠，在底部重新开启 Legend 图例以供对照
    aspect.ratio = 0.4,         # 完美的正方形多线图排版
    strip.background = element_rect(fill = "gray95"),
    strip.text = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 9, color = "black")
  )

p.guild_micha_binned
