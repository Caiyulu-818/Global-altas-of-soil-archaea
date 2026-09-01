#双violin####
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggpubr) # 必须加载 ggpubr
#library(plyr)   # 必须加载 plyr

# ==========================================
# 【第一部分】：定义分裂小提琴函数
# ==========================================
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- grid::polylineGrob(both$x, both$y, id = rep(1:(nrow(quantiles)/2), each = 2),
                                                                   default.units = "native",
                                                                   gp = grid::gpar(col = both$colour, lwd = both$linewidth, lty = both$linetype, alpha = both$alpha))
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             } else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           }
)

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, 
        show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

# ==========================================
# 【第二部分】：数据清洗、排序与 Type 分组
# ==========================================
scale_factor <- 3
sampleall.arc<-read.csv("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/submit/data/supp/Supplementary Table 1.csv")
paddy<-read.csv("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/sampleinfo/sample.paddy.soil_filter1.csv")
a.shannonguild3<-read.csv("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/shannonguild3.csv")
a.shannonguild4<-merge(a.shannonguild3,paddy,by.x="Sample.y",by.y="Sample")
shannon_guild4=a.shannonguild4


write.csv(shannon_guild4,"/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/shannon_guild4.csv")
shannon_guild4<-read.csv("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/shannon_guild4.csv")
plot_data <- shannon_guild4 %>%
  dplyr::select(Ecosystem.y, Shannon.TD, Menhinick_Coef) %>%
  na.omit() %>%
  filter(!Ecosystem.y%in%c("Bare Land","Shrubland")) %>% 
  dplyr::mutate(Menhinick_Coef_scaled = Menhinick_Coef * scale_factor) %>%
  tidyr::pivot_longer(
    cols = c(Shannon.TD, Menhinick_Coef_scaled),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  dplyr::mutate(Metric = dplyr::recode(Metric,
                                       "Shannon.TD" = "Taxonomic (Shannon)",
                                       "Menhinick_Coef_scaled" = "Functional (Menhinick)"))

# ✨ 添加 Type 分组 (优化写法，更整洁)
plot_data <- plot_data %>% 
  dplyr::mutate(type = ifelse(Ecosystem.y%in% c("Paddy soil","Wetland"), "Flooded Ecosystem.y", "Upland Ecosystem.y"))

plot_data$Metric <- factor(plot_data$Metric, levels = c("Taxonomic (Shannon)", "Functional (Menhinick)"))

# 严谨的分组计算与降序排序
plot_data <- plot_data %>%
  dplyr::group_by(Ecosystem.y) %>%
  dplyr::mutate(
    shannon_median = median(Value[Metric == "Taxonomic (Shannon)"], na.rm = TRUE),
    menhinick_median = median(Value[Metric == "Functional (Menhinick)"], na.rm = TRUE),
    diff_gap = shannon_median - menhinick_median,
    diy=diff_gap/2+ menhinick_median
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Ecosystem.y = fct_reorder(Ecosystem.y, diff_gap, .desc = TRUE))
x<-unique(plot_data[,c(1,8)])
x
# ==========================================
# 【第三部分】：提取坐标，计算“阶梯括号”中点与分面 X 坐标
# ==========================================
dodge_offset <- 0.05

mean_data <- plot_data %>%
  dplyr::group_by(Ecosystem.y, Metric, type) %>%  # 必须带上 type 以便后续作图
  dplyr::summarise(mean_val = median(Value, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Metric, values_from = mean_val) %>%
  dplyr::rename(
    Shannon_mean = `Taxonomic (Shannon)`,
    Menhinick_mean_neg = `Functional (Menhinick)`
  ) %>%
 
  
  # 🚨 核心魔法：按 type 分组计算分面后的相对 X 坐标
  dplyr::group_by(type) %>%
  dplyr::mutate(
    x_base = as.integer(droplevels(Ecosystem.y)), # 重新生成分面内的 1, 2, 3 坐标
    x_left = x_base - dodge_offset,  
    x_right = x_base + dodge_offset  
  ) %>%
  dplyr::ungroup()
mean_data
# ==========================================
# 【第四部分】：绘制带“中点连线”的终极分面图
# ==========================================
# 假设你的颜色向量，这里给一个示范，请换回你的 metric_colors
metric_colors <- c("Taxonomic (Shannon)" = "#bb7a8c", "Functional (Menhinick)" = "#509296")

plot_data %>% filter(type=="Flooded Ecosystem.y") %>% distinct() %>% 
  ggplot( aes(x = Ecosystem.y, y = Value, 
                      fill = Metric, color = Metric, 
                      group = interaction(Metric, Ecosystem.y))) +
  
  # 🌟 图层 1：分裂小提琴
  geom_split_violin(
    aes(alpha = Metric), 
    width = 0.9, trim = TRUE, scale = "width"
  ) +  
  
  # 🌟 图层 2：添加误差棒与均值点
  stat_summary(aes(color = Metric),
               fun.data = "mean_sd", geom = "errorbar", 
               width = 0.1, position = position_dodge(width = 0.3), 
               alpha = 0.8, linewidth = 0.6) +
  stat_summary(aes(color = Metric),
               fun = "mean", geom = "point", 
               position = position_dodge(width = 0.25), size = 1.5) +


# ==========================================
# 【第四部分】：绘制带箱线图的终极分面图
# ==========================================
metric_colors <- c("Taxonomic (Shannon)" = "#bb7a8c", "Functional (Menhinick)" = "#509296")



#qiayao####
# 在绘图前确保您的 type 已经 factor 化
plot_data$type <- factor(plot_data$type, levels = rev(c("Flooded Ecosystem.y", "Upland Ecosystem.y")))


# 确保已加载 gghalves 包
library(gghalves)
library(ggeasy)

# 在绘图前确保您的 type 已经 factor 化
plot_data$type <- factor(plot_data$type, levels = rev(c("Flooded Ecosystem.y", "Upland Ecosystem.y")))

ggplot(plot_data, aes(x = Ecosystem.y, y = Value, fill = Metric, color = Metric)) +
  
  # 🌟 核心修改 1：添加右侧半小提琴图 (精准复刻参考图的右侧密度分布)
  # geom_half_violin(
  #   side = "r",                             # 强制分布图在右侧
  #  position = position_dodge(width = 1), # 与箱图保持相同的错开宽度
  # alpha = 0.6,                            # 设置较高的透明度，作为背景衬托
  #color = NA,                             # 去除小提琴边框，视觉更干净
  #trim = TRUE
  #) +
  
  # 🌟 核心修改 2：保留箱线图的上下误差棒 (T型横线)
#stat_boxplot(
# geom = "errorbar", 
#width = 0.15,                           # 稍微调窄，配合稍后细长的箱体
#linewidth = 0.5,                        
#position = position_dodge(width = 0.8)  
#) +

# 🌟 3. 核心魔法：双三角形沙漏箱图
geom_boxplot(
  width = 0.4,                           
  linewidth = 0.3,                        
  notch = TRUE,                           # 开启凹槽
  notchwidth = 0,                         # ⚠️ 核心魔法参数：中位数处宽度设为0，形成两个完美的三角形
  outlier.shape = NA,                     # 彻底隐藏所有散点
  # fill = "white",                         # 填白凸显线条轮廓
  alpha = 0.6,                            
  position = position_dodge(width = 0.8)  
) +
  # ==========================================
# 双 Y 轴与样式设置 (保留您的完整逻辑)
# ==========================================
scale_y_continuous(
  name = "Shannon index",
  sec.axis = sec_axis(~ . / scale_factor, name = "Guild Functional Menhinick")
) +
  facet_wrap(~ type, scales = "free_x", space = "free_x") +
  
  scale_color_manual(values = metric_colors) +
  scale_fill_manual(values = metric_colors) +
  
  theme_bw() +
  xlab(NULL) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",                
    legend.title = element_blank(),
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 30, size = 12, color = "black", hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.y.left = element_text(color = "black", size = 15, margin = ggplot2::margin(r = 10)),
    axis.title.y.right = element_text(color = "black", size = 15, margin = ggplot2::margin(l = 10), angle = -90)
  ) + 
  easy_all_text_size(size=20)

# 导出为 PDF (保持了你原有的路径和背景透明设置)
ggsave("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/f3b_boxplotall_zhui.pdf", height = 3, width = 6, bg = "transparent")


#线性拟合

library(dplyr)
library(purrr)
library(ggplot2)

# 1. 数据准备
p.guild1 <- shannon_guild4 %>% 
  dplyr::select("Menhinick_Coef", "Shannon.TD", "Ecosystem.y")
p.guild1$Ecosystem.y<-gsub("Artificial surfaces","Artificial Surfaces",p.guild1$Ecosystem.y)

# ==========================================
# 2. 批量计算每个 Ecosystem.y 的一阶线性回归公式与 P 值
# ==========================================
eq_data <- p.guild1%>%
  group_by(Ecosystem.y) %>%
  dplyr::filter(!Ecosystem.y %in% c("Shrubland", "Bare Land")) %>%
  dplyr::mutate(type = ifelse(Ecosystem.y %in% c("Paddy soil", "Wetland"), "Type 2", "Type 1")) %>%
 # filter(type=="Type 1") %>% 
  group_modify(~ {
    # ✨ 核心修正 1：拟合标准的一阶简单线性模型 (y ~ x)
    model <- lm(Menhinick_Coef ~ Shannon.TD, data = .x)
    
    # 提取系数
    coefs <- coef(model)
    
    # 提取 P 值
    f_stat <- summary(model)$fstatistic
    pval <- if (is.null(f_stat)) NA_real_ else pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)
    
    # ✨ 核心修正 2：只提取截距 (b0) 和 1次项斜率 (b1)，彻底删除 b2
    b0 <- coefs[1]  # 截距
    b1 <- coefs[2]  # x 的系数
    
    # 格式化 P 值
    p_text <- ifelse(pval < 0.001, "p < 0.001", sprintf("p = %.3f", pval))
    
    # ✨ 核心修正 3：拼接一阶线性公式文本 (y = b0 + b1x)
    eq_text <- sprintf("y = %.2f %s %.2fx\n%s",
                       b0,
                       ifelse(b1 >= 0, "+", "-"), abs(b1),
                       p_text)
    
    # 返回带有公式文本的数据框
    tibble(
      Ecosystem.y = unique(.x$Ecosystem.y),
      eq_text = eq_text,
      b0 = b0,
      b1 = b1,
      pval = pval
    )
  }) %>%
  ungroup()

# ==========================================
# 3. 绘图并映射各个分面的公式

p.guild2<-p.guild1 %>%
  dplyr::filter(!Ecosystem.y %in% c("Shrubland", "Bare Land")) %>%
  dplyr::mutate(type = ifelse(Ecosystem.y %in% c("Paddy soil", "Wetland"), "Type 2", "Type 1")) 
p.guild2$Ecosystem.y<-factor(p.guild2$Ecosystem.y,levels = c("Tundra","Forest","Grassland","Agricultural Land","Artificial Surfaces","Wetland","Paddy soil"))
eq_data$Ecosystem.y  <-factor(eq_data$Ecosystem.y ,levels = c("Tundra","Forest","Grassland","Agricultural Land","Artificial Surfaces","Wetland","Paddy soil"))
# ✨ 核心修复：按 type 的字母表顺序 (Type 1 -> Type 2) 物理排序数据行
 # dplyr::arrange(type, Ecosystem.y) %>%filter(type=="Type 1") %>% 
  ggplot(p.guild2,aes(x = Shannon.TD, y = Menhinick_Coef, fill = Ecosystem.y, color = Ecosystem.y)) +
  geom_point(alpha = 0.2, size = 1.5) + 
  
  # ✨ 核心修正 4：平滑线图层回归最纯粹的简单线性公式 y ~ x
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, alpha = 0.6, linewidth = 0.6) + 
    facet_wrap(~Ecosystem.y, ncol = 7,scales = "free_x") + 
  # 使用 geom_text 映射计算好的线性公式
  geom_text(
    data = eq_data, 
    aes(x = Inf, y = 1, label = eq_text), 
    inherit.aes = FALSE,           # 不继承主图的颜色映射，保持公式为黑色
    hjust = 1.05, vjust = -0.2,    # 定位在右下角
    size = 2, color = "black"      
  ) + scale_color_manual(values = eco_color_map) +
    scale_fill_manual(values = eco_color_map)+
  
  
 # facet_wrap(~Ecosystem.y, ncol = 7) + 
  
  theme_bw() +
  xlab("Shannon index") + 
  ylab("Guild Functional Menhinick") + 
  theme(
    panel.grid = element_blank(),
    legend.position = "none",  
    aspect.ratio = 1.2,
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black")
  ) + 
  ggeasy::easy_all_text_size(size = 10) # + 
# 如果需要使用自定义颜色，把下面两行取消注释

ggsave("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/ok/f3c_new.pdf",height = 4,width = 8,bg="transparent")



