library(tidyverse)
library(ggforce)
library(patchwork)

data <- data.frame(
  category = c("Metabolism", "Lipid metabolism", "AMPK signaling pathway", 
               "PPAR signaling pathway", "Lipid biosynthesis proteins"),
  value = c(37, 27, 18, 16, 12)
)

petals = 5
petal_angle = 360/petals

data <- data |>
  mutate(petal = row_number(),
         theta0 = petal * petal_angle) |>
  reframe(theta = theta0 + c(0, -petal_angle/2,  0, 
                             petal_angle/2, 0),
          r = value^0.5 * c(0, 0.6, 1, 0.6, 0), 
          .by = c(category, value, petal, theta0))
label_data <- data %>%
  group_by(category) %>%
  slice_max(r, n = 1) %>%
  ungroup()

p1 <- ggplot(data, aes(theta, r, group = petal, fill = category)) +
  ggforce::stat_bspline(geom = "area", n = 1000) +
  geom_text(data = label_data, aes(label = value)) +
  labs(fill = NULL) +
  coord_radial() +
  theme_void() +
  theme(legend.position = "none")
p1
p2 <- ggplot(data, aes(theta, r, group = petal, fill = category)) +
  ggforce::stat_bspline(geom = "area", n = 1000) +
  geom_text(data = label_data, aes(label = value)) +
  labs(fill = NULL) +
  ggsci::scale_fill_npg() + #指定Nature配色
  coord_radial() +
  theme_void() +
  theme(legend.position = "none")
p1 + p2
ggsave("Petal plot2.pdf", width = 6, height = 3)


env.driver<-read.csv("/Users/caiyulu/Desktop/MAGcode/sediment/mag4.0/new_MAG/arc/01.arcinfo/F3/FD/guild.driver/combined_explained_variation_taxall.csv")
env.driver<-read.csv("/Users/caiyulu/Desktop/MAGcode/sediment/mag4.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/TD_fdcombined_explained_variation_taxall.csv")

env.driver$Variable<-gsub("PCNM.*","Space",env.driver$Variable)

env.driver1 <- env.driver %>%
  dplyr::group_by(guild, Variable) %>% 
  # ✨ 核心修复：强制使用 dplyr::summarise，并带上 .groups = "drop" 自动解绑
  dplyr::summarise(I.perc = sum(I.perc, na.rm = TRUE), .groups = "drop") %>% 
  # 2. ✨ 核心补全步骤：
  # 自动提取 Variable 的所有唯一类型，并为每个 guild 补齐缺项
  tidyr::complete(
    guild, 
    Variable = unique(env.driver$Variable), 
    fill = list(I.perc = 0)  # 缺失的组合填充为 0
  ) %>% 
  # 推荐使用这种最经典稳妥的 left_join 写法
  dplyr::left_join(keep0130, by = c("Variable" = "keep"))
keep0130<-read_xlsx("/Users/caiyulu/Desktop/MAGcode/sediment/mag4.0/new_MAG/arc/01.arcinfo/F3/FD/guild.driver/keep0130.xlsx") %>% 
  filter(abbreviation!="Population density")
keep0130$cluster[is.na(keep0130$cluster)]<-0
#keep0130 <-keep0130 %>% mutate(keep="Space") %>%mutate(abbreviation="Space",type="Geo-distance")
template <- env.driver1 %>%
  dplyr::select(abbreviation, type) %>%
  dplyr::distinct() %>%
  dplyr::arrange(type, abbreviation) %>% 
  dplyr::mutate(
    petal = row_number(),                 # 1 到 28 的绝对固定卡位
    petals_total = n(),                   # 总数 (28)
    petal_angle = 360 / petals_total,     # 每个卡位占用的角度 (360/28)
    theta0 = petal * petal_angle,       # 中心轴线角度
    I.perc=100 )

env.driver2<-env.driver1 %>% dplyr::select(guild,abbreviation,I.perc,type)
#env.driver2<-merge(env.driver2,keep0130[,c(4,5)],by ="abbreviation",all.y=TRUE)

petals = 22
petal_angle = 360/petals
env.driver1$abbreviation<-gsub("NDVImean","NDVI",env.driver1$abbreviation)
env.driver1$abbreviation<-gsub("Nppmean","NPP",env.driver1$abbreviation)

g1rose<- env.driver1 %>% filter(guild=="Guild 1") 
g2rose<- env.driver1 %>%filter(guild=="Guild 2") 
g3rose<- env.driver1 %>%filter(guild=="Guild 3") 
g4rose<- env.driver1 %>%filter(guild=="Guild 4") 
g5rose<- env.driver1 %>% filter(guild=="Guild 5") 
# 按照你要求的逆时针逻辑定义顺序####
# 顺序：植被 -> 地形 -> 土壤 -> 空间 -> 气候
fixed_order <- c(
  "NDVI", "NPP", "Evenness",                # 植被/生产力 (橘粉)
  "Northness", "Slope", "Elevation", "Eastness",    # 地形 (深蓝)
  "cfv", "Soil salinity ece", "Soil biomass",       # 土壤物理 (绿色)
  "ocs", "SOC", "TN", "pH", "ocd",                  # 土壤化学 (绿色)
  "Sand content", "Clay", "cec", "bdod",            # 土壤质地 (绿色)
  "Space",                                          # 空间 (浅蓝)
  "MAT", "MAP"                                      # 气候 (红色)
)

# 计算每个花瓣应占的角度总分（根据变量总数平分 360 度）
total_vars <- length(fixed_order)+1
petal_angle <- 2 * pi / total_vars 

g4rose.plot <- g4rose %>%
  # 1. 补全所有缺失变量，确保 space 等 0 值项也在列表里
  tidyr::complete(abbreviation = fixed_order, fill = list(I.perc = 0)) %>%
  
  # 2. 重新关联 type 标签（防止 complete 丢掉类型）
 # left_join(distinct(keep0130[,c("abbreviation", "type")]), by = "abbreviation") %>%
  
  # 3. ✨ 核心：锁定 Factor Level
  mutate(abbreviation = factor(abbreviation, levels = fixed_order)) %>%
  
  # 4. 按照这个固定的顺序排序，生成固定的位置编号
  arrange(abbreviation) %>%
  mutate(petal = as.numeric(abbreviation),
         theta0 = petal * petal_angle) %>%
  
  # 5. 生成贝塞尔曲线所需的坐标点
  reframe(theta = theta0 + c(0, -petal_angle/2,  0, petal_angle/2, 0),
          r = I.perc^0.4 * c(0, 0.6, 1, 0.6, 0), 
          .by = c(abbreviation, I.perc, petal, theta0, type)) %>%
  
  # 6. 统一比例尺：你可以手动设一个最大值，比如 100^0.4
  mutate(max_r = max(r, na.rm = TRUE)) 

# 生成对应的标签数据
label_data <- g1rose.plot %>%
  filter(I.perc > 0) %>%
  group_by(abbreviation) %>%
  slice_max(r, n = 1, with_ties = FALSE) %>%
  ungroup()

gd1 <- ggplot(g1rose.plot, aes(theta, r, group = petal, fill = type)) +
  # 绘制虚线背景轴
  geom_segment(data = g1rose.plot%>% filter(I.perc > 0),
               aes(x = theta0, xend = theta0, y = 0, yend = r), 
               inherit.aes = FALSE, color = "grey50", linetype = "dotted", linewidth = 0.5) +
  
  ggforce::stat_bspline(geom = "area", n = 1000) +
  
  # 标签位置锁定
  geom_text(data = label_data, 
            aes(x = theta0, y = r, label = abbreviation), 
            size = 4.5, color = "black", fontface = "bold") +
  
  scale_fill_manual(values = my_npg_colors) +
  
  # start = 0 通常是正右方或正上方，你可以微调这个值 (单位是弧度)
  coord_radial(start = 0, expand = FALSE) + 
  
  theme_void() +
  theme(aspect.ratio = 0.98, legend.position = "none")
gd4
gd1
gd2
gd5
gd3
#不固定顺序####
#left_join(distinct(keep0130[,c(4,5)]),by ="abbreviation",keep = TRUE) %>% 
g5rose.plot<-g5rose %>% 
  dplyr::select(guild,abbreviation,I.perc,type) %>% 
  dplyr::arrange(type, I.perc) %>%  
 dplyr:: mutate(petal = row_number(),
         theta0 = petal * petal_angle)  %>% 
  reframe(theta = theta0 + c(0, -petal_angle/2,  0, 
                             petal_angle/2, 0),
          r = I.perc^0.4 * c(0, 0.6, 1, 0.6, 0), 
          .by = c(abbreviation, I.perc, petal, theta0)) %>% 
  left_join(distinct(keep0130[,c(4,5)]),by ="abbreviation") %>% 
  #dplyr::arrange(I.perc,descring=TRUE) %>% 
  mutate(max_r=max(r))
label_data <- g5rose.plot %>%
  # ✨ 修复 1：确保 abbreviation 是包含所有水平的 Factor
  mutate(abbreviation = as.factor(abbreviation)) %>%
  # ✨ 核心修复：强制保留所有可能的组合，即使数据框里没有对应的行
  group_by(abbreviation, .drop = FALSE) %>% 
  # ✨ 关键：添加 with_ties = FALSE，强制只取一行
  slice_max(r, n = 1, with_ties = FALSE) %>%
  ungroup()
# 将 NPG 的前 5 个颜色与你的具体分类严格绑定
my_npg_colors <- c(
  "Climate"                   = "#E64B35FF", # 红色
  "Geo-distance"              = "#4DBBD5FF", # 浅蓝
  "Soil properties"           = "#00A087FF", # 绿色
  "Topography"                = "#3C5488FF", # 深蓝
  "Vegetation & Productivity" = "#F39B7FFF"  # 橘粉色
)

gd5<-ggplot(g5rose.plot, aes(theta, r, group = petal, fill = type)) +
  geom_segment(data =g5rose.plot, 
               aes(x = theta0, xend = theta0, y = 0, yend = max_r-1),
               inherit.aes = FALSE, color = "grey80", linetype = "dotted", linewidth = 0.3) +
  ggforce::stat_bspline(geom = "area", n = 1000) +
 
  # 🌟 图层 3：花瓣顶部的文字标签
  geom_text(data = label_data, 
            aes(x = theta0, y = r, label=abbreviation),
                #label = paste0(abbreviation,"\n(",I.perc,"%)")), 
            size = 5, color = "black",lineheight = 0.8) +
  
 # geom_text(data = label_data, #aes(label=abbreviation))+
  #          aes(label = paste0(abbreviation,"\n(",I.perc,"%)"))) +
  labs(fill = NULL) +
  scale_fill_manual(values = my_npg_colors)+
#  ggsci::scale_fill_npg() + #指定Nature配色
  coord_radial() +
  theme_void() +theme(aspect.ratio = 0.99,legend.position = "none")
fd.driver
Td.driver
gd1
gd2
gd3
gd4
gd5
gdall<-plot_grid(gd1,gd2,gd3,gd4,gd5,nrow = 1)
gdall
TD.FD.driver<-plot_grid(fd.driver,Td.driver,nrow = 1)
ggsave(gdall,filename="/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/guild.driver/driverflower/gdall0328_noIPER.pdf",
       height = 6,width = 9,bg="transparent")
ggsave(gd1,filename="/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0//new_MAG/arc/01.arcinfo/F3/FD/guild.driver/driverflower/grose1_fix_nolabel_r.png",
       height = 6,width = 10,bg="transparent")
 ggsave(gd2,filename="/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/guild.driver/driverflower/grose2x_fix_nolabel_r.png",
       height = 6,width = 10,bg="transparent")
ggsave(gd3,filename="/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/guild.driver/driverflower/grose3_fix_nolabel_r.png",
       height = 6,width = 8,bg="transparent")
ggsave(gd4,filename="/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/guild.driver/driverflower/grose4_fix_nolabel_r.png",
       height = 6,width = 8,bg="transparent")
ggsave(gd5,filename="/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/guild.driver/driverflower/grose5_fix_nolabel_r.png",
       height = 6,width = 8,bg="transparent")

library(dplyr)
library(ggplot2)
library(ggforce)

library(dplyr)
library(ggplot2)
library(ggforce)
library(patchwork)

