#物种丰度
arc.mag.abundance <-
  readRDS(
    "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F4/arc.mag.abundance.rds"
  )
data <- arc.mag.abundance
rownames(data) <- NULL
data <- data %>% column_to_rownames(var = "Genome")
data <- sapply(data, as.numeric)
data <-
  data[, colSums(data > 0 |
                   is.na(data), na.rm = TRUE) != nrow(data)]
data <- apply(data, 2, function(x)
  x / sum(x))#归一化处理，每个样本中物种丰度的总和为1
rownames(data) <- arc.mag.abundance$Genome
paddy <-
  read.csv(
    "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/sampleinfo/sample.paddy.soil_filter1.csv"
  )

library(tidyverse)
library(vegan)

# MAG abundance matrix
otu <- t(data)

# 保留有metadata的sample
otu <- otu[rownames(otu) %in% s1meta$Sample1, ]

# metadata
meta <- s1meta %>%
 # filter(!Ecosystem %in% c("Bare Land","Shrubland")) %>%
  mutate(type = ifelse(
    Ecosystem %in% c("Paddy soil", "Wetland"),
    "Flooded Ecosystem",
    "Upland Ecosystem"
  ))

# 保证顺序一致
meta <- meta[match(rownames(otu), meta$Sample1), ]


rownames(meta) <- meta$Sample1

#meta<-meta %>% filter(!Ecosystem %in% c("Bare Land","Shrubland"))

#otu<-otu[rownames(otu) %in%meta$Sample1,]

#meta <- meta[match(rownames(otu), meta$Sample1), ]


# 检查
all(rownames(otu) == rownames(meta))

beta.bc <- vegdist(otu,
                   method = "bray")


set.seed(123)

nmds <- metaMDS(otu,
                distance = "bray",
                k = 2,
                trymax = 100)

nmds$stress


nmds_df <- as.data.frame(scores(nmds, display = "sites"))

nmds_df$Sample <- rownames(nmds_df)


nmds_df <- nmds_df %>%
  left_join(meta %>%
              select(Sample1, type, Ecosystem),
            by = c("Sample" = "Sample1"))
adonis_mag <- adonis2(
  beta.bc ~ Ecosystem,
  data=meta,
  permutations=999
)

adonis_mag

ggplot(
  nmds_df, #%>% #filter(!Ecosystem%in%c("Bare Land","Shrubland")),
  aes(
    NMDS1,
    NMDS2,
    color=Ecosystem
  )
)+
  xlim(-5,5)+ylim(-5,5)+
  geom_point(
    size=3,
    alpha=0.4
  )+
 # stat_ellipse(
  #  aes(group=type),
   # level=0.95
  #)+
  scale_color_manual(
    values=eco_color_map
  )+
  theme_classic()+
  labs(
    color=NULL,
    x="NMDS1",
    y="NMDS2"
  )+
  annotate(
    "text",
    x=0.02,
    y=0.02,
    label=paste0(
      "PERMANOVA\n",
      "R² = ",
      round(adonis_mag$R2[1],3),
      "\nP < 0.001\n",
      "Stress = ",
      round(nmds$stress,3)
    ),
    hjust=0,
    size=5
  )

#beta commun####
beta_disp_mag <- betadisper(
  beta.bc,
  meta$Ecosystem
)
anova(beta_disp_mag)

disp_mag_df1 <- data.frame(
  Sample = names(beta_disp_mag$distances),
  Distance = beta_disp_mag$distances
) 


disp_mag_df1 <- disp_mag_df1 %>%
  left_join(
    meta %>%
      select(Sample1,type,Ecosystem),
    by=c("Sample"="Sample1")
  )%>% filter(!Ecosystem%in%c("Bare Land","Shrubland"))
disp_mag_df1$type <- factor(
  disp_mag_df1$type,
  levels=c(
    "Upland Ecosystem",
    "Flooded Ecosystem"
  )
)

eco_order_mag <- disp_mag_df1 %>%
  group_by(type, Ecosystem) %>%
  summarise(
    mean_distance=mean(Distance),
    .groups="drop"
  ) %>%
  arrange(type, mean_distance)


disp_mag_df1$Ecosystem <- factor(
  disp_mag_df1$Ecosystem,
  levels=c(
    eco_order_mag$Ecosystem
  )
)

upland_mag_disp <- disp_mag_df1 %>%
  filter(type=="Upland Ecosystem")


kruskal.test(
  Distance ~ Ecosystem,
  data=upland_mag_disp
)

flood_mag_disp <- disp_mag_df1 %>%
  filter(type=="Flooded Ecosystem")


wilcox.test(
  Distance ~ Ecosystem,
  data=flood_mag_disp
)
sig_mag_df <- data.frame(
  type = c(
    "Upland Ecosystem",
    "Flooded Ecosystem"
  ),
  label = c(
    "Kruskal-Wallis\nP < 0.001",
    "Wilcoxon\nP = 0.003"
  ),
  x = c(
    3,
    1.5
  ),
  y = c(
    max(disp_mag_df1$Distance)*1.15,
    max(disp_mag_df1$Distance)*1.15
  )
)
sig_mag_df$type <- factor(
  sig_mag_df$type,
  levels=c(
    "Upland Ecosystem",
    "Flooded Ecosystem"
  )
)

#mag paired####
run_pair_adonis <- function(group1, group2){
  
  sub_meta <- meta %>%
    filter(
      Ecosystem %in% c(group1, group2)
    )
  
  
  sub_otu <- otu[
    rownames(otu) %in% sub_meta$Sample1,
  ]
  
  
  sub_meta <- sub_meta[
    match(
      rownames(sub_otu),
      sub_meta$Sample1
    ),
  ]
  
  
  sub_beta <- vegdist(
    sub_otu,
    method="bray"
  )
  
  
  result <- adonis2(
    sub_beta ~ Ecosystem,
    data=sub_meta,
    permutations=999
  )
  
  
  data.frame(
    Comparison=paste(
      group1,
      "vs",
      group2
    ),
    R2=result$R2[1],
    F=result$F[1],
    P=result$`Pr(>F)`[1]
  )
}
other_ecosystems <- c(
  "Forest",
  "Grassland",
  "Wetland",
  #"Bare Land",
  "Tundra",
  "Artificial Surfaces"
  #"Shrubland"
)


mag_ag_results <- lapply(
  other_ecosystems,
  function(x){
    run_pair_adonis(
      "Agricultural Land",
      x
    )
  }
) %>% unlist()


mag_ag_results
sig_mag_pair <- data.frame(
  group1 = c(
    "Agricultural Land",
    "Agricultural Land",
    "Agricultural Land",
    "Agricultural Land"
  ),
  group2 = c(
    "Forest",
    "Grassland",
    "Tundra",
    "Artificial Surfaces"
  ),
  y.position = c(
    0.80,
    0.85,
    0.89,
    0.95
  ),
  label = c(
    "***",
    "***",
    "***",
    "***"
  ),
  type = "Upland Ecosystem"
)
sig_mag_pair$type<-factor(
  sig_mag_pair$type,
  levels=c(
    "Upland Ecosystem")
)
library(ggpubr)
library(ggeasy)
ggplot(
  disp_mag_df1,
  aes(
    Ecosystem,
    Distance,
    fill=Ecosystem,
    color=Ecosystem
  )
)+
  geom_boxplot(
    width=0.6,
    outlier.shape=NA,
    alpha=0.5
  )+
  facet_grid(
    ~type,
    scales="free_x",
    space="free_x"
  )+
  scale_color_manual(
    values=eco_color_map
  )+
  scale_fill_manual(
    values=eco_color_map
  )+geom_text(
    data=sig_mag_df,
    aes(
      x=x,
      y=y,
      label=label
    ),
    inherit.aes=FALSE,
    size=4
  )+stat_pvalue_manual(
    sig_mag_pair%>% 
      filter(type=="Upland Ecosystem"),
    label="label",
    xmin="group1",
    xmax="group2",
    y.position="y.position",
    tip.length=0.01,
    size=5,
    inherit.aes=FALSE
  )+
  theme_classic()+
  theme(
    axis.text.x=element_text(
      angle=45,
      hjust=1
    ),
    strip.background=element_blank(),
    legend.position="none"
  )+
  labs(
    x=NULL,
    y="(Beta-dispersion) Distance to centroid"
  )+easy_all_text_size(size=12)+easy_all_text_color(color="black")
ggsave("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/ok/f3beta_community.pdf",height = 6,width = 8,bg="transparent")


#functional####
library(tidyverse)
library(vegan)
library(ape)
library(ggplot2)

ko.beta <- as.matrix(cluster_KO_all)

ko_rel <- sweep(
  ko.beta,
  1,
  rowSums(ko.beta),
  "/"
)
ko_beta <- vegdist(
  ko_rel,
  method="bray"
)
# NMDS
set.seed(123)

ko_nmds <- metaMDS(
  ko_rel,
  distance = "bray",
  k = 2,
  trymax = 100
)

# 查看stress
ko_nmds$stress

ko_nmds_df <- as.data.frame(scores(ko_nmds,
                                   display = "sites"))

ko_nmds_df$Sample <- rownames(ko_nmds_df)

meta_ko <- s1meta %>%
  filter(Sample1 %in% ko_nmds_df$Sample)
meta_ko <- meta_ko[match(ko_nmds_df$Sample,
                         meta_ko$Sample1),]



ko_nmds_df <- ko_nmds_df %>%
  left_join(meta_ko,
            by = c("Sample" = "Sample1")) %>% 
  filter(!Ecosystem%in%c("Bare Land","Shrubland"))

ggplot(ko_nmds_df,
         #filter(!Ecosystem%in%c("Bare Land","Shrubland")),
       aes(NMDS1,
           NMDS2,
           color = Ecosystem)) +
  geom_point(size = 3,
             alpha = 0.6) +scale_color_manual(values = eco_color_map)+
  annotate(
    "text",
    x = max(ko_nmds_df$NMDS1)-0.5,
    y = max(ko_nmds_df$NMDS2),
    label = paste0(
      "PERMANOVA\n",
      "R² = 0.131\n",
      "P < 0.001\n",
      "Stress = 0.136"
    ),
    hjust=0,
    size=5
  )+
 # stat_ellipse(aes(group = Ecosystem),
              # level = 0.95) +
  theme_classic() +
  easy_all_text_size(size = 15)+
  easy_all_text_color(color = "black")+
  labs(x = paste0("NMDS1"),
       y = "NMDS2")
ggsave("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/ok/f3beta_nmds1.pdf",height = 5,width = 8,bg="transparent")

adonis_ko <- adonis2(
  ko_beta ~ Ecosystem,
  data=meta_ko,
  permutations=999
)

adonis_ko
s1meta<-read.csv('/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/paper/写作参考/nature/submit/Supplementary Table/Supplementary Table 1.csv')
meta_ko <- s1meta %>%
  filter(
    Sample1 %in% ko_pcoa_df$Sample
  )
meta_ko <- meta_ko[
  match(
    rownames(ko_rel),
    meta_ko$Sample1
  ),
]
adonis_ko <- adonis2(
  ko_beta ~ Ecosystem,
  data=meta_ko,
  permutations=999
)

adonis_ko



#paired####


library(dplyr)
library(vegan)


# 筛选 upland / non-flooded ecosystems
library(vegan)
library(dplyr)


run_pair_adonis <- function(group1, group2){
  
  # metadata筛选
  sub_meta <- meta_ko %>%
    filter(
      Ecosystem %in% c(group1, group2)
    )
  
  
  # 提取对应KO矩阵
  sub_otu <- ko_rel[
    rownames(ko_rel) %in% sub_meta$Sample1,
  ]
  
  
  # 保证顺序一致
  sub_meta <- sub_meta[
    match(
      rownames(sub_otu),
      sub_meta$Sample1
    ),
  ]
  
  
  # Bray-Curtis
  sub_beta <- vegdist(
    sub_otu,
    method="bray"
  )
  
  
  # PERMANOVA
  result <- adonis2(
    sub_beta ~ Ecosystem,
    data=sub_meta,
    permutations=999
  )
  
  
  # 提取结果
  data.frame(
    Comparison=paste(
      group1,
      "vs",
      group2
    ),
    R2=result$R2[1],
    F=result$F[1],
    P=result$`Pr(>F)`[1]
  )
}
other_ecosystems <- c(
  "Forest",
  "Grassland",
  "Wetland",
 # "Bare Land",
  "Tundra",
  "Artificial Surfaces"
  #"Shrubland"
)

ag_results <- lapply(
  other_ecosystems,
  function(x){
    run_pair_adonis(
      "Agricultural Land",
      x
    )
  }
)


ag_results <- bind_rows(
  ag_results
)


ag_results


library(vegan)
library(dplyr)
library(ggplot2)


# 定义系统
meta_ko <- meta_ko %>%
  mutate(
    type = ifelse(
      Ecosystem %in% c("Paddy soil", "Wetland"),
      "Flooded Ecosystem",
      "Upland Ecosystem"
    )
  )


# beta dispersion
beta_disp <- betadisper(
  ko_beta,
  meta_ko$Ecosystem
)


# 提取distance
disp_df1 <- data.frame(
  Sample = meta_ko$Sample1,
  Distance = beta_disp$distances,
  Ecosystem = meta_ko$Ecosystem,
  type = meta_ko$type
) %>% filter(!Ecosystem%in%c("Bare Land","Shrubland"))


disp_df1$type <- factor(
  disp_df1$type,
  levels=c(
    "Upland Ecosystem",
    "Flooded Ecosystem"
  )
)


library(dplyr)

ecosystem_order <- disp_df1 %>%
  group_by(Ecosystem) %>%
  summarise(
    mean_distance = mean(Distance, na.rm = TRUE)
  ) %>%
  arrange(mean_distance) %>%
  pull(Ecosystem)


disp_df1$Ecosystem <- factor(
  disp_df1$Ecosystem,
  levels = ecosystem_order
)

upland_disp <- disp_df1 %>%
  filter(type == "Upland Ecosystem")


anova_upland <- aov(
  Distance ~ Ecosystem,
  data = upland_disp
)

kruskal.test(
  Distance ~ Ecosystem,
  data = upland_disp
)
summary(anova_upland)

flood_disp <- disp_df1 %>%
  filter(type == "Flooded Ecosystem")


wilcox.test(
  Distance ~ Ecosystem,
  data = flood_disp
)


library(ggplot2)

sig_df <- data.frame(
  type = c(
    "Upland Ecosystem",
    "Flooded Ecosystem"
  ),
  label = c(
    "Kruskal-Wallis\nP < 0.001",
    "Wilcoxon\nP = 0.001"
  ),
  x = c(
    3,
    1.5
  ),
  y = c(
    max(disp_df1$Distance)*1.15,
    max(disp_df1$Distance)*1.15
  )
)
sig_df$type<-factor(
  sig_df$type,
  levels = c(
    "Upland Ecosystem",
    "Flooded Ecosystem"
  )
)

disp_df1$type <- factor(
  disp_df1$type,
  levels = c(
    "Upland Ecosystem",
    "Flooded Ecosystem"
  )
)

library(FSA)

upland_disp <- disp_df1 %>%
  filter(type=="Upland Ecosystem")


dunn_upland <- dunnTest(
  Distance ~ Ecosystem,
  data=upland_disp,
  method="bonferroni"
)

dunn_upland


sig_pair <- data.frame(
  group1 = c(
    "Agricultural Land",
    "Agricultural Land",
    "Agricultural Land"
  ),
  group2 = c(
    "Forest",
    "Grassland",
    "Artificial Surfaces"
  ),
  p.adj = c(
    0.007255312,
    1.410065e-09,
    1.996035e-18
  ),
  y.position = c(
    0.45,
    0.55,
    0.65
  ),
  type = "Upland Ecosystem"
)


sig_pair$label <- c(
  "**",
  "***",
  "***"
)
sig_pair$type<-factor(
  sig_pair$type,
  levels = c(
    "Upland Ecosystem",
    "Flooded Ecosystem"
  )
)

ggplot(disp_df1,
       aes(Ecosystem,
           Distance,
           fill = Ecosystem,
           color = Ecosystem)) +
  geom_boxplot(width = 0.6,
               outlier.shape = NA,
               alpha = 0.5) +
  
  facet_grid(~ type,
             scales = "free_x",
             space = "free_x") +
  
  geom_text(
    data = sig_df,
    aes(x = x,
        y = y,
        label = label),
    inherit.aes = FALSE,
    size = 4
  ) +
  
  stat_pvalue_manual(
    sig_pair,
    label = "label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    size = 5,
    inherit.aes = FALSE
  ) +
  
  scale_color_manual(values = eco_color_map) +
  scale_fill_manual(values = eco_color_map) +
  
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45,
                               hjust = 1),
    strip.background = element_blank(),
    legend.position = "none"
  ) +
  labs(x = NULL,
       y = "(Beta-dispersion) Distance to centroid") +
  easy_all_text_size(size = 15) +
  easy_all_text_color(color = "black")

ggsave("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/ok/f3beta_fun1.pdf",height = 6,width = 8,bg="transparent")


library(tidyverse)
library(ggeasy)


# taxonomic beta dispersion
tax_df <- disp_mag_df1 %>%
  select(
    Ecosystem,
    type,
    Distance
  ) %>%
  mutate(
    Metric="Taxonomic composition"
  )


# functional beta dispersion
fun_df <- disp_df1 %>%
  select(
    Ecosystem,
    type,
    Distance,
  ) %>%
  mutate(
    Metric="Functional composition"
  )


# 合并
beta_all <- bind_rows(
  tax_df,
  fun_df
)


# factor顺序
beta_all$type <- factor(
  beta_all$type,
  levels=c(
    "Upland Ecosystem",
    "Flooded Ecosystem"
  )
)


beta_all$Metric <- factor(
  beta_all$Metric,
  levels=c(
    "Taxonomic composition",
    "Functional composition"
  )
)
metric_colors <- c(
  "Taxonomic composition"="#D34564",
  "Functional composition"="#396095"
)

ggplot(
  beta_all,
  aes(
    x=Ecosystem,
    y=Distance,
    fill=Metric,
    color=Metric
  )
)+
  
  
  geom_boxplot(
    width=0.35,
    linewidth=0.35,
    notch=FALSE,
    notchwidth=0,
    outlier.shape=NA,
    alpha=0.65,
    position=position_dodge(width=0.7)
  )+
  
  
  facet_wrap(
    ~type,
    scales="free_x",
    space="free_x"
  )+
  
  
  scale_fill_manual(
    values=metric_colors
  )+
  
  scale_color_manual(
    values=metric_colors
  )+
  
  
  theme_bw()+
  
  theme(
    panel.grid=element_blank(),
    
    legend.position="top",
    
    legend.title=element_blank(),
    
    strip.background=
      element_rect(
        fill="gray90",
        color="black"
      ),
    
    strip.text=
      element_text(
        size=12,
        face="bold"
      ),
    
    axis.text.x=
      element_text(
        angle=45,
        hjust=1,
        size=12,
        color="black"
      ),
    
    axis.text.y=
      element_text(
        size=12,
        color="black"
      ),
    
    axis.title.y=
      element_text(
        size=15
      )
  )+
  
  
  labs(
    x=NULL,
    y="Beta-dispersion\n(Distance to centroid)"
  )+
  
  easy_all_text_size(size=18)



# 合并数据
plot_beta <- bind_rows(
  disp_mag_df1 %>%
    mutate(Metric="Taxonomic composition"),
  
  disp_df1 %>%
    mutate(Metric="Functional composition")
)


# 计算每个生态系统两个指标的median，用于排序
order_df <- plot_beta %>%
  group_by(
    type,
    Ecosystem,
    Metric
  ) %>%
  summarise(
    med = mean(Distance, na.rm=TRUE),
    .groups="drop"
  ) %>%
  pivot_wider(
    names_from = Metric,
    values_from = med
  ) %>%
  mutate(
    diff =
      `Taxonomic composition` -
      `Functional composition`
  ) %>%
  arrange(desc(diff))


order_levels <- order_df$Ecosystem

plot_beta$Ecosystem <- factor(
  plot_beta$Ecosystem,
  levels = unique(order_levels)
)
plot_beta$Metric<-factor( plot_beta$Metric,levels = c( "Taxonomic composition",
                                                       "Functional composition"))
ggplot(
  plot_beta,
  aes(
    Ecosystem,
    Distance,
    color=Metric
  )
)+ geom_boxplot(
  width=0.35,
  linewidth=0.35,
  notch=FALSE,
  notchwidth=0,
  outlier.shape=NA,
  alpha=0.65,
  position=position_dodge(width=0.7)
)+
  
  
  # geom_point(
  #  width=0.15,
  # size=1.5,
  #alpha=0.4
  #)+
  
  stat_summary(
    fun=mean,
    geom="point",
    size=5,
    shape=18,
    position=position_dodge(width=0.7)
  )+
  
  facet_grid(
    ~type,
    scales="free_x",
    space="free_x"
  )+
  
  scale_color_manual(
    values=c(
      "Taxonomic composition"="#D34564",
      "Functional composition"="#396095"
    )
  )+
  
  theme_classic()+
  
  theme(
    axis.text.x=
      element_text(
        angle=45,
        hjust=1
      ),
    strip.background=element_blank(),
    legend.position="top"
  )+
  
  labs(
    x=NULL,
    y="Beta-dispersion\n(Distance to centroid)",
    color=NULL
  )+easy_all_text_size(size=12)+easy_all_text_color(color = "black")


ggplot(plot_beta,
       aes(Ecosystem,
           Distance,
           color = Metric)) +
  stat_summary(
    fun.data = median_hilow,
    geom = "errorbar",
    width = 0.1,
    position = position_dodge(width = 0.6)
  ) +
  stat_summary(fun = mean,
               geom = "segment",
               aes(xend = after_stat(x),
                   yend = after_stat(y))) + stat_summary(
                     fun = mean,
                     geom = "point",
                     size = 4,
                     shape=18,
                     position = position_dodge(width = 0.6)
                   ) + scale_color_manual(values = c(
                     "Taxonomic composition" = "#D34564",
                     "Functional composition" = "#396095"
                   )) +  facet_grid(~ type,
                                    scales = "free_x",
                                    space = "free_x") +
  theme_classic() +
  
  theme(
    axis.text.x =
      element_text(angle = 45,
                   hjust = 1),
    strip.background = element_blank(),
    legend.position = "top"
  ) +
  
  labs(x = NULL,
       y = "Beta-dispersion\n(Distance to centroid)",
       color = NULL) + easy_all_text_size(size = 10) + easy_all_text_color(color = "black")


ggsave("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/f3b_boxplotall_supply1.pdf", height = 4, width = 10, bg = "transparent")


#
library(tidyverse)
library(ggeasy)


#===========================
# Taxonomic diversity
#===========================

tax_df <- shannon_guild4 %>%filter(!Ecosystem.y%in%c("Bare Land","Shrubland")) %>% 
  select(
    Ecosystem = Ecosystem.y,
    Value = Shannon.TD
  )


# 按mean排序
tax_order <- tax_df %>%
  group_by(Ecosystem) %>%
  summarise(
    mean_value = mean(Value, na.rm=T)
  ) %>%
  arrange(desc(mean_value)) %>%
  pull(Ecosystem)


tax_df$Ecosystem <- factor(
  tax_df$Ecosystem,
  levels = tax_order
)



ggplot(
  tax_df,
  aes(
    Ecosystem,
    Value,
    color=Ecosystem
  )
)+
  
  # mean ± 95% CI
  stat_summary(
    fun.data = mean_cl_normal,
    geom="errorbar",
    width=0.15,
    linewidth=0.8
  )+
  
  stat_summary(
    fun=mean,
    geom="point",
    shape=18,
    size=4
  )+scale_color_manual(values =eco_color_map )+
  
  
  theme_classic()+
  
  theme(aspect.ratio = 0.5,
    axis.text.x =
      element_text(
        angle=45,
        hjust=1
      ),legend.position="none"
  )+
  
  labs(
    x=NULL,
    y="Taxonomic diversity\n(Shannon index)"
  )+
  
  easy_all_text_size(size=15)+
  easy_all_text_color(color="black")

ggsave("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/ok/f3_td.pdf",height = 6,width = 8,bg="transparent")

#===========================
# Functional guild diversity
#===========================


fun_df <- shannon_guild4 %>%filter(!Ecosystem.y%in%c("Bare Land","Shrubland")) %>% 
  select(
    Ecosystem = Ecosystem.y,
    Value = Menhinick_Coef
  )


# 按mean排序

fun_order <- fun_df %>%
  group_by(Ecosystem) %>%
  summarise(
    mean_value = mean(Value, na.rm=T)
  ) %>%
  arrange(desc(mean_value)) %>%
  pull(Ecosystem)


fun_df$Ecosystem <- factor(
  fun_df$Ecosystem,
  levels = fun_order
)



ggplot(
  fun_df,
  aes(
    Ecosystem,
    Value,color=Ecosystem
  )
)+
  
  
  stat_summary(
    fun.data = mean_cl_normal,
    geom="errorbar",
    width=0.15,
    linewidth=0.8
  )+
  
  
  stat_summary(
    fun=mean,
    geom="point",
    shape=18,
    size=4
  )+
  
  
  theme_classic()+scale_color_manual(values =eco_color_map )+
  
  
  theme(aspect.ratio = 0.5,
    axis.text.x =
      element_text(
        angle=45,
        hjust=1
      ),legend.position="none"
  )+
  
  labs(
    x=NULL,
    y="Functional guild diversity\n(Menhinick index)"
  )+
  
  easy_all_text_size(size=15)+
  easy_all_text_color(color="black")

ggsave("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/FDdriver/ok/f3_gfm.pdf",height = 6,width =8 ,bg="transparent")

