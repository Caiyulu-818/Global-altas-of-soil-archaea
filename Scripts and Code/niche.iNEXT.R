#niche.####
library(vegan)
library(spaa)
library(ggplot2)
library(gghalves)
library(tidyr)

# 输入 OTU/MAG abundance matrix
otu <- arc.mag.abundance %>%
  column_to_rownames(var="Genome")


## =========================
## 1. Niche overlap
## =========================

# 生态位重叠指数
# 方法：
# levins
# schoener
# morisita
# pianka
# petraitis
# czech

over <- niche.overlap(
  otu_table1,
  method = "levins"
)


## Bootstrap niche overlap

nob1 <- niche.overlap.boot(
  otu_table1,
  method = "pianka"
)


## =========================
## 2. Niche width
## =========================

# 计算生态位宽度
niche_width <- niche.width(
  mat = t(otu_table),
  method = "levins"
) %>%
  t() %>%
  as.data.frame()


## 处理 NA 和 Inf

niche_width_clean <- niche_width %>%
  rownames_to_column("ID") %>%
  mutate(
    V1 = ifelse(
      is.infinite(V1) | is.na(V1),
      NA,
      V1
    )
  ) %>%
  mutate(
    niche_width_norm = V1 / ncol(otu_table)
  ) %>%
  drop_na(niche_width_norm)


## 添加 cluster 信息

niche_width_clean <- merge(
  niche_width_clean,
  arc1645,
  by.x="ID",
  by.y="genome",
  all.y=TRUE
)

colnames(niche_width_clean)[2] <- "niche_width index"


## =========================
## 3. Specialist / Generalist 分类
## =========================

# 四分位数划分

q25 <- quantile(
  niche_width_clean$niche_width_norm,
  0.25,
  na.rm = TRUE
)

q75 <- quantile(
  niche_width_clean$niche_width_norm,
  0.75,
  na.rm = TRUE
)


niche_width_labeled <- niche_width_clean %>%
  mutate(
    Category = case_when(
      niche_width_norm <= q25 ~ "Specialist",
      niche_width_norm >= q75 ~ "Generalist",
      TRUE ~ "Not significant"
    )
  )


write.csv(
  niche_width_labeled,
  "niche_width_labeled.csv"
)


## =========================
## 4. Cluster 间 niche width 比较
## =========================

niche_width_labeled$Cluster <- paste0(
  "Cluster ",
  niche_width_labeled$Cluster
)


# cluster 顺序

cluster_orderniche <- niche_width_labeled %>%
  group_by(Cluster) %>%
  summarise(
    mean_niche =
      mean(
        niche_width_norm,
        na.rm = TRUE
      )
  ) %>%
  arrange(desc(mean_niche)) %>%
  pull(Cluster)


# Kruskal-Wallis

kruskal_niche <- kruskal.test(
  niche_width_norm ~ Cluster,
  data = niche_width_labeled
)


## Violin plot

ggplot(
  niche_width_labeled,
  aes(
    x=factor(Cluster,
             cluster_orderniche),
    y=log(niche_width_norm),
    fill=Cluster
  )
)+
geom_violin(
  width=0.5,
  color=NA
)+
geom_boxplot(
  width=0.5,
  color=NA
)+
geom_jitter(
  width=0.15,
  size=1.5,
  alpha=0.1
)+
theme_bw()


## =========================
## 5. Specialist / Generalist 分类组成
## Sankey
## =========================

library(ggsankey)


niche_tax <- niche_width_labeled %>%
  group_by(
    Cluster,
    phylum,
    class,
    Category
  ) %>%
  summarise(
    count=n()
  )


niche_tax_sankey <- niche_width_labeled %>%
  make_long(
    Cluster,
    class,
    Category
  )


ggplot(
  niche_tax_sankey,
  aes(
    x=x,
    next_x=next_x,
    node=node,
    next_node=next_node,
    fill=node,
    label=node
  )
)+
geom_alluvial(
  flow.alpha=0.6
)+
geom_alluvial_text(
  size=3
)+
theme_alluvial()


## =========================
## 6. Specialist / Generalist 分类组成
## Donut plot
## =========================

library(ggforce)
library(patchwork)


donut_by_class <- function(
    df,
    cluster_name,
    category_name,
    inner_r=0.52,
    outer_r=1
){

d <- df %>%
  filter(
    Cluster==cluster_name,
    Category==category_name
  ) %>%
  group_by(class) %>%
  summarise(
    Freq=sum(count)
  ) %>%
  arrange(desc(Freq))


d <- d %>%
  mutate(
    frac=Freq/sum(Freq),
    start_angle=cumsum(frac)-frac,
    end_angle=cumsum(frac)
  )


ggplot(d)+
geom_arc_bar(
 aes(
 x0=0,
 y0=0,
 r0=inner_r,
 r=outer_r,
 start=start_angle*2*pi,
 end=end_angle*2*pi,
 fill=class
 )
)+
coord_fixed()+
theme_void()

}


# 每个 Cluster 绘制 Generalist + Specialist

clusters <- unique(niche_tax$Cluster)


plots <- lapply(
  clusters,
  function(cl){
    
    p1 <- donut_by_class(
      niche_tax,
      cl,
      "Generalist"
    )
    
    p2 <- donut_by_class(
      niche_tax,
      cl,
      "Specialist"
    )
    
    p1+p2
  }
)


final_plot <- wrap_plots(
  plots,
  ncol=2
)


#richness across cluster####

library(vegan)


## =========================
## 1. Observed richness
## =========================


S1_obs <- specnumber(
  t(otu_table1)
)

S1_obs <- data.frame(
  Sample=rownames(t(otu_table1)),
  Richness=as.vector(S1_obs)
)%>%
mutate(
  Cluster="Cluster 1"
)



S2_obs <- specnumber(
  t(otu_table2)
)

S2_obs <- data.frame(
  Sample=rownames(t(otu_table2)),
  Richness=as.vector(S2_obs)
)%>%
mutate(
  Cluster="Cluster 2"
)



S3_obs <- specnumber(
  t(otu_table3)
)

S3_obs <- data.frame(
  Sample=rownames(t(otu_table3)),
  Richness=as.vector(S3_obs)
)%>%
mutate(
  Cluster="Cluster 3"
)



S4_obs <- specnumber(
  t(otu_table4)
)

S4_obs <- data.frame(
  Sample=rownames(t(otu_table4)),
  Richness=as.vector(S4_obs)
)%>%
mutate(
  Cluster="Cluster 4"
)



S5_obs <- specnumber(
  t(otu_table5)
)

S5_obs <- data.frame(
  Sample=rownames(t(otu_table5)),
  Richness=as.vector(S5_obs)
)%>%
mutate(
  Cluster="Cluster 5"
)


richnessall <- rbind(
  S1_obs,
  S2_obs,
  S3_obs,
  S4_obs,
  S5_obs
)


## richness 分布比较

ggplot(
  richnessall,
  aes(
    x=Cluster,
    y=Richness
  )
)+
geom_violin(
  aes(fill=Cluster)
)+
geom_boxplot(
  width=0.4
)+
geom_jitter(
  width=0.2,
  alpha=0.2
)+
theme_bw()



# =========================
# 2. iNEXT rarefaction / extrapolation
# =========================

library(iNEXT)


## abundance matrix → incidence frequency

to_incidence_vec <- function(
    otu_tab,
    transpose=TRUE
){

mat <- if(transpose){
  t(otu_tab)
}else{
  as.matrix(otu_tab)
}


pa <- (mat>0)*1


T_units <- nrow(pa)


inc_freq <- colSums(pa)


c(
T_units,
inc_freq
)

}



## 构建 cluster 输入

incidence_list <- list(

"Cluster 1" =
to_incidence_vec(
otu_table1
),

"Cluster 2" =
to_incidence_vec(
otu_table2
),

"Cluster 3" =
to_incidence_vec(
otu_table3
),

"Cluster 4" =
to_incidence_vec(
otu_table4
),

"Cluster 5" =
to_incidence_vec(
otu_table5
)

)



## 设置统一外推终点

T_vec <- sapply(
incidence_list,
function(v)v[1]
)


endpoint_unified <- 
2*max(T_vec)



## iNEXT

out <- iNEXT(
incidence_list,
q=0,
datatype="incidence_freq",
se=TRUE,
conf=0.95,
nboot=300,
endpoint=endpoint_unified
)



## =========================
## 3. Rarefaction curve
## =========================


p_richness <- ggiNEXT(
out,
type=1
)+
xlab(
"Number of sampling units"
)+
ylab(
"Species richness (q = 0)"
)


ggsave(
"richness_iNEXT.pdf",
p_richness,
height=6,
width=8
)



## =========================
## 4. Coverage-based curve
## =========================


p_cov <- ggiNEXT(
out,
type=3
)+
xlab(
"Sample coverage"
)+
ylab(
"Species richness (q = 0)"
)


ggsave(
"richness_iNEXT_cov.pdf",
p_cov,
height=6,
width=8
)



## =========================
## 5. Chao2 richness estimation
## =========================

Chao2 <- estimateD(
  incidence_list,
  q=0,
  datatype="incidence_freq"
)


write.csv(
Chao2,
"Chao2_cluster.csv"
)