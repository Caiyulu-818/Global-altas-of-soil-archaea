############################################################
# Functional clustering and ecological characterization
# of soil archaeal MAGs
############################################################


############################
# 1. Functional annotation processing
############################

library(tidyverse)
library(cluster)
library(factoextra)
library(ape)


# Read KO functional annotation table
arc.f <- read.csv(
  "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/02.arcfunction/atc1645all.csv"
)


# Summarize KO counts for each genome
arc.f1 <- arc.f %>%
  group_by(Ecosystem, genome, count) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  mutate(genome = str_remove(genome, ".\\d+$")) %>%   # Remove bin suffix
  distinct()


# Add genome metadata
arc1322bin <- read.csv(
  "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/sampleinfo/arc1322sampleinfo的副本.csv"
)

arc.f1 <- merge(
  arc.f1,
  arc1322bin,
  by = "genome",
  all.x = TRUE
)


arc1322 <- read.csv(
  "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/sampleinfo/arcsampleinfo1322.csv"
)



############################
# 2. Sequencing depth and KO abundance visualization
############################


# Remove potential outliers
arc.f2_filtered <- arc.f2 %>%
  filter(Paired_seqs / 10000000 < 75 | count < 3)


# Relationship between sequencing depth and KO abundance
ggplot(
  na.omit(arc.f2),
  aes(
    x = Paired_seqs / 10000000,
    y = count / 1000,
    color = Ecosystem.x
  )
) +
  geom_point(
    size = 1,
    alpha = 0.6
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  labs(
    x = "Paired reads (million reads)",
    y = "KO count (thousand)"
  ) +
  theme_classic() +
  scale_color_manual(values = dc) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggsave(
  "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/arc_figure/f2/Sf3.ecosystem.pdf"
)



############################
# 3. Functional profile matrix construction
############################


# Construct genome × KO abundance matrix

arc.ko <- arc.f[, c(2, 21, 23)] %>%
  distinct() %>%
  group_by(genome, KO) %>%
  summarise(
    count = sum(count),
    .groups = "drop"
  ) %>%
  spread(
    key = KO,
    value = count,
    fill = 0
  ) %>%
  as.data.frame()


rownames(arc.ko) <- arc.ko[,1]

arc.ko <- arc.ko[,-1]



############################
# 4. Hierarchical clustering based on KO profiles
############################


# Calculate Euclidean distance
dist_matrix <- dist(
  arc.ko,
  method = "euclidean"
)


# Hierarchical clustering using Ward's method
hc <- hclust(
  dist_matrix,
  method = "ward.D2"
)


# Export clustering tree
hc.tree <- as.phylo(hc)

write.tree(
  hc.tree,
  "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/arc_figure/function.hc.tree4.tree"
)



############################
# 5. Determine optimal cluster number
############################


# Cut dendrogram into five functional clusters
clusters <- cutree(
  hc,
  k = 5
)



# Silhouette analysis
sil_width <- c()

for(k in 2:10){

  pam_fit <- pam(
    arc.ko.num,
    k = k
  )

  sil_width[k] <- pam_fit$silinfo$avg.width

}


plot(
  2:10,
  sil_width[2:10],
  type = "b",
  pch = 19,
  frame = FALSE,
  xlab = "Number of clusters",
  ylab = "Average silhouette width"
)



# Gap statistic analysis

set.seed(123)

gap_stat <- clusGap(
  arc.ko.num,
  FUN = hcut,
  K.max = 10,
  B = 50
)



# Maximum increment method

gap_values <- gap_stat$Tab[, "gap"]

increments <- diff(gap_values)

k_max_increase <- which.max(increments) + 1



# One standard error rule

SE_values <- gap_stat$Tab[, "SE.sim"]

gap_k <- gap_stat$Tab[, "gap"]

k_oneSE <- max(
  which(
    gap_k >= (
      max(gap_k) -
        SE_values[which.max(gap_k)]
    )
  )
)



fviz_gap_stat(gap_stat)



############################
# 6. Assign functional clusters
############################


arc1645$Cluster <- clusters

arc.ko$Cluster <- arc1645$Cluster



# Extract genomes belonging to each functional cluster

cluster1_genome <- arc.ko %>% filter(Cluster == 1)
cluster2_genome <- arc.ko %>% filter(Cluster == 2)
cluster3_genome <- arc.ko %>% filter(Cluster == 3)
cluster4_genome <- arc.ko %>% filter(Cluster == 4)
cluster5_genome <- arc.ko %>% filter(Cluster == 5)



############################
# 7. Extract KO profiles for each cluster
############################


gathered_ko_list <- list()


for(i in 1:5){

  cluster_data <- arc.ko %>%
    filter(Cluster == i)


  cluster_data_gathered <- cluster_data %>%
    select(-Cluster) %>%
    gather(
      key = "KO",
      value = "count",
      -genome
    ) %>%
    mutate(
      Cluster = i
    )


  gathered_ko_list[[paste0("Cluster_", i)]] <-
    cluster_data_gathered

}



final_gathered_ko <- bind_rows(
  gathered_ko_list
) %>%
  filter(count != 0)



# Cluster size

cluster.type <- final_gathered_ko %>%
  select(genome, Cluster) %>%
  distinct() %>%
  group_by(Cluster) %>%
  summarise(
    genome_number = n()
  )



############################
# 8. Functional characterization of clusters
############################


arc1645info <- all.ko.arc.info[, c(2:13)] %>%
  distinct()


final_gathered_ko_info <- merge(
  final_gathered_ko,
  arc1645info,
  by = "genome"
)



# Extract each functional guild

cluster1_ko <- final_gathered_ko_info %>% filter(Cluster == 1)
cluster2_ko <- final_gathered_ko_info %>% filter(Cluster == 2)
cluster3_ko <- final_gathered_ko_info %>% filter(Cluster == 3)
cluster4_ko <- final_gathered_ko_info %>% filter(Cluster == 4)
cluster5_ko <- final_gathered_ko_info %>% filter(Cluster == 5)



############################
# 9. Generate sample-level KO abundance matrices
############################


cluster5_KO_matrix <- cluster5_ko[, c(1:3,5)] %>%
  distinct() %>%
  group_by(Sample, KO) %>%
  summarise(
    count = sum(count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = KO,
    values_from = count,
    values_fill = 0
  ) %>%
  na.omit() %>%
  column_to_rownames("Sample")



write_csv(
  cbind(
    sample = rownames(cluster4_KO_matrix),
    cluster4_KO_matrix
  ),
  "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/cluster.diversity/cluster4_KO_matrix.csv"
)



############################
# 10. Taxonomic distribution of functional clusters
############################


arc.mag.abundance <- read.csv(
  "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/arc16/mag_mapping/arc.magabundance.csv"
)


arc.mag.abundance <- as.data.frame(
  t(arc.mag.abundance)
) %>%
  rownames_to_column(
    var = "Genome"
  )


colnames(arc.mag.abundance) <- arc.mag.abundance[2,]

arc.mag.abundance <- arc.mag.abundance[-c(1:2),]


arc.mag.abundance <- arc.mag.abundance[
  ,
  !duplicated(colnames(arc.mag.abundance))
]



# Merge abundance and functional cluster information

arc.mag.abundance1 <- merge(
  arc.mag.abundance,
  arc.ko[, c(1,3921)],
  by.x = "Genome",
  by.y = "genome"
)



write.csv(
  arc.mag.abundance1,
  "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/Rscript&sheet/arc.mag.abundance.cluster3.csv"
)



arc.mag.abundanceinfo <- merge(
  arc.mag.abundance1,
  arc1645info,
  by.x = "Genome",
  by.y = "genome"
)



cluster1_tax <- arc.mag.abundanceinfo %>% filter(Cluster == 1)
cluster2_tax <- arc.mag.abundanceinfo %>% filter(Cluster == 2)
cluster3_tax <- arc.mag.abundanceinfo %>% filter(Cluster == 3)
cluster4_tax <- arc.mag.abundanceinfo %>% filter(Cluster == 4)
cluster5_tax <- arc.mag.abundanceinfo %>% filter(Cluster == 5)



############################
# 11. Summarize taxonomic composition
############################


cluster_info <- final_gathered_ko_info[, c(1,4,10:15)] %>%
  distinct() %>%
  group_by(
    Cluster,
    phylum,
    class
  ) %>%
  summarise(
    count = n(),
    .groups = "drop"
  )



write.csv(
  final_gathered_ko_info_samplerename,
  "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/Rscript&sheet/finalgathered5_fun.cluster.csv"
)


write.csv(
  cluster_info,
  "/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/02.arcfunction/arc1645/function.cluster.csv"
)