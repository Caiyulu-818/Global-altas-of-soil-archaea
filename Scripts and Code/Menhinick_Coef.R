#物种丰度
arc.mag.abundance<-readRDS("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F4/arc.mag.abundance.rds")
data <-arc.mag.abundance
rownames(data)<-NULL
data<-data%>%column_to_rownames(var="Genome")
data <- sapply(data, as.numeric)
data <- data[, colSums(data > 0 | is.na(data), na.rm = TRUE) != nrow(data)]
data<- apply(data, 2, function(x)x/sum(x))#归一化处理，每个样本中物种丰度的总和为1
rownames(data) <- arc.mag.abundance$Genome

#功能矩阵
#table <- merge(Gia,data,by="genus")  #为了使gia和otu的genus名字对应起来
#table_fr<-merge(arc.ko,data,by.x="genome",by.y="Genome")
Gia_table <- arc.ko%>%column_to_rownames(var="genome") #Genes
# 删除所有行和为零或全NA的行
Gia_table <- Gia_table[rowSums(Gia_table > 0 | is.na(Gia_table), na.rm = TRUE) > 0, ]
# 删除所有列和为零或全NA的列
Gia_table <- Gia_table[, colSums(Gia_table > 0 | is.na(Gia_table), na.rm = TRUE) > 0]

#genome -guild 信息
arc1645_cluster<-arc1645 %>%as.data.frame() %>%  dplyr::select(genome,Cluster)

Num_spe=dim(data)[1]
Num_samp=dim(data)[2]
# 定义加权Jaccard距离函数
weighted_jaccard_distance <- function(XI, XJ) {
  min_values <- pmin(XI, XJ)
  max_values <- pmax(XI, XJ)
  weighted_dist <- sum(min_values) / sum(max_values)
  distance <- 1 - weighted_dist
  return(distance)
}

dist_matrix <- matrix(0, Num_spe, Num_spe) # 创建一个距离矩阵
for (i in 1:(Num_spe-1)) {
  for (j in (i+1):Num_spe) {
    XI <- as.numeric(as.vector(Gia_table[i,-1]))
    XJ <- as.numeric(as.vector(Gia_table[j,-1]))
    dist <- weighted_jaccard_distance(XI, XJ)
    dist_matrix[i, j] <- dist
    dist_matrix[j, i] <- dist
  }
}
dij <- as.matrix(dist_matrix)
#
TD <- numeric(Num_samp) # 创建一个长度为Num_samp的空数值向量TD
FD <- numeric(Num_samp) # 创建一个长度为Num_samp的空数值向量FD
FR <- numeric(Num_samp) # 创建一个长度为Num_samp的空数值向量FR
nFR <- numeric(Num_samp) # 创建一个长度为Num_samp的空数值向量FR
q=1 #这里 q=1进行计算

for (m in 1:Num_samp) { # 遍历每个样本
  otu_vector <- data[, m] # 获取第i个样本的物种丰度向量
  otu_matrix <- otu_vector %*% t(otu_vector) # 计算物种丰度向量的外积矩阵
  diag(otu_matrix) <- 0 # 将外积矩阵的对角线元素置为0
  TD[m] <- sum(otu_matrix^q) # 计算物种丰度的q次幂的和（TD）
  FD[m] <- sum(otu_matrix^q * dij) # 计算物种间距离乘以物种丰度的q次幂的和（FD）
  FR[m] <- TD[m] - FD[m] # 计算差值（FR）
  nFR[m] <- FR[m]/TD[m]
  print(m)
}


# 假设您的物种丰度表为 otu_table5，距离矩阵为 dij5
# genome-guild 对应表为 arc1645_cluster (需包含 genome 和 Cluster 两列)

# ---------------------------------------------------------
# 关键步骤 1：确保 Guild 信息与丰度表的行（Genome）完全一一对应
# ---------------------------------------------------------
# 提取 otu_table5 的行名 (即 Genome ID)
genome_ids <- rownames(data)

# 使用 match 函数将 arc1645_cluster 中的 Cluster 按照 otu_table5 的顺序对齐
# 结果是一个长度等于 Num_spe 的向量，里面按顺序存着每个物种对应的 Guild 编号
guild_vector <- arc1645_cluster$Cluster[match(genome_ids, arc1645_cluster$genome)]

# ---------------------------------------------------------
# 关键步骤 2：初始化结果向量
# ---------------------------------------------------------
TD <- numeric(Num_samp)      # 经典分类学多样性 (Simpson)
FD_base <- numeric(Num_samp) # 基础功能距离 (原版 FD / Rao's Q)
mGFD <- numeric(Num_samp)    # ★ 新增：Menhinick 调整后的 GFD ★
N_gen <- numeric(Num_samp)   # 保存每个样本的 N_genome，方便质控检查
S_gui <- numeric(Num_samp)   # 保存每个样本的 S_guild，方便质控检查

# ---------------------------------------------------------
# 关键步骤 3：遍历样本计算 mGFD
# ---------------------------------------------------------
for (m in 1:Num_samp) { 
  otu_vector <- data[, m] # 获取第 m 个样本的物种相对丰度向量
  
  # a. 筛选出该样本中实际存在（丰度 > 0）的物种索引
  present_idx <- which(otu_vector > 0)
  
  # b. 计算 N_genome（该样本中存在的基因组数量）
  N_genome <- length(present_idx)
  
  # 如果该样本里有物种存在，才进行计算（防止除以0的报错）
  if (N_genome > 0) {
    
    # c. 计算 S_guild（这些存在的物种，一共涵盖了几个 Guild？）
    present_guilds <- unique(guild_vector[present_idx])
    S_guild <- length(present_guilds)
    
    # d. ★ 计算 Menhinick 解耦惩罚系数 ★
    Menhinick_penalty <- S_guild / sqrt(N_genome)
    
    # e. 计算基础的矩阵和距离 (与您原代码一致)
    otu_matrix <- otu_vector %*% t(otu_vector) 
    diag(otu_matrix) <- 0 
    
    TD[m] <- sum(otu_matrix)               # Simpson 多样性
    FD_base[m] <- sum(otu_matrix * dij)   # 基础的 Rao's Q (分子后半部分)
    
    # f. ★ 计算最终的 mGFD ★
    mGFD[m] <- Menhinick_penalty * FD_base[m]
    
    # 记录参数方便后续作图和检查
    N_gen[m] <- N_genome
    S_gui[m] <- S_guild
    
  } else {
    # 如果样本为空，全部记为 0
    TD[m] <- 0
    FD_base[m] <- 0
    mGFD[m] <- 0
    N_gen[m] <- 0
    S_gui[m] <- 0
  }
  
  # 打印进度
  if (m %% 10 == 0) print(paste("Processed sample:", m, "/", Num_samp))
}

# ---------------------------------------------------------
# 关键步骤 4：将结果整合为 Dataframe 方便后续分析与绘图
# ---------------------------------------------------------
results_df <- data.frame(
  Sample = colnames(data),
  N_genome = N_gen,
  S_guild = S_gui,
  Menhinick_Coef = S_gui / sqrt(N_gen), # 避免除以0
  TD = TD,
  FD_base = FD_base,
  mGFD = mGFD
)

# 查看前几行结果
head(results_df)



