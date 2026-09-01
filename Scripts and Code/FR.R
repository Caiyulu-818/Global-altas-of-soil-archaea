library(dplyr)
library(tidyr)
library(readxl)

Gia <- read_xlsx("./gia.xlsx", sheet = 2)
Gia <- Gia[rowSums(Gia[,-1])!=0, ]
saveRDS(arc.mag.abundance,"arc.mag.abundance.rds")
saveRDS(arc.ko,"arc.ko.rds")
arc.mag.abundance<-readRDS("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F4/arc.mag.abundance.rds")
data <-arc.mag.abundance
rownames(data)<-NULL
data<-data%>%column_to_rownames(var="Genome")
data <- sapply(data, as.numeric)
data <- data[, colSums(data > 0 | is.na(data), na.rm = TRUE) != nrow(data)]
data<- apply(data, 2, function(x)x/sum(x))#归一化处理，每个样本中物种丰度的总和为1
rownames(data) <- arc.mag.abundance$Genome

data1<-as.data.frame(t(cluster1_tax))
data1 <- sapply(data1, as.numeric)
data1 <- data1[, colSums(data1 > 0 | is.na(data1), na.rm = TRUE) != nrow(data1)]
otu_table1 <- apply(data1, 2, function(x)x/sum(x))#归一化处理，每个样本中物种丰度的总和为1
rownames(otu_table1)<-colnames(cluster1_tax)

data2<-as.data.frame(t(cluster2_tax))
data2 <- sapply(data2, as.numeric)
data2 <- data2[, colSums(data2 > 0 | is.na(data2), na.rm = TRUE) != nrow(data2)]
otu_table2 <- apply(data2, 2, function(x)x/sum(x))#归一化处理，每个样本中物种丰度的总和为1
rownames(otu_table2)<-colnames(cluster2_tax)

data3<-as.data.frame(t(cluster3_tax))
data3 <- sapply(data3, as.numeric)
data3 <- data3[, colSums(data3 > 0 | is.na(data3), na.rm = TRUE) != nrow(data3)]
otu_table3 <- apply(data3, 2, function(x)x/sum(x))#归一化处理，每个样本中物种丰度的总和为1
rownames(otu_table3)<-colnames(cluster3_tax)

data4<-as.data.frame(t(cluster4_tax))
data4 <- sapply(data4, as.numeric)
data4 <- data4[, colSums(data4 > 0 | is.na(data4), na.rm = TRUE) != nrow(data4)]
otu_table4 <- apply(data4, 2, function(x)x/sum(x))#归一化处理，每个样本中物种丰度的总和为1
rownames(otu_table4)<-colnames(cluster4_tax)

data5<-as.data.frame(t(cluster5_tax))
data5 <- sapply(data5, as.numeric)
data5 <- data5[, colSums(data5 > 0 | is.na(data5), na.rm = TRUE) != nrow(data5)]
otu_table5 <- apply(data5, 2, function(x)x/sum(x))#归一化处理，每个样本中物种丰度的总和为1
rownames(otu_table5)<-colnames(cluster5_tax)

#ko matrix
Gia1<-cluster1_genome%>%column_to_rownames(var="genome") #Genes
# 删除所有行和为零或全NA的行
Gia1 <- Gia1[rowSums(Gia1 > 0 | is.na(Gia1), na.rm = TRUE) > 0, ]
# 删除所有列和为零或全NA的列
Gia1 <- Gia1[, colSums(Gia1 > 0 | is.na(Gia1), na.rm = TRUE) > 0]

#ko matrix
Gia2<-cluster2_genome%>%column_to_rownames(var="genome") #Genes
# 删除所有行和为零或全NA的行
Gia2 <- Gia2[rowSums(Gia2 > 0 | is.na(Gia2), na.rm = TRUE) > 0, ]
# 删除所有列和为零或全NA的列
Gia2 <- Gia2[, colSums(Gia2 > 0 | is.na(Gia2), na.rm = TRUE) > 0]


Gia3<-cluster3_genome%>%column_to_rownames(var="genome") #Genes
# 删除所有行和为零或全NA的行
Gia3 <- Gia3[rowSums(Gia3 > 0 | is.na(Gia3), na.rm = TRUE) > 0, ]
# 删除所有列和为零或全NA的列
Gia3 <- Gia3[, colSums(Gia3 > 0 | is.na(Gia3), na.rm = TRUE) > 0]

Gia4<-cluster4_genome%>%column_to_rownames(var="genome") #Genes
# 删除所有行和为零或全NA的行
Gia4 <- Gia4[rowSums(Gia4 > 0 | is.na(Gia4), na.rm = TRUE) > 0, ]
# 删除所有列和为零或全NA的列
Gia4 <- Gia4[, colSums(Gia4 > 0 | is.na(Gia4), na.rm = TRUE) > 0]

Gia5<-cluster5_genome%>%column_to_rownames(var="genome") #Genes
# 删除所有行和为零或全NA的行
Gia5 <- Gia5[rowSums(Gia5 > 0 | is.na(Gia5), na.rm = TRUE) > 0, ]
# 删除所有列和为零或全NA的列
Gia5 <- Gia5[, colSums(Gia5 > 0 | is.na(Gia5), na.rm = TRUE) > 0]


#对应信息
#table <- merge(Gia,data,by="genus")  #为了使gia和otu的genus名字对应起来
#table_fr<-merge(arc.ko,data,by.x="genome",by.y="Genome")
Gia_table <- arc.ko%>%column_to_rownames(var="genome") #Genes
# 删除所有行和为零或全NA的行
Gia_table <- Gia_table[rowSums(Gia_table > 0 | is.na(Gia_table), na.rm = TRUE) > 0, ]
# 删除所有列和为零或全NA的列
Gia_table <- Gia_table[, colSums(Gia_table > 0 | is.na(Gia_table), na.rm = TRUE) > 0]
otu_table <- data
class(otu_table)
otu_table <- sapply(otu_table, as.numeric)
otu_table <- otu_table[, colSums(otu_table > 0 | is.na(otu_table), na.rm = TRUE) != nrow(otu_table)]
otu_table <- apply(otu_table, 2, function(x)x/sum(x))#归一化处理，每个样本中物种丰度的总和为1
rownames(otu_table) <- arc.mag.abundance$Genome


Num_spe=dim(otu_table5)[1]
Num_samp=dim(otu_table5)[2]
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
    XI <- as.numeric(as.vector(Gia5[i,-1]))
    XJ <- as.numeric(as.vector(Gia5[j,-1]))
    dist <- weighted_jaccard_distance(XI, XJ)
    dist_matrix[i, j] <- dist
    dist_matrix[j, i] <- dist
  }
}
dij5 <- as.matrix(dist_matrix)
#
TD <- numeric(Num_samp) # 创建一个长度为Num_samp的空数值向量TD
FD <- numeric(Num_samp) # 创建一个长度为Num_samp的空数值向量FD
FR <- numeric(Num_samp) # 创建一个长度为Num_samp的空数值向量FR
nFR <- numeric(Num_samp) # 创建一个长度为Num_samp的空数值向量FR
q=1 #这里 q=1进行计算

for (m in 1:Num_samp) { # 遍历每个样本
  otu_vector <- otu_table5[, m] # 获取第i个样本的物种丰度向量
  otu_matrix <- otu_vector %*% t(otu_vector) # 计算物种丰度向量的外积矩阵
  diag(otu_matrix) <- 0 # 将外积矩阵的对角线元素置为0
  TD[m] <- sum(otu_matrix^q) # 计算物种丰度的q次幂的和（TD）
  FD[m] <- sum(otu_matrix^q * dij5) # 计算物种间距离乘以物种丰度的q次幂的和（FD）
  FR[m] <- TD[m] - FD[m] # 计算差值（FR）
  nFR[m] <- FR[m]/TD[m]
  print(m)
}
sample= colnames(otu_table5)
fr5 <- data.frame(sample=sample,TD=TD,FD=FD,FR=FR,nFR=nFR)

fr<-fr%>%mutate(Group="Total")
fr1<-fr1%>%mutate(Group="Cluster 1")
fr2<-fr2%>%mutate(Group="Cluster 2")
fr3<-fr3%>%mutate(Group="Cluster 3")
fr4<-fr4%>%mutate(Group="Cluster 4")
fr5<-fr5%>%mutate(Group="Cluster 5")

allfr<-rbind(fr,fr1,fr2,fr3,fr4,fr5)
library(dplyr)

summary_stats <- na.omit(allfr) %>%
  group_by(Group) %>%
  summarize(
    Median = median(nFR),  # 请将 nFR 替换为您要计算的数值列名
    Mean = mean(nFR),     # 请将 nFR 替换为您要计算的数值列名
    .groups = "drop"
  ) %>%
  arrange(desc(Group))

# 打印结果
print(summary_stats)
write.csv(summary_stats,"./allfrsummary_stats.csv",row.names = F)