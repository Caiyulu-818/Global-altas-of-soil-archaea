#dbrda####
# 加载R包 Load the package
suppressWarnings(suppressMessages(library(vegan)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(ggpubr)))
suppressWarnings(suppressMessages(library(ggrepel)))
suppressWarnings(suppressMessages(library(rdacca.hp)))
#install.packages("rdacca.hp")
library(adespatial)
library(geosphere)
require(ade4)
require(SoDA)
library(tidyverse)
setwd("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD")
paddy<-read.csv("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/02.arcfunction/arc1645/sample.paddy.soil.csv")
paddy<- paddy %>%
  mutate(Ecosystem = ifelse(sample == "SRR13926678", "Paddy soil", Ecosystem))
# 基于CRAN安装R包，检测没有则安装
#p_list =c("vegan","ggplot2","ggpubr","ggrepel","rdacca.hp","vegan","psych","reshape2")
#for(p in p_list){if(!requireNamespace(p)){install.packages(p)}
 # library(p, character.only =TRUE, quietly =TRUE, warn.conflicts =FALSE)}

## dbRDA plot using R softwaredbRDA图R语言实战
# load data
# 读入物种数据，以细菌 OTU 水平丰度表为例
otu = readRDS("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/otu.rds")

#otu <- data.frame(t(arc.mag.abundance2_long_wide)) %>% na.omit()
#write_rds(otu,"/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FD/otu.rds")
# 分组数据
# group data
#matadata <- read.table(paste("data/group_data.txt",sep=""), header=T, row.names=1, sep="\t", comment.char="")
#otu = otu[rownames(otu) %in% rownames(matadata), ]

# 读取环境数据
# confounding factors
env.info.cluster<-readRDS("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F4/env.info.cluster.rds")

env = env.info.cluster
env = na.omit(env)
rownames(env)<-NULL
env = env[env$site %in% rownames(otu), ]
env <- env[!duplicated(env$site), ] 
rownames(env) =env$site
rownames = rownames(env)%>%as.data.frame()
otu = otu%>%filter(rownames(otu)%in%rownames(env))

# 根据原理一步步计算 db-RDA
# Calculate db-RDA step by step according to the principle
# 计算样方距离，以 Bray-curtis 距离为例，详情 ?vegdist
# Calculate the sample distance, taking Bray-curtis distance as an example, details ?vegdist
dis_bray <- vegdist(otu, method = 'bray')

# PCoA 排序，这里通过 add = TRUE校正负特征值，详情 ?cmdscale
# PCoA sorting, here add = TRUE is used to correct negative eigenvalues, details ?cmdscale
pcoa <- cmdscale(dis_bray, k = nrow(otu) - 1, eig = TRUE, add = TRUE)

# 提取 PCoA 样方得分（坐标）
# Extract PCoA sample scores (coordinates)
pcoa_site <- pcoa$point

# db-RDA，环境变量与 PCoA 轴的多元回归
# db-RDA, multiple regression of environmental variables and PCoA axes
# 通过 vegan 包的 RDA 函数 rda() 执行，详情 ?rda
# Execute via the RDA function rda() of the vegan package, details ?rda
# 排序行名
pcoa_site <- pcoa_site[order(rownames(pcoa_site)), ]
env <- env[order(rownames(env)), ]

###经纬度fnction####加入准备空间变量
env_geo<-env[,c(15,14)]###经纬度###邻域矩阵主坐标 (PCNM) 85、86。MEM分析根据 GFS 的地理坐标生成一组正交变量85 ；这些变量在 db-RDA 中作为解释变量来模拟社区数据中的空间结构
pcnmenv_geo<-pcnm(distm(env_geo,fun=distGeo))$vectors
pcnmenv_geo<-as.data.frame(pcnmenv_geo)
db_rda_geo <- rda(pcoa_site ~ ., data = pcnmenv_geo)
#anova(db_rda_geo, by = "terms")  # 检查显著性
env<-data.frame(env[,c(1:13)],pcnmenv_geo)#所有env
write.csv(env,"env_pcnm.csv")
#开始rda 分析
db_rda <- rda(pcoa_site, env, scale = FALSE)
# summary(db_rda)
anova.cca(db_rda)#用置换检验（Permutation test）检验整个RDA模型是否显著0.001 
# Step 7: forward.sel 筛选重要变量
p.R2a <- RsquareAdj(db_rda)$adj.r.squared#0.1784489
print(p.R2a)
p.env.fwd <- forward.sel(pcoa_site, env, adjR2thresh=p.R2a,
                         nperm=999)#提取 调整后的 R²，避免过拟合影响模型解释率
write.csv(p.env.fwd,"p.env.fwd.csv")
p.env.fwd<-read.csv("p.env.fwd.csv",row.names = 1)
env.sign_p <- sort(p.env.fwd$order)
env.sign.tax1<- env[,c(env.sign_p)]###和群落显著相关环境因子 done
# Step 8: 构建简化 db-RDA 模型
db_rda_final1 <- rda(pcoa_site ~ ., data =env.sign.tax1, scale = FALSE)

write_rds(db_rda_final1,file = "db_rda_finalall_tax.rds")
db_rda_final1<-read_rds("db_rda_finalall_tax.rds")
# Step 9: 模型统计与解释度
anova.cca(db_rda_final1)#0.001
RsquareAdj(db_rda_final1)#$r.squared 0.3332969 $adj.r.squared[1] 0.2651817

# Step 10: envfit 被动拟合环境因子投影（用于箭头图）
envfit_tax1 <- envfit(db_rda_final1, env.sign.tax1, permutations = 999)
r <- as.matrix(envfit_tax1$vectors$r)
p <- as.matrix(envfit_tax1$vectors$pvals)
env.p <- cbind(r,p)
colnames(env.p) <- c("r2","p-value")
fittax1 <- as.data.frame(env.p)
fittax1$p.adj = p.adjust(fittax1$`p-value`, method = 'BH')
# KK
write.csv(as.data.frame(fittax1),file="./rdaenvfit_magtaxaall.csv")
fittax1<-read.csv("./rdaenvfit_magtaxaall.csv",row.names=1)
fittax1_filter<-fittax1%>%filter(p.adj<0.05)
env.sign.tax1_filter<-env.sign.tax1[,colnames(env.sign.tax1)%in%rownames(fittax1_filter)]
# Step 11:利用rdacca.hp计算每个环境因子的效应
# Use rdacca.hp to calculate the effect of each environmental factor
#dis_bray <- as.dist(distance_mat)
#相对解释度
cap.hp = rdacca.hp(dis_bray, env.sign.tax1_filter, method = 'dbRDA', type = 'R2', scale = FALSE)

#step 12 写出结果
# 1. 准备 rdacca.hp 的结果
expl_var <- cap.hp$Hier.part%>%as.data.frame()
colnames(expl_var) <- c("Unique", "Average.share", "Individual", "I.perc")
expl_var$Absolute_explained_variation <- expl_var$I.perc * cap.hp$Total_explained_variation
expl_var$Total_explained_variation<-cap.hp$Total_explained_variation
# 2. 将变量名列作为显式列
expl_var$Variable <- rownames(expl_var)
fittax1_filter$Variable <- rownames(fittax1_filter)

# 3. 合并两个表（按 Variable 合并）
combined_result <- merge(fittax1_filter, expl_var,
                         by = "Variable", all.x = TRUE)

# 4. 排序可选（按解释度或 p 值）
combined_result <- combined_result[order(combined_result$p.adj, decreasing = FALSE), ]

# 5. 输出结果
write.csv(combined_result, "combined_explained_variation_taxall.csv", row.names = FALSE)

print("done")