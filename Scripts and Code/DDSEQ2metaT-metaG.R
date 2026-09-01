############################################################
## Function group differential analysis
## Metagenome + Metatranscriptome
############################################################


############################################################
## 1. Packages
############################################################

library(tidyverse)
library(DESeq2)
library(clusterProfiler)
library(UpSetR)
library(Rtsne)
library(vegan)
library(pheatmap)
library(ComplexHeatmap)
library(data.table)
library(cowplot)
library(ggplot2)



############################################################
## Part I. Metagenome functional group differential analysis
############################################################


setwd(
"/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/02.arcfunction/arc1645/"
)


############################################################
## 1. Prepare KO annotation
############################################################


kr <- read.delim(
"KOannotation.tsv",
stringsAsFactors = FALSE
)


## split k141 annotation

k141 <- kr[grepl("^k141",kr[,1]),]

k141$contig <- k141$gene

k141$gene <-
gsub(
"(k141)(.*bin\\.\\d+)",
"\\2\\1",
k141$gene
)


no141 <- kr[!grepl("^k141",kr[,1]),]

no141$contig <- no141$gene

no141$gene <-
gsub(
"_genomic",
"_genomick141",
no141$gene
)


krbind <- rbind(
k141,
no141
)


sep_krbind <-
separate(
krbind,
col="gene",
into=c("genome",NA),
sep="k"
)


sep_krbind$genome <-
gsub(
"_genomic",
"",
sep_krbind$genome
)


sep_krbind_sum <-
sep_krbind %>%
group_by(genome,contig) %>%
summarise(
count=n(),
.groups="drop"
)


write.csv(
sep_krbind_sum,
"sep_krbind_sumarc.csv"
)



############################################################
## 2. Construct KO abundance matrix
############################################################


final_gathered_ko_info <-
read.csv(
"/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/Rscript&sheet/finalgathered3_fun.cluster.csv"
)


ko_wide <-
final_gathered_ko_info %>%
select(genome,KO,count) %>%
distinct() %>%
spread(
key=genome,
value=count,
fill=0
) %>%
column_to_rownames("KO")


ko_wide <-
ko_wide[rowSums(ko_wide!=0)>0,]



group_ko <-
final_gathered_ko_info %>%
select(genome,Cluster) %>%
unique()



############################################################
## 3. DESeq2 differential KO analysis
############################################################


dds <-
DESeqDataSetFromMatrix(
countData=ko_wide,
colData=group_ko,
design=~Cluster
)


dds <- estimateSizeFactors(
dds,
type="poscounts"
)


dds <-
DESeq(
dds,
fitType="local"
)


dds$Cluster <-
as.factor(dds$Cluster)



############################################################
## 4. Pairwise cluster comparison
############################################################


cluster_levels <-
c("1","2","3","4","5")


comparison_list <-
combn(
cluster_levels,
2,
simplify=FALSE
)


setwd(
"/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/01.arcinfo/F3/FC"
)


for(comp in comparison_list){

contrast_name <-
paste0(
"Cluster_",
comp[2],
"_vs_",
comp[1]
)


res <-
results(
dds,
contrast=c(
"Cluster",
comp[2],
comp[1]
)
)%>%
as.data.frame()


res_filtered <-
res[
!is.na(res$padj) &
res$padj<0.05 &
abs(res$log2FoldChange)>1,
]


write.csv(
res_filtered,
paste0(
contrast_name,
".csv"
)
)

}



############################################################
## 5. Extract cluster-specific enriched KO
############################################################


get_enriched_genes <- function(
cluster_id,
all_res_list,
direction="up",
lfc_cutoff=1
){

cluster_id <- as.character(cluster_id)

others <-
setdiff(
as.character(1:5),
cluster_id
)

gene_sets <- list()


for(target_id in others){

forward <-
paste0(
"res",
cluster_id,
"_",
target_id
)


reverse <-
paste0(
"res",
target_id,
"_",
cluster_id
)



if(forward %in% names(all_res_list)){

res <- all_res_list[[forward]]

genes <-
rownames(
subset(
res,
padj<0.05 &
log2FoldChange>lfc_cutoff
)
)


}else{


res <- all_res_list[[reverse]]

genes <-
rownames(
subset(
res,
padj<0.05 &
log2FoldChange< -lfc_cutoff
)
)

}


gene_sets[[target_id]] <- genes

}


Reduce(
intersect,
gene_sets
)

}



############################################################
## 6. KEGG enrichment
############################################################


up_cluster1 <-
get_enriched_genes(
"1",
all_res_list
)


up_cluster2 <-
get_enriched_genes(
"2",
all_res_list
)


up_cluster3 <-
get_enriched_genes(
"3",
all_res_list
)


up_cluster4 <-
get_enriched_genes(
"4",
all_res_list
)


up_cluster5 <-
get_enriched_genes(
"5",
all_res_list
)



enrichment_results5 <-
enricher(
up_cluster5,
TERM2GENE=kegg_level[c(7,7)],
TERM2NAME=kegg_level[c(7,6)],
minGSSize=1
)



############################################################
## 7. KO overlap analysis (UpSet)
############################################################


df_wide <-
final_gathered_ko_info %>%
select(KO,Cluster,count) %>%
group_by(KO,Cluster) %>%
summarise(
count=sum(count),
.groups="drop"
) %>%
spread(
Cluster,
count,
fill=0
)


rownames(df_wide)<-
df_wide$KO

df_wide$KO<-NULL


df_binary<-df_wide

df_binary[df_binary>0]<-1



############################################################
## 8. KO t-SNE
############################################################


df_features_norud <-
df_wide[
!duplicated(df_wide),
]


tsne_result <-
Rtsne(
as.matrix(df_features_norud),
perplexity=500
)


############################################################
## 9. Functional pathway classification
############################################################


## Carbon fixation
Carbon_fix <- c(...)

## Methane
Methane <- c(...)

## Carbon decomposition
Carbon_decomp <- c(...)

## Nitrogen
Nitrogen <- c(...)

## Sulfur
Sulfur <- c(...)

## Phosphorus
Phosphorus <- c(...)

## IAA
IAA_pathway <- c(...)

## Pollutant
Pollutant <- c(...)


pathway_class <-
c(
Carbon_fix,
Methane,
Carbon_decomp,
Nitrogen,
Sulfur,
Phosphorus,
IAA_pathway,
Pollutant
)



############################################################
## 10. MAG functional heatmap
############################################################


mat_meta <-
sep_krbind_metabolism2 %>%
select(Cluster,fun,count) %>%
pivot_wider(
names_from=fun,
values_from=count,
values_fill=0
)%>%
column_to_rownames("Cluster")%>%
as.matrix()



pheatmap(
mat_meta,
cluster_rows=FALSE,
cluster_cols=FALSE,
scale="row"
)




############################################################
## Part II. Metatranscriptome functional analysis
############################################################



############################################################
## 1. Read transcriptome abundance
############################################################


bind5 <-
fread(
"/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/transcriptomic/trans.gene.quant_TAX.csv"
)


bind5 <-
bind5 %>%
select(-matches("MT"))



############################################################
## 2. Extract metabolic genes
############################################################


sep_trans_metabolism <-
bind5 %>%
filter(
KO %in% Metabolism_gene1$KO
)



expr_orf_wide <-
bind5[,c(1,4:1278,1281)]



############################################################
## 3. Convert TPM to long format
############################################################


setDT(expr_orf_wide)


sample_cols <-
setdiff(
names(expr_orf_wide),
c("Name","genome")
)


expr_long <-
melt(
expr_orf_wide,
id.vars=c("Name","genome"),
measure.vars=sample_cols,
variable.name="sample",
value.name="TPM"
)


expr_long[is.na(TPM),TPM:=0]



############################################################
## 4. Add cluster information
############################################################


arc1645 <-
read.csv(
"arc1645.clusterinfo.csv"
)


colnames(arc1645)[2]<-"genome"



expr_long <-
expr_long[
arc1645[,c(2,3)],
on="genome"
]



############################################################
## 5. Relative expression calculation
############################################################


cluster_median <-
expr_long[
,
.(median_TPM=median(TPM[TPM>0],na.rm=TRUE)),
by=.(Cluster,sample)
]


expr_orf_relative <-
cluster_median[
expr_long,
on=.(Cluster,sample)
]


expr_orf_relative[
median_TPM>0,
RelativeExpr:=TPM/median_TPM
]



############################################################
## 6. KO level aggregation
############################################################


expr_orf_global <-
expr_orf_relative %>%
group_by(Name,genome,Cluster)%>%
summarise(
mean_relative=mean(RelativeExpr,na.rm=TRUE),
.groups="drop"
)%>%
left_join(
koannotation.arc[,c(1,2)],
by=c("Name"="gene")
)%>%
group_by(KO,Cluster)%>%
summarise(
mean_relative=sum(mean_relative,na.rm=TRUE),
.groups="drop"
)



############################################################
## 7. Transcript functional heatmap
############################################################


mat_trans <-
expr_orf_global %>%
select(Cluster,fun,mean_relative)%>%
pivot_wider(
names_from=fun,
values_from=mean_relative,
values_fill=0
)%>%
column_to_rownames("Cluster")%>%
as.matrix()



pheatmap(
mat_trans,
cluster_rows=FALSE,
cluster_cols=FALSE,
scale="row"
)
