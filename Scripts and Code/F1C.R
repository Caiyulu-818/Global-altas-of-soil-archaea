#key_gene_heamap ####
key<-read.csv("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/02.arcfunction/arc_fun_ko/ko_p2_filled_key_t_info.csv",header = FALSE)
novel<-rep("a",1649)
key<-key%>%mutate(novel)
key$novel[key$V14=="s__"]<-"uSGBs"
key$novel[!key$V14=="s__"]<-"kSGBs"
table(key$novel)
key_c<-key[,c(2,10,21:107)]

key_c<-key_c[c(-1,-2,-4),]
colnames(key_c)<-key_c[1,]
colnames(key_c)[1]<-"genome"
colnames(key_c)[2]<-"class"
colnames(key_c)[89]<-"info"
key_c<-key_c[-1,]
write.csv(key_c,"/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/02.arcfunction/arc_fun_ko/key_c.csv")
key_c<-read.csv("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/02.arcfunction/arc_fun_ko/key_c.csv")

# 根据class分组，统计每列中不为0的genome数目占总genome的比例####
library(dplyr)
library(dplyr)
result_df2<- key_c %>%
  group_by(class) %>%
  reframe(
    # 对列3到87（即基因组列）进行操作
    across(3:87, ~ {
      total_non_zero <- sum(. != 0)  # 计算该列的非零基因组数量
      
      # 如果该列没有非零元素，直接返回0的占比
      if (total_non_zero == 0) {
        ksgb_count  <- 0
        usgb_count  <- 0
      } else {
        ksgb_count <- sum(. != 0 & info == "kSGBs")  # 计算ksgb在非零基因组中的数量
        usgb_count <- sum(. != 0 & info == "uSGBs")  # 计算usgb在非零基因组中的数量
        #ksgb_ratio <- ksgb_count / total_non_zero  # ksgb在非零基因组中的占比
        #usgb_ratio <- usgb_count / total_non_zero  # usgb在非零基因组中的占比
      }
      
      # 返回一个数据框，包含ksgb_ratio和usgb_ratio作为两个独立列
      c(
        ksgb_count = ksgb_count, 
        usgb_count = usgb_count
      )
    })
  )%>%mutate(
    info = rep(c("kSGB", "uSGB"), length.out = n())  # 交替生成 "usgb" 和 "ksgb"
  )
result_df2<-result_df2[,c(1,87,2:86)]


key_c_total_long<-gather(result_df1 ,key="gene",value = "Ratio of MAGs",-class)

key_c_total_long1<-merge(key_c_total_long,key_fun,by="gene",all.x = TRUE)
key_c_total_long1<-key_c_total_long1%>%filter(!fun=="SOX system")

#add metabolism  ####
fun_key_ko<-read.csv("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/02.arcfunction/arc_fun_ko/fun_ko_gene_r1.csv")
fun_key_ko1<-fun_key_ko%>%separate_rows(ko,sep = "/\\s*")
ntri_common<-c("narI, narV_ntri","nxrB","nxrA")
fun_key_ko1<-fun_key_ko1%>%filter(!X3%in%ntri_common)
fun_key_ko<-fun_key_ko[,c(1,2)]
fun_key_ko<-unique(fun_key_ko)
#key_c_total_long1$fun<-gsub("H4MPT Methyl-branch","WLP",key_c_total_long1$fun)

key_c_total_long2<-merge(key_c_total_long1,fun_key_ko,by="fun",all.x = TRUE)
table(key_c_total_long2$fun)
colnames(key_c_total_long2)[6]<-"Metabolism"
#合并相同行####
#key_c<-apply(key_c, 2, )
key_c<-key_c%>%mutate_at(vars(2:105), as.numeric)
key_c_total<-key_c%>%
  group_by(class)%>%
  summarise(across(.cols = everything(), .fns = sum))
 # mutate_at(vars(matches("^KO")), as.numeric)
key_c_total_long<-gather(key_c_total,key="Gene",value = "Gene_number",-class)
result_df4<-gather(result_df2,key="Gene",value = "Ratio",-c(class,info))
#画图 ####
table(key_c_total_long3$Ratio.of.MAGs)
colnames(key_c_total_long1)[5]<-"Gene"
ggplot(key_c_total_long1)
key_c_total_long2[is.na(key_c_total_long2$MAGs_discrete)]<-"0"
key_c_total_long3$MAGs_discrete <- cut(
  key_c_total_long3$`Ratio of MAGs`,
  breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
  labels = c("0", "<=0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", "0.8-1")
)

key_c_total_long3<-read.csv("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/02.arcfunction/arc_fun_ko/key_c_total_long3.csv")
key_c_total_long2<-read.csv("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/02.arcfunction/arc_fun_ko/key_c_total_long2_r1.csv")
# 重新排序Gene列的因子级别，根据Metabolism和fun进行排序
key_c_total_long2$fun<- reorder(
  key_c_total_long2$fun,
  key_c_total_long2$Metabolism # 使用每个分组的第一个值进行排序
)



table(key$V18)
key_c_total_long3<-key_c_total_long2%>%filter(!gene%in%ntri_common)
table(key$V18)
#gene_order
library(tidyverse)
library(ggplot2)
library(cowplot)

#delete fchA
key_c_total_long3$Ratio.of.MAGs<-key_c_total_long3%>%filter(!gene=="K01500")

key_c_total_long3<-read.csv("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/02.arcfunction/arc_fun_ko/key_c_total_long3.csv")

gene_order<-read.csv("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/02.arcfunction/gene_order.csv")
result_df6<-merge(result_df5,gene_order,by="gene",all.x=TRUE)
o<-gene_order$gene
o
result_df6$gene<-factor(result_df6$gene,levels = rev(o))
key_c_total_long4$gene<-factor(key_c_total_long4$gene,levels =rev(o))
table(key_c_total_long3$class)
result_df6$fun<-factor(result_df6$fun,levels = c(" Dicarboxylate-hydroxybutyrate cycle",
                                                               "3-Hydroxypropionate bi-cycle",
                                                               "Hydroxypropionate-hydroxybutylate cycle",
                                                               "CBB",
                                                               "rTCA",
                                                               "RGP +Serine=>Pyruvate",
                                                               "WLP",
                                                               "H4MPT Methyl-branch",
                                                               "Methanogenesis",
                                                               "Denitrification",
                                                               "Dissimilatory nitrate reduction",
                                                               "Nitrification",
                                                               "Assimilatory nitrate reduction",
                                                               "Nitrogen fixation",
                                                               "Dissimilatory sulfate reduction and oxidation",
                                                               "Assimilatory sulfate reduction",
                                                               "Sulfur reductase ",
                                                               "IAA pathway"))

key_c_total_long4$class<-gsub("c__","",key_c_total_long4$class)
result_df5$class<-factor(result_df5$class,levels = rev(c("Iainarchaeia",
"Micrarchaeia",
"B1Sed10-29",
"Nanoarchaeia",
"SpSt-1190",
"Nanosalinia",
"Aenigmatarchaeia",
"LC30",
"Lokiarchaeia",
"Thorarchaeia",
"Methanomethylicia",
"Thermoprotei",
"Thermoprotei_A",
"Bathyarchaeia",
"Nitrososphaeria",
"Nitrososphaeria_A",
"EX4484-205",
"Methanobacteria",
"Hydrothermarchaeia",
"Hadarchaeia",
"Thermococci",
"E2",
"SW-10-69-26",
"EX4484-6",
"Thermoplasmata",
"Themoplasmata_A",
"Methanonatronarchaeia",
"Halobacteria",
"Methanomicrobia",
"Syntropharchaeia",
"Bog-38",
"Methanosarcinia",
"Methanocellia")))
m<-table(key_c_total_long4$class)
m

arc1645_p<-arc1645[,c("phylum","class")]
arc1645_p<-unique(arc1645_p)
arc1645_p$class<-gsub("c__","",arc1645_p$class)
key_c_total_long4<- merge( key_c_total_long3,arc1645_p,by="class",all.x = TRUE)
library(cowplot)
library(ggplot2)
library(ggforce)
library(dplyr)
library(tidyr)
library(forcats)
library(cowplot)
result_df3<-gather(result_df2,key="gene",value = "count",fwdA..fmdA:DDC.)
result_df4<-spread(result_df3,key="info",value = "count",fill=0)
result_df5<-result_df4%>%filter(kSGB==0&!uSGB==0)
result_df5$gene<-gsub("\\.\\.", ", ",result_df5$gene)
result_df5$gene<-gsub("\\.$", "",result_df5$gene)

result_df5$class<-gsub("c__","",result_df5$class)


library(ggsci)
col
result_df5%>%filter(!is.na(gene))%>%ggplot(aes(x= fct_infreq(gene), y = fct_infreq(class))) +
  geom_point(aes(size=uSGB,color=fun),shape=20,stroke = 0.4 ) + 
  #scale_fill_manual(values = mycol)+
  scale_color_manual(values =mycol35)+
    #scale_color_gradient(low = "white", high = "black")+
  #scale_fill_gradient(low = "white", high = "black")+# 使用灰度颜色映射
  # 使用灰度颜色映射
 #scale_color_manual(values = c("<=0.25" = "grey70", "0.26-0.5" = "grey50", "0.51-0.75" = "grey30", "0.76-1" = "black") ) + 
 # scale_fill_manual(values = c("0"="white","<=0.25" = "grey80", "0.26-0.5" = "grey50", "0.51-0.75" = "grey30", "0.76-1" = "black") ) +
  scale_size_continuous(breaks = c(1, 5, 10, 15, 20, 27),limits = c(0.01, 27),name = "Gene only from uSGBs") +
  #theme_minimal() +
 # scale_fill_aaas(name = "") +
 # scale_shape_manual(name = "",values = c(21,16,16))+
  theme_cowplot(font_size = 7) +
  xlab("Gene") + 
  ylab("Class") +
 # coord_fixed(ratio = 0.8)+
  #coord_flip()+
 #facet_grid(.~fun,scales = "free",space = "free") +
 #facet_grid(fun~.,scales = "free",space = "free") +
  theme(#aspect.ratio = 0.7,
       strip.text.y  = element_text(angle = 0,face = "bold",),
       axis.title = element_text(face = "bold"),
        axis.text.x= element_text(angle = 45,
                                  hjust = 1),
        legend.position = "top")

ggsave("/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/02.arcfunction/arc_fun_ko/arc_1645_key_usgbonly_order.pdf",height=10,width = 12)
write.csv(key_c_total_long3,"/Users/caiyulu/Desktop/MAGcode/sediment/mag2.0/new_MAG/arc/02.arcfunction/arc_fun_ko/key_c_total_long3.csv")



