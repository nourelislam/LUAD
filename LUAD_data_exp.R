if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")


BiocManager::install("tidyverse")
install.packages("xlsx")
install.packages("gplots")
install.packages("ggplot2")
install.packages("RSQLite")

library("RSQLite")
library(ggplot2)
library(gplots)
library(dplyr)
library(tidyverse)
library(DESeq2)
library("reshape2") 
library(ComplexHeatmap)
library("pheatmap")
library("fission")
library(RColorBrewer)
library(xlsx)

head(df_new)

rownames(df_new) = df_new$X
# df_new = df_new[-1]
colnames(Similar_MI_RF_features) = "ensembl_id"
MI_RF_data = filter(df_new, X %in% Similar_MI_RF_features$ensembl_id)
rownames(MI_RF_data) = MI_RF_data$X
MI_RF_data = select(MI_RF_data, -c("X","symbol"))
five_common = filter(df_new, symbol %in% c("SFTPC", "EMP2","AGER","EPAS1","CAV1"))
rownames(five_common) = five_common$X
five_common = five_common[-1:-2]

metadata = data.frame(row.names = colnames(five_common), tumor_type = colnames(five_common))
metadata$tumor_type <- sub("PrimaryTumor", "T", metadata$tumor_type)
metadata$tumor_type <- sub("SolidTissueNormal", "N", metadata$tumor_type)

#########
df_new = df_new[-1:-2]
countData = df_new

all(colnames(countData) %in% rownames(metadata))
all(colnames(countData) == rownames(metadata))

##### Deseq2 whole diiferential gene expressoin ######
dds = DESeqDataSetFromMatrix( countData = countData , colData = metadata , design = ~ tumor_type)
dds.run = DESeq(dds)
#resultsNames(dds.run)
### direct results or specifying teh contrast (to make a res object based on two specific conditions/treatment)
res = results(dds.run, contrast = c("tumor_type", "T", "N"), alpha = 0.05, lfcThreshold = 1)
summary (res)
# remove nulls
res = res[complete.cases(res), ]
#summary(res.WT.33.cont)
res.df = as.data.frame(res)
plotMA(res, ylim=c(-1,1)) 

res.degs = res.df[res.df$padj< 0.05 & abs(res.df$log2FoldChange) >log2(2),]
write.xlsx(res.degs, file="DEGs.xlsx", sheetName="DEGs")
#write_excel_csv(res.degs, path = "visualization/res.degs.xlsx",quote=F)
countData$geneid = rownames(countData)
similar_all_genes_new = select(similar_all_genes_new, -c("_id", "_score"))
similar_all = filter(countData, geneid %in% similar_all_genes_new$query)
similar_all = left_join(similar_all, similar_all_genes_new, by= c("geneid" = "query"))
rownames(similar_all) = similar_all$symbol
similar_all = select(similar_all, -c( "symbol"))

common_all_genes = select(common_all_genes, -c("_id", "_score"))
common = filter(countData, geneid %in% common_all_genes$query)
common = left_join(common, common_all_genes, by= c("geneid" = "query"))
rownames(common) = common$symbol
common = select(common, -c("geneid", "symbol"))

X = arrange(res.degs, padj)
X = X[0:100,]
X = filter(countData, geneid %in% X$geneid)
X = X[-550]

x = filter(X, geneid %in% similar_all_genes_new$query)

### heatmap vis ####
colors = list(tumor_type = c(PrimaryTumor = "#1B9E77", SolidTissueNormal = "#D95F02"))
pheatmap(common, scale = "row",  cluster_rows = F,
         show_colnames = F,fontsize_row = 5,
         annotation_col = ann.col, annotation_colors = colors, rev(brewer.pal(4, "RdBu")))

pheatmap(common.Omar, scale = "row", cluster_rows = F, show_colnames = F,
         fontsize_row = 7, annotation_col = ann.col, annotation_colors = colors, 
         rev(brewer.pal(4, "RdBu")))

similar_all = as.matrix(similar_all)
X = as.matrix(X)
heatmap.2(X, scale = "row", col=colorRampPalette(c("blue","white","red")),
          trace ="none", main = "similar_all_genes")
################

degs_upto_2 = res.df[res.df$padj< 0.05 & abs(res.df$log2FoldChange) >log2(2) & abs(res.df$log2FoldChange) < log2(4),] ### abs function for absolute value "up and down" 
degs_upto_4 = res.df[res.df$padj< 0.05 & abs(res.df$log2FoldChange) >log2(4) & abs(res.df$log2FoldChange) < log2(16),]
degs_upto_6 = res.df[res.df$padj< 0.05 & abs(res.df$log2FoldChange) >log2(16) & abs(res.df$log2FoldChange) < log2(64),]
degs_upto_10 = res.df[res.df$padj< 0.05 & abs(res.df$log2FoldChange) >log2(64),]

##### log fold change of Omar's normalization ###
res.degs$geneid = rownames(res.degs)
log.similar_all_genes_new = filter(res.degs, geneid %in% similar_all_genes_new$query)
log.similar_rfe_RF_genes_new = filter(res.degs, geneid %in% similar_rfe_RF_genes_new$query)
log.similar_rfe_MI_genes_new = filter(res.degs, geneid %in% similar_rfe_MI_genes_new$query)
log.similar_MI_RF_genes_new = filter(res.degs, geneid %in% similar_MI_RF_genes_new$query)
log.common.Omar = filter(res.degs, geneid %in% common_all_genes$query)
log_Omar_sorted = arrange(log.common.Omar, log2FoldChange)
write_excel_csv(log_Omar_sorted, path = "Omar's normalization/log.common.Omar.sorted.xlsx",quote=F)

write.table(log.similar_MI_RF_genes_new, file="Omar's normalization/log.common", sep="\t", quote=F)
##### log fold change of Nour's normalization ###
log.similar_all_genes_Nour = filter(res.degs, geneid %in% similar_all_genes_new$query)
log.similar_rfe_RF_genes_Nour = filter(res.degs, geneid %in% similar_rfe_RF_genes_new$query)
log.similar_rfe_MI_genes_Nour = filter(res.degs, geneid %in% similar_rfe_MI_genes_new$query)
log.similar_MI_RF_genes_Nour = filter(res.degs, geneid %in% similar_MI_RF_genes_new$query)
log.common_Nour = filter(res.degs, geneid %in% common_all_genes$query)
log_Nour_sorted = arrange(log.common_Nour, log2FoldChange)

write_excel_csv(log_Nour_sorted, path = "Nour's normalization/log.common.Nour.sorted.xlsx",quote=F)
write.table(log.similar_MI_RF_genes_Nour, file="Nour's normalization/log.MI_RF", sep="\t", quote=F)

### log fold change of the original data #####
log.similar_all_genes_original = filter(res.degs, geneid %in% similar_all_genes_new$query)
log.similar_rfe_RF_genes_original = filter(res.degs, geneid %in% similar_rfe_RF_genes_new$query)
log.similar_rfe_MI_genes_original = filter(res.degs, geneid %in% similar_rfe_MI_genes_new$query)
log.similar_MI_RF_genes_original = filter(res.degs, geneid %in% similar_MI_RF_genes_new$query)
log.common.original = filter(res.degs, geneid %in% common_all_genes$query)

write.table(log.similar_MI_RF_genes_original, file="original data/log.MI_RF", sep="\t", quote=F)

write_excel_csv(log.common_Nour, path = "Nour's normalization/log.common.Nour.xlsx",quote=F)

## normalization to all data ###
dds <- DESeqDataSetFromMatrix(countData = countData, colData = metadata, design = ~ tumor_type)
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts = as.data.frame(normalized_counts)
normalized_counts$GeneID = rownames(normalized_counts)

common.Omar = filter(normalized_counts, GeneID %in% similar_all_genes_new$query)
common.Omar = left_join(common.Omar, similar_all_genes_new, by= c("GeneID" = "query"))
rownames(common.Omar) = common.Omar$symbol
common.Omar = select(common.Omar, -c( "GeneID" ,"symbol"))
transposed = t(common.Omar)
transposed = as.data.frame(transposed)
transposed$type = rownames(transposed)
metadata$sample_type = rownames(metadata)
joining = left_join(transposed, metadata, by= c("type" = "sample_type")) 
joining = joining[-13]
vis_data <- melt(joining, id = "tumor_type")
names(vis_data)[2] = "Genes"

### Box plots ###
ggplot(vis_data, aes(x = tumor_type, y = value, fill = Genes)) +  # ggplot function
  geom_boxplot(outlier.shape = NA)+
  labs(x = "Sample_type", y = "Normalized_counts")+
  theme(legend.position = "none")+
  facet_wrap(~Genes, scale="free")+
  geom_jitter(shape=20, position=position_jitter(0.3))
### MI_RF Vis ##
similar_MI_RF_genes_new = select(similar_MI_RF_genes_new, -c("_id", "_score"))
MI_RF = filter(normalized_counts, GeneID %in% similar_MI_RF_genes_new$query)
MI_RF = left_join(MI_RF, similar_MI_RF_genes_new, by= c("GeneID" = "query"))
rownames(MI_RF) = MI_RF$symbol
MI_RF = select(MI_RF, -c( "GeneID" ,"symbol"))
MI_RF = MI_RF[-c(2,4,6,7,8,9,10,11,12,13,16),]
transposed = t(MI_RF)
transposed = as.data.frame(transposed)
transposed$type = rownames(transposed)
metadata$sample_type = rownames(metadata)
joining = left_join(transposed, metadata, by= c("type" = "sample_type")) 
joining = joining[-8]
vis_data <- melt(joining, id = "tumor_type")
names(vis_data)[2] = "Genes"

ggplot(vis_data, aes(x = tumor_type, y = value, fill = Genes)) +  # ggplot function
  geom_boxplot(outlier.shape = NA)+
  labs(x = "Sample_type", y = "Normalized_counts")+
  theme(legend.position = "none")+
  facet_wrap(~Genes, scale="free")+
  geom_jitter(shape=18, position=position_jitter(0.3))

### rfe_MI Vis ###
similar_rfe_MI_genes_new = select(similar_rfe_MI_genes_new, -c("_id", "_score"))
rfe_MI = filter(normalized_counts, GeneID %in% similar_rfe_MI_genes_new$query)
rfe_MI = left_join(rfe_MI, similar_rfe_MI_genes_new, by= c("GeneID" = "query"))
rownames(rfe_MI) = rfe_MI$symbol
rfe_MI = select(rfe_MI, -c( "GeneID" ,"symbol"))
transposed = t(rfe_MI)
transposed = as.data.frame(transposed)
transposed$type = rownames(transposed)
metadata$sample_type = rownames(metadata)
joining = left_join(transposed, metadata, by= c("type" = "sample_type")) 
joining = joining[-13]
vis_data <- melt(joining, id = "tumor_type")
names(vis_data)[2] = "Genes"

ggplot(vis_data, aes(x = tumor_type, y = value, fill = Genes)) +  # ggplot function
  geom_boxplot(outlier.shape = NA)+
  labs(x = "Sample_type", y = "Normalized_counts")+
  theme(legend.position = "none")+
  facet_wrap(~Genes, scale="free")+
  geom_jitter(shape=18, position=position_jitter(0.3))

### rfe_RF Vis ###
similar_rfe_RF_genes_new = select(similar_rfe_RF_genes_new, -c("_id", "_score"))
rfe_RF = filter(normalized_counts, GeneID %in% similar_rfe_RF_genes_new$query)
rfe_RF = left_join(rfe_RF, similar_rfe_RF_genes_new, by= c("GeneID" = "query"))
rownames(rfe_RF) = rfe_RF$symbol
rfe_RF = select(rfe_RF, -c( "GeneID" ,"symbol"))
rfe_RF = rfe_RF[-c(7,11,18,26,29,30,31,33,34,47,40,42),]
transposed = t(rfe_RF)
transposed = as.data.frame(transposed)
transposed$type = rownames(transposed)
metadata$sample_type = rownames(metadata)
joining = left_join(transposed, metadata, by= c("type" = "sample_type")) 
joining = joining[-39]
#last_25 = joining[,25:51]
vis_data <- melt(joining, id = "tumor_type")
names(vis_data)[2] = "Genes"

ggplot(vis_data, aes(x = tumor_type, y = value, fill = Genes)) +  # ggplot function
  geom_boxplot(outlier.shape = NA)+
  labs(x = "Sample_type", y = "Normalized_counts")+
  theme(legend.position = "none")+
  facet_wrap(~Genes, scale="free")+
  geom_jitter(shape=18, position=position_jitter(0.3))

###################

joining = left_join(transposed, metadata, by= c("type" = "sample_type")) 
joining = joining[-51]
last_25 = select(last_25, -c("tumor_type"))
joining = select(joining, -c(colnames(last_25)))
vis_data <- melt(joining, id = "tumor_type")
names(vis_data)[2] = "Genes"

ggplot(vis_data, aes(x = tumor_type, y = value, fill = Genes)) +  # ggplot function
  geom_boxplot(outlier.shape = NA)+
  labs(title="RFE_RF common genes",
       x = "Sample_type", y = "Normalized_counts")+
  facet_wrap(~Genes, scale="free")+
  geom_jitter(shape=18, position=position_jitter(0.3))


countData$name = rownames(countData)
vis.five = df_new[-1]
vis.five = vis.five[-1]
ann.col = metadata[-2]
colnames(ann.col) = "Sample_type"
normalized_counts = normalized_counts[-2,]

metadata = data.frame(row.names = colnames(five_RFE_RF), tumor_type = colnames(five_RFE_RF))
metadata$tumor_type <- gsub("[0-9]", "", metadata$tumor_type)
metadata$tumor_type <- gsub("[:.:]", "", metadata$tumor_type)

x <- res.degs %>% 
  rownames_to_column(var="samplename") %>% 
  as_tibble()



















