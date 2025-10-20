rm(list = ls())

setwd("/data/nas1/yuanyt/project/05.YQQL-0602-6/")
if (! dir.exists("./12_GSEA")){
  dir.create("./12_GSEA")
}
setwd("./12_GSEA")


library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(limma)
library(tidyverse)
library(lance)

## 01 Read data----------
kegggmt <- read.gmt("../../../pipeline/GSEA/c2.cp.kegg.v2023.1.Hs.symbols.gmt") #Read the gmt file
data <- read.csv("../00_Rowdata/tumor.tcga.csv",header=TRUE, row.names=1, check.names=FALSE)
riskScore <- readRDS("../08_LASSO/riskScore.rds")
riskScore$sample <- rownames(riskScore)
riskScore$risk <- ifelse(riskScore$riskScore > median(riskScore$riskScore),"High Risk","Low Risk")



dat_final <- data[,riskScore$sample]

set.seed(1)
riskScore$risk <- factor(riskScore$risk,levels = c("High Risk","Low Risk"))
colData<-data.frame(riskScore) 
dds<-DESeqDataSetFromMatrix(countData = round(dat_final),colData=colData,design = ~risk)
dds = dds[rownames(counts(dds)) > 1,]
dds<-DESeq(dds)
dds

res <- results(dds, contrast = c("risk","High Risk","Low Risk"))
resOrdered <- res[order(res$padj),]
DEG <- as.data.frame(resOrdered) %>% na.omit() 
saveRDS(DEG,file='DEG.rds')
DEG <- readRDS("DEG.rds")

risk_gene_symbol <- rownames(DEG)
risk_gene <- DEG  

gene <- bitr(risk_gene_symbol,
             fromType = "SYMBOL",
             toType = c("ENTREZID"),
             OrgDb = "org.Hs.eg.db")  
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)  
colnames(gene) <- c("GeneSymbol","ENTREZID")  
risk_gene$GeneSymbol <- rownames(risk_gene)  
data_all <- risk_gene %>% dplyr::inner_join(.,gene,by="GeneSymbol")    
data_all_sort <- data_all %>% arrange(desc(log2FoldChange))  
genelist <- data_all_sort$log2FoldChange #Extract the foldchanges from largest to smallest
names(genelist) <- data_all_sort$GeneSymbol #Add the corresponding ENTREZID to the foldchange extracted above
write.csv(genelist,file="genelist.csv")


###KEGG enrichment
KEGG<-GSEA(genelist,TERM2GENE = kegggmt,pvalueCutoff = 0.05) #GSEA Analysis
sortKEGG<-data.frame(KEGG)
sortKEGG<-sortKEGG[order(sortKEGG$p.adjust, decreasing = F),]#Sort by enrichment score from high to low
kegg_result <- KEGG[KEGG$p.adjust<0.05 & abs(KEGG$NES)>1] ###
write.csv(kegg_result,file="kegg_result.csv")


paths <- rownames(sortKEGG[c(1:5),])#Select the channel ID to be displayed
p <-  gseaplot2(KEGG,c(1:5),base_size =10,
                color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
                rel_heights = c(1.5, 0.3, 0.5),
                title="KEGG")

png("KEGG.png",w=8,h=5,units = "in",res = 300)
print(p)
dev.off()
pdf("KEGG.pdf",w=8,h=5)
print(p)
dev.off()









