rm(list = ls())

setwd("/data/nas1/yuanyt/project/05.YQQL-0602-6/")
if (! dir.exists("./04_Veen")){
  dir.create("./04_Veen")
}
setwd("./04_Veen")


library(clusterProfiler)
DEG <- read.csv("../01_DEGs/DEG_sig.csv",header = T,row.names = 1)
brown <- read.csv("../03_WGCNA/brown.modgene.csv",header = T,row.names = 1)
blue <- read.csv("../03_WGCNA/blue.modgene.csv",header = T,row.names = 1)
wg <- rbind(brown,blue)

candi <- intersect(DEG$GeneSymbol,wg$modgene)%>%as.data.frame()
colnames(candi) <- "symbol"
write.csv(candi,file = 'candi.csv')



library(ggvenn)
mydata<-list('DEG'=DEG$X,'WGCNA'=wg$modgene)
pdf('01.venn.pdf',w=5,h=5)
ggvenn(mydata,c('DEG','WGCNA'),
       fill_color = c("#ffb2b2","#87CEEB"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 4,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#87CEEB"),
       text_color = 'black')
dev.off()
png('01.venn.png',w=400,h=400)
ggvenn(mydata,c('DEG','WGCNA'),
       fill_color = c("#ffb2b2","#87CEEB"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 4,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#87CEEB"),
       text_color = 'black')
dev.off()


###Read data——————
data <- read.csv("../01_DEGs/DEG_sig.csv", header=TRUE, row.names=1, check.names=FALSE)

####  GO enrichment
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
intersect <- candi
gene_transform <- bitr(intersect$symbol,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 1,
                qvalueCutoff = 1,
                readable = TRUE)
write.csv(ego,file = "GO.csv")

go_result <- data.frame(ego)

go_result <- go_result[,c(1,2,3,9,6)]#Only the data in columns 1,2,3,9 and 6 are needed
go_result$geneID <- str_replace_all(go_result$geneID,"/",",")#Replace the gene names in the geneID column with commas (/)
names(go_result) <- c("Category","ID","term","Genes","adj_pval")#Rename the column names
go_result <- go_result[order(go_result$adj_pval),]
go_result <- go_result[c(1:10),]

#####Obtain the differential genes and the corresponding FCS
logfc <- data[candi$symbol,]
logfc <- data.frame(gene = rownames(logfc),logFC=logfc$log2FoldChange)
colnames(logfc) <- c('ID', 'logFC')
row.names(logfc)=logfc[,1]
circ <- circle_dat(go_result, logfc) #Pay attention to the column names

go_bar <- GOCircle(circ, nsub = 10,
                   zsc.col = c('firebrick3', 'white', 'royalblue3'),
                   label.fontface='bold')
ggsave('02.GO.png',go_bar,width =16,height = 10)
ggsave('02.GO.pdf',go_bar,width =16,height = 10)




## KEGG enrichment（Bubble chart）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.csv(kk@result,file = "KEGG.csv")

kk_enrich <- kk@result
kk_enrich <- data.frame(Category = "KEGG",
                        ID = kk_enrich$ID,
                        Term = kk_enrich$Description, 
                        Genes = gsub("/", ", ",kk_enrich$geneID), 
                        adj_pval =kk_enrich$p.adjust)

circ <- circle_dat(kk_enrich, logfc) #Pay attention to the column names

interest_term_kk <- head(kk_enrich$Term, 5)
interest_gene_kk <- head(circ[,c("genes", "logFC")], 100)
interest_gene_kk <- unique(interest_gene_kk)
chord <- chord_dat(data = circ, genes=interest_gene_kk,
                   process = interest_term_kk )

# Change the variables of the GOChord program to make the image display more aesthetically pleasing
# (1) Display the enrichment pathway in a column in the image: 185 rows with ncol = 1
# (2) Show KEGG or other results: : Line 185 "GO Terms" change to "KEGG"
# (3) Make the Log2FC caption appear on the left and displayed vertically: 202 lines legend.direction = "vertical"
####trace(GOChord, edit = TRUE)

#trace(GOChord,edit=T)
p<-GOChord(
  data = chord,
  title = 'KEGG Pathways',
  space = 0.02,#Term spacing
  limit = c(0,1),
  gene.order = 'logFC',
  gene.space = 0.2,
  gene.size = 4,
  lfc.col = c('#D80305','white','#104680'), #Up-regulate gene color
  ribbon.col=brewer.pal(length(unique(head(kk_enrich$Term, 5))), "Set3"),#GO term color Settings
  process.label = 30
)
p
pdf('03.KEGG.pdf',w=14,h=15.5,family='Times')
print(p)
dev.off()

png('03.KEGG.png',w=14,h=15.5,family='Times',units='in',res=600)
print(p)
dev.off()





