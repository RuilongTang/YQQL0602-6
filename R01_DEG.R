


rm(list = ls())

setwd("")
if (! dir.exists("./01_DEGs")){
  dir.create("./01_DEGs")
}
setwd("./01_DEGs")


## 01 Read data----------
data <- read.csv("../00_Rowdata/dat.tcga.csv",header=TRUE, row.names=1, check.names=FALSE)
group <- read.csv("../00_Rowdata/tcga.group.csv",header = T,row.names=1)
dat.final <- data
# 02 Difference analysis----------
library(DESeq2)

##02-1 dentification of differentially expressed genes
# Construct the dds matrix
colData<-group
colData$group<-factor(colData$group,levels = c("Control","CRC"))
dds<-DESeqDataSetFromMatrix(countData = round(dat.final),colData=colData,design = ~group)
dds = dds[rownames(counts(dds)) > 1,]

# Extract the standardized data------
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds,normalized=T) 
## Normalize the original dds
dds<-DESeq(dds)

## Display dds information
dds

## Extract the DESeq2 analysis results
## Use the Result function to extract the difference analysis results
## Define the extracted difference analysis results as the variable "res"
## contrast: Defines who is compared with whom
res =results(dds, contrast = c("group","CRC","Control"))


## The result res is sorted by the pvalue using the order () function.
res =res[order(res$pvalue),]
head(res)
summary(res)
## Save all output results
write.csv(res,file="All_results.csv")
saveRDS(res, "DEG_raw.rds")
res <- readRDS("DEG_raw.rds")

table(res$pvalue<0.05)    # Show a significant difference in the number of genes
##FALSE  TRUE 
##15631  2392

# Extract and save the results of significant differences
## Obtain differentially expressed genes with padj< 0.05 and expression multiples greater than 1 or less than -1 after taking 2 as the logarithm.
## Use the subset () function to filter the required results into the variable DEG
## Usage:subset(x,...) x represents objects,... For the parameters or conditions of screening


DEG <- subset(res, pvalue < 0.05 & abs(log2FoldChange) >0.5 )
DEG<-as.data.frame(res)
DEG<-na.omit(DEG)
## Use dim to view the dimension and scale of this result
dim(DEG)
head(DEG)
## Add the "change" column
logFC_cutoff<-0.5
DEG$change=as.factor(
  ifelse(DEG$pvalue<0.05&abs(DEG$log2FoldChange)>logFC_cutoff,
         ifelse(DEG$log2FoldChange>logFC_cutoff,'UP','DOWN'),'NOT'))
table(DEG$change)
#  DOWN   NOT    UP  
# 231 15996   748
sig_diff <- subset(DEG,
                   DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) >= logFC_cutoff)
## 979
DEG_write <- cbind(GeneSymbol=rownames(DEG), DEG)
write.csv(DEG_write, file = "DEG_all.csv")
sig_diff_write <- cbind(GeneSymbol=rownames(sig_diff), sig_diff)
write.csv(sig_diff_write, file = "DEG_sig.csv")


## Volcano map------
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)
dat_rep<-DEG[rownames(DEG)%in%
               c(head(rownames(subset(sig_diff,sig_diff$log2FoldChange>0.5)),10),
                 head(rownames(subset(sig_diff,sig_diff$log2FoldChange< -0.5)),10)),]
volcano_plot<- ggplot(data = DEG, 
                      aes(x = log2FoldChange,
                          y = -log10(pvalue), 
                          color =change)) +
  scale_color_manual(values = c("#87CEFA", "gray","#F08080")) +
  scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0,1,5,10,20,50, 100,200)) +
  geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = c(-0.5,0.5),
             lty = 4,
             col = "black",  ####"darkgray"
             lwd = 0.6)+
  geom_hline(yintercept = -log10(0.05),
             lty = 4,
             col = "black",  ####"darkgray"
             lwd = 0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold",
                                   color="black",
                                   family = "Times",
                                   size=13),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.title.x = element_text(face = "bold",
                                    color = "black",
                                    size = 15),
        axis.title.y = element_text(face = "bold",
                                    color = "black",
                                    size = 15)) +
  geom_label_repel(
    data = dat_rep,
    aes(label = rownames(dat_rep)),
    max.overlaps = 20,
    size = 3,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log (Fold Change)",
       y = "-log10 (Pvalue)")
volcano_plot
ggsave('volcano.png', volcano_plot,width = 8, height = 7,dpi = 600)
ggsave('volcano.pdf', volcano_plot,width = 8, height = 7,dpi = 600)
dev.off()

####Complex heat map-------------------------
library(ComplexHeatmap)


# rt<-read.csv("../00_Rowdata/fpkm.tcga.csv",header=TRUE, row.names=1, check.names=FALSE)
rt<-normalized_counts
sig_diff1 <- sig_diff[order(sig_diff$padj),]
heat <- rt[head(rownames(sig_diff1),20),]
heat<-log2(heat+1)


# mat <- as.matrix(heat)
mat <- t(scale(t(heat)))
mat[mat < (-2)] <- (-2)
mat[mat > 2] <- 2

group$group <- factor(group$group, levels = c("Control", "CRC"))
group <- group[order(group$group),]


pdf('02.heatmap.pdf',  w=6,h=6,family='Times')
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(3, "cm")) %v%
  HeatmapAnnotation(Group = group$group, col = list(Group = c("CRC" = "#B72230", "Control" = "#104680"))) %v%
  Heatmap(mat, 
          row_names_gp = gpar(fontsize = 9),
          show_column_names = F,
          show_row_names = F,
          ####show_colnames = FALSE,
          name = "expression", 
          ### cluster_cols = F,
          cluster_rows = T,
          height = unit(6, "cm"),
          #cluster_columns = FALSE,
          ###cluster_rows = FALSE,
          col = colorRampPalette(c("#0A878D", "white","#D80305"))(100))
dev.off()

png('02.heatmap.png',w=6,h=6,units='in',res=600,family='Times')
densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",height = unit(3, "cm")) %v%
  HeatmapAnnotation(Group = group$group, col = list(Group = c("CRC" = "#B72230", "Control" = "#104680"))) %v%
  Heatmap(mat, 
          row_names_gp = gpar(fontsize = 9),
          show_column_names = F,
          show_row_names = F,
          ###show_colnames = FALSE,
          name = "expression", 
          ###cluster_cols = F,
          cluster_rows = T,
          height = unit(6, "cm"),
          #cluster_columns = FALSE,
          ###cluster_rows = FALSE,
          col = colorRampPalette(c("#0A878D", "white","#D80305"))(100))
dev.off()


