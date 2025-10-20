rm(list = ls())

setwd("/data/nas1/yuanyt/project/05.YQQL-0602-6/")
if (! dir.exists("./15_checkpoint")){
  dir.create("./15_checkpoint")
}
setwd("./15_checkpoint")

library(tidyverse)
library(lance)
library(readxl)
## 01 Read data ----------
data <- read.csv("../00_Rowdata/fpkm.tumor.csv",header=TRUE, row.names=1, check.names=FALSE)
riskScore <- readRDS("../08_LASSO/riskScore.rds")
riskScore$sample <- rownames(riskScore)
riskScore$risk <- ifelse(riskScore$riskScore > median(riskScore$riskScore),"High Risk","Low Risk")
high.sample<-riskScore$sample[which(riskScore$risk=='High Risk')]
low.sample<-riskScore$sample[which(riskScore$risk=='Low Risk')]
data <- log2(data+1)
fpkm_data <- data[,riskScore$sample]


dim(fpkm_data)
## Extract the test point genes
imm_gene <- read_xlsx('/data/nas1/liky/database/Immune_check.xlsx')
# checkpoint <- data.frame(V1=c('PD-L1','PD-1','TIM3','CTLA4','B7-H4',
#                               'B7-H3','BTLA','VISTA','IDO1','PSGL-1',
#                               'LAG3','PD-L2','OX40','IDO2','TNFRSF8',
#                               'CD27','ICOS','TNFRSF18','TIGIT','TNFRSF9',
#                               'TNFRSF14','TNFSF4','CD28','LGALS9','CD70','CD80',
#                               'TNFSF15','NRP1','BTNL2','HHLA2','ICOSLG',
#                               'CD40','TNFSF9','TNFSF14','CD86','KIR3DL1',
#                               'CD200','ADORA2A','TNFRSF25','CD244',
#                               'CD48','LAIR1','CD40LG','TMIGD2','CD200R1',
#                               'TNFSF18','CD44','CD160'))

checkpoint <- imm_gene$symbol
length(checkpoint)
checkpoint_dat <- fpkm_data[which(rownames(fpkm_data)%in%checkpoint),]
dim(checkpoint_dat)
checkpoint_dat$gene<-rownames(checkpoint_dat)

library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
head(checkpoint_dat[,1:3])
violin_dat <- gather(checkpoint_dat, key=sample, value='expr', -c("gene"))
head(violin_dat)
violin_dat$group <- ifelse(violin_dat$sample %in% high.sample,
                           "High Risk", "Low Risk") 
head(violin_dat)
colnames(violin_dat)
violin_dat$group <- factor(violin_dat$group, levels = c("Low Risk","High Risk"))
#Difference analysis --------------
library(rstatix)
stat.test<-violin_dat%>%
  group_by(gene)%>%
  wilcox_test(expr ~ group)%>%
  adjust_pvalue(method = 'fdr')
# stat.test$p.adj<-ifelse(stat.test$p<0.001,"***",
#                         ifelse(stat.test$p<0.01,"**",
#                                ifelse(stat.test$p<0.05,"*",'ns')))
write.csv(stat.test,file = '01.checkpoint_result.csv',row.names = F)

DE.checkpoint <- subset(stat.test,stat.test$p < 0.05)
##40
violin_dat<-violin_dat[violin_dat$gene%in%DE.checkpoint$gene,]
violin_dat$gene <- factor(violin_dat$gene,levels = DE.checkpoint$gene)


violin_plot <- ggplot(violin_dat, aes(x=gene, 
                                      y=expr,
                                      fill=group)) +
  #geom_violin(trim=F,color="black") + #To draw a violin picture, use "color=" to set the color of the outline of the violin picture. (# No outline can be set to white. Below, setting the background to white actually indicates no outline.)
  #If "trim" is TRUE(the default value), the tail of the violin will be trimmed to the data range. If it is FALSE, do not trim the tail.
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               aes(fill=group),
               position=position_dodge(0.9),
               outlier.shape = NA)+ #Draw the box plot. Here, width=0.1 controls the width of the box plot in the violin plot
  scale_fill_manual(values= c("#658fcb", "#f38687"), name = "Group")+
  labs(title="Immune Checkpoint", x="", y = "Expression level",size=20) +
  stat_compare_means(data = violin_dat,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  # stat_pvalue_manual(DE.checkpoint,
  #                    x ="gene",
  #                    y.position = 7,
  #                    size = 3.2,
  #                    color = "black",
  #                    family = "Times",
  #                    label = "p.adj")+
  ylim(c(0,9))+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=8),
        legend.title = element_text(face = "bold", size = 10),
        text = element_text(family = 'Times'),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())#+
#facet_wrap(~gene,scales = "free",nrow = 2) #Separate or combine them together, two lines or one line
violin_plot
ggsave('01.checkpoint.pdf',violin_plot,w=10,h=6)
ggsave('01.checkpoint.png',violin_plot,w=10,h=6)


#Correlation analysis of prognostic genes and immune checkpoints---------------
library(Hmisc)
library(dplyr)


gene <- readRDS("../08_LASSO/cox_gene.rds")
gene <- fpkm_data[gene,]%>%t%>%as.data.frame()

de_check <- checkpoint_dat[DE.checkpoint$gene,-607]%>%t%>%as.data.frame()


####Genes + checkpoints/samples
nc <- cbind(gene,de_check)
nc <- as.matrix(nc)

#####Correlation analysis
m=rcorr(nc,type = 'spearman')$r[1:ncol(gene),(ncol(gene)+1):ncol(nc)]
p=rcorr(nc,type = 'spearman')$P[1:ncol(gene),(ncol(gene)+1):ncol(nc)]


tmp <- p
tmp <- ifelse(tmp >0.05,"ns",ifelse(tmp >0.01,"*",ifelse(tmp >0.001,"**",ifelse(tmp >0.0001,"***","****"))))


cor <- m
cor <- round(cor,digits = 3)
textMatrix = paste(cor,"\n",tmp, sep = "")
dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)
textMatrix <- t(textMatrix)
write.csv(textMatrix,file = 'correlaion.csv')
m1 <- t(m)

library(WGCNA)
pdf(file = '02.correlation.pdf',w=7,h=10)
par(mar = c(5, 5, 3, 3));
labeledHeatmap(Matrix = m1, 
               xLabels = colnames(m1), 
               yLabels = rownames(m1), 
               cex.lab = 0.8, 
               # ySymbols = colnames(MEs_col)[-20], 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50, 0.8),
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.8, 
               zlim = c(-1,1),
               # main = paste("correlation")
)
dev.off()
png(file = '02.correlation.png',w=7,h=10,units = 'in',res = 600)
par(mar = c(5, 5, 3, 3));
labeledHeatmap(Matrix = m1, 
               xLabels = colnames(m1), 
               yLabels = rownames(m1), 
               cex.lab = 0.8, 
               # ySymbols = colnames(MEs_col)[-20], 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50, 0.8),
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.8, 
               zlim = c(-1,1),
               # main = paste("correlation")
)
dev.off()











