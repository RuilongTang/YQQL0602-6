rm(list = ls())

setwd("/data/nas1/yuanyt/project/05.YQQL-0602-6/")
if (! dir.exists("./18_MSI")){
  dir.create("./18_MSI")
}
setwd("./18_MSI/new_checkpoint")


###Read data -------------------------
data <- read.csv("../../00_Rowdata/fpkm.tumor.csv",header=TRUE, row.names=1, check.names=FALSE)
data <- log2(data+1)
checkpoint <- c("PDCD1","CD274","CTLA4")
checkpoint_dat <- data[which(rownames(data)%in%checkpoint),]


survival <- read.csv("../../02_ssGSEA/icd.group.csv",header = T,row.names = 1)
survival <- survival[,c(1,2,3)]
riskScore <- readRDS("../../08_LASSO/riskScore.rds")
riskScore$sample <- rownames(riskScore)
riskScore <- merge(riskScore,survival,by= "sample")
riskScore$risk <- ifelse(riskScore$riskScore > median(riskScore$riskScore),"High Risk","Low Risk")

phenotype_coad <- read_tsv(file = '../../10_progmodel/TCGA-COAD.GDC_phenotype.tsv.gz')
phenotype_read <- read_tsv(file = '../../10_progmodel/TCGA-READ.GDC_phenotype.tsv.gz')
coad <- phenotype_coad[,colnames(phenotype_coad) %in% colnames(phenotype_read)]
read <- phenotype_read[,colnames(phenotype_read) %in% colnames(phenotype_coad)]

phenotype <- rbind(coad,read)
msi <- data.frame(sample = phenotype$submitter_id.samples,
                  MSI=phenotype$microsatellite_instability)
msi <- merge(msi, riskScore, by='sample')
msi1 <- na.omit(msi)
high.sample<-msi1$sample[which(msi1$MSI=='YES')]
low.sample<-msi1$sample[which(msi1$MSI=='NO')]



checkpoint_dat1 <- checkpoint_dat[,msi1$sample]
checkpoint_dat1$gene<-rownames(checkpoint_dat1)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
violin_dat <- gather(checkpoint_dat1, key=sample, value='expr', -c("gene"))
head(violin_dat)
violin_dat$group <- ifelse(violin_dat$sample %in% high.sample,
                           "MSI-H", "MSI-L/MSS") 
head(violin_dat)
colnames(violin_dat)
violin_dat$group <- factor(violin_dat$group, levels = c("MSI-H", "MSI-L/MSS"))
#Difference analysis --------------
library(rstatix)
stat.test<-violin_dat%>%
  group_by(gene)%>%
  wilcox_test(expr ~ group)%>%
  adjust_pvalue(method = 'fdr')


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
  scale_fill_manual(values= c("#f38687","#658fcb"), name = "Group")+
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
  # ylim(c(0,9))+
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
ggsave('01.checkpoint.pdf',violin_plot,w=6,h=5)
ggsave('01.checkpoint.png',violin_plot,w=6,h=5)

