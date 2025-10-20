rm(list = ls())

setwd("/data/nas1/yuanyt/project/05.YQQL-0602-6/")
if (! dir.exists("./20_IC50")){
  dir.create("./20_IC50")
}
setwd("./20_IC50")


###读取数据-------------------
dat.tcga <- read.csv("../00_Rowdata/fpkm.tumor.csv",header=TRUE, row.names=1, check.names=FALSE)
dat.tcga <- log2(dat.tcga+1)
riskScore <- readRDS("../08_LASSO/riskScore.rds")
riskScore$sample <- rownames(riskScore)
riskScore$risk <- ifelse(riskScore$riskScore > median(riskScore$riskScore),"High Risk","Low Risk")
group<-data.frame(sample=riskScore$sample,group = riskScore$risk)
dat.tcga <- dat.tcga[,riskScore$sample]


library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
library(reshape2)
library(ggpubr)
drug <- c("oxaliplatin","xeloda","lrinotecan","bevacizum","Cetuximab")

th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
dir='/data/nas1/yuanyt/pipeline/GDSC/DataFiles/Training Data/'
GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
GDSC2_Res <- exp(GDSC2_Res) 

testExpr<- as.matrix(dat.tcga)
dim(testExpr)  

GDSC2_Res <- GDSC2_Res[,c('Oxaliplatin_1089',"Oxaliplatin_1806","Irinotecan_1088")]
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData = testExpr,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 5, 
              printOutput = TRUE, 
              removeLowVaringGenesFrom = 'rawData' )
result <- read.csv('calcPhenotype_Output/DrugPredictions.csv',row.names = 1)
drug <- log2(result+1)
drug$sample <-rownames(drug) 




### 差异分析--------
drug_plot <- drug
drug_plot <- merge(group,drug_plot,by='sample')
drug_plot2 <- tidyr::gather(drug_plot,Drug,IC50,-c('sample','group'))

library(rstatix)
library(ggplot2)
library(ggpubr)

res.drug <- drug_plot2%>%
  group_by(Drug)%>%
  wilcox_test(IC50 ~ group) %>%
  adjust_pvalue(method = 'BH')%>%
  add_significance('p')

write.csv(res.drug,file = 'res.drug.csv')
colnames(drug_plot2)

violin<-drug_plot2
drug_plot <- ggplot(violin, aes(x=group,
                                y=IC50,
                                fill=group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.3,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#FF6A6A", "#20B2AA"), name = "Group")+
  labs(title="", x="", y = "IC50",size=20) +
  stat_compare_means(data = violin,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,label.x = 1.4) +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=0,hjust=,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  facet_wrap(~Drug,scales = "free",nrow = 1) +
  guides(fill='none')
drug_plot
ggsave(filename = '01.Risk.drug.pdf',drug_plot,w=8,h=4)
ggsave(filename = '01.Risk.drug.png',drug_plot,w=8,h=4)


###MSI MSS---------------------
msi <- read.csv("../18_MSI/MSI.csv",header = T,row.names = 1)
group <- data.frame(sample = rownames(msi),group = msi$group)

### 差异分析--------
drug_plot <- drug[group$sample,]
drug_plot <- merge(group,drug_plot,by='sample')
drug_plot2 <- tidyr::gather(drug_plot,Drug,IC50,-c('sample','group'))

library(rstatix)
library(ggplot2)
library(ggpubr)

res.drug <- drug_plot2%>%
  group_by(Drug)%>%
  wilcox_test(IC50 ~ group) %>%
  adjust_pvalue(method = 'BH')%>%
  add_significance('p')

write.csv(res.drug,file = '02.msi.drug.csv')
colnames(drug_plot2)

violin<-drug_plot2
drug_plot <- ggplot(violin, aes(x=group,
                                y=IC50,
                                fill=group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.3,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#FF6A6A", "#20B2AA"), name = "Group")+
  labs(title="", x="", y = "IC50",size=20) +
  stat_compare_means(data = violin,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,label.x = 1.4) +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=0,hjust=,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  facet_wrap(~Drug,scales = "free",nrow = 1) +
  guides(fill='none')
drug_plot
ggsave(filename = '02.MSI.drug.pdf',drug_plot,w=8,h=4)
ggsave(filename = '02.MSI.drug.png',drug_plot,w=8,h=4)






