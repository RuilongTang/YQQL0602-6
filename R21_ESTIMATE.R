rm(list = ls())

setwd("/data/nas1/yuanyt/project/05.YQQL-0602-6/")
if (! dir.exists("./21_ESTIMATE")){
  dir.create("./21_ESTIMATE")
}
setwd("./21_ESTIMATE")



###Read data-------------------------
res.cibersort <- readRDS("../14_CIBERSORT/cibersort.rds")
dat.cibersort <- res.cibersort %>% tibble::column_to_rownames(var = "cell_type")
###Delete rows that are all zeros
dat.cibersort <- dat.cibersort[which(rowSums(abs(dat.cibersort)) > 0),]
riskScore <- readRDS("../08_LASSO/riskScore.rds")
riskScore$sample <- rownames(riskScore)


#Sort out the data --------------------------
group <- riskScore
group$label <- ifelse(group$riskScore > median(group$riskScore),"High Risk","Low Risk")
Low.sample<-group$sample[which(group$label=='Low Risk')]
High.sample<-group$sample[which(group$label=='High Risk')]


group_estimate<-group$label%>%as.factor()
design<-model.matrix(~0 + group_estimate)
rownames(design)<-group$sample
colnames(design)<-levels(group_estimate)
design<-as.data.frame(design)

Low<-rownames(design)[which(design$`Low Risk`==1)]
High<-rownames(design)[which(design$`High Risk`==1)]
length(Low)
length(High)

#estimate analysis -----------------------
#install.packages("estimate", repos="http://R-Forge.R-project.org")
library(estimate)
expr_train <- data

write.table(expr_train,
            'expr.txt',
            col.names = T,
            row.names = T,
            quote = F, sep="\t")
#Generate expr_train.gct
filterCommonGenes(input.f = './expr.txt', 
                  output.f = 'expr_train.gct', 
                  id = 'GeneSymbol')
# [1] "Merged dataset includes 9872 genes (540 mismatched)."

# Generate train_purity.gct
estimateScore('expr_train.gct', 'train_purity.gct', platform="illumina")
es_score <- read.table('train_purity.gct', skip = 2, header = T, check.names = F)
colnames(es_score)<-gsub('.','-',colnames(es_score),fixed = T)
immu_score <- es_score[,3:length(es_score)]
rownames(immu_score) <- es_score$NAME
write.csv(immu_score,file = "01.es_score.csv")

###Combine the risk score, immunity score and matrix score ---------
immu_score <- read.csv("01.es_score.csv",header = T,row.names = 1, check.names=FALSE)
immu_score <- t(immu_score) %>% data.frame()
immu_score$sample <- rownames(immu_score)

heat.group <- merge(immu_score, riskScore, by='sample')
heat.group$risk <- ifelse(heat.group$riskScore > median(heat.group$riskScore),"High Risk","Low Risk")


####Draw a heat map -----------------
colnames(heat.group)
heat.group<-dplyr::select(heat.group,c('risk','riskScore','StromalScore','ImmuneScore','ESTIMATEScore','sample'))
rownames(heat.group) <- heat.group$sample
heat.group <- heat.group[,-6]
heat.group <- heat.group[order(heat.group$riskScore),]



rt_dat<-dat.cibersort
rt_dat <- rt_dat[,rownames(heat.group)]
rt_dat <- data.frame(t(scale(t(rt_dat))))
colnames(rt_dat)<-gsub('.','-',colnames(rt_dat),fixed = T)
rt_dat[rt_dat < (-2)] <- (-2)
rt_dat[rt_dat > 2] <- 2



ann_colors<-list(risk = c("High Risk" = "#2E8B57","Low Risk" = "#F5F5F5"),
  riskScore = c("white","#2E8B57"),
  StromalScore = c("white","#7B68EE"),
  ImmuneScore = c("white","#DB7093"),
  ESTIMATEScore = c("white","#F4A460")
)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
hp <- pheatmap::pheatmap(rt_dat,
                         color = bluered(100),
                         border_color = NA,
                         annotation_col = heat.group,
                         annotation_colors = ann_colors,
                         labels_row = NULL,
                         clustering_method = 'ward.D2',
                         show_rownames = T,
                         show_colnames = F,
                         fontsize_col = 5,
                         cluster_cols = F,
                         cluster_rows = T)
png("heatmap.png",w=12,h=8,units = "in",res = 300)
print(hp)
dev.off()
pdf("heatmap.pdf",w=12,h=8)
print(hp)
dev.off()

