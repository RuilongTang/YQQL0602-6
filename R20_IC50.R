rm(list = ls())

setwd("/data/nas1/yuanyt/project/05.YQQL-0602-6/")
if (! dir.exists("./20_IC50")){
  dir.create("./20_IC50")
}
setwd("./20_IC50")


###Read data -------------------
dat.tcga <- read.csv("../00_Rowdata/fpkm.tumor.csv",header=TRUE, row.names=1, check.names=FALSE)
dat.tcga <- log2(dat.tcga+1)
riskScore <- readRDS("../08_LASSO/riskScore.rds")
riskScore$sample <- rownames(riskScore)
riskScore$risk <- ifelse(riskScore$riskScore > median(riskScore$riskScore),"High Risk","Low Risk")
dat.tcga <- dat.tcga[,riskScore$sample]

###Calculate IC50
library(pRRophetic)
library(ggplot2)
set.seed(12345)
drug<-read.table(file = '/data/nas1/luchunlin/pipeline/Medicinal_Sensity/drugs.txt',sep='\t',header=F)
a <- data.frame(row.names=riskScore$sample,riskScore$sample)
cnt<-1
while (cnt < 139) {
  predictedPtype <- pRRopheticPredict(as.matrix(dat.tcga), drug[cnt,],selection=1)
  Tipifarnib<-data.frame(predictedPtype)
  colnames(Tipifarnib)<-drug[cnt,]
  a<-cbind(a,Tipifarnib)
  cnt = cnt + 1
}
write.csv(a,'IC50.csv')

a <- read.csv("IC50.csv",header = TRUE,row.names = 1)
colnames(a)[1] <- "id"

result <- a[,-1]
# write it in the form of a function for easy invocation
removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
result<-removeColsAllNa(result)
na_flag <- apply(is.na(result), 2, sum)
result <- result[, which(na_flag == 0)]
dim(result)

medicinal_result <- t(result)
####Calculate according to the grouped drawing
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)

high_group <- riskScore$sample[which(riskScore$risk =="High Risk")]
low_group <- riskScore$sample[which(riskScore$risk =="Low Risk")]


pvalue = padj = log2FoldChange <- matrix(0, nrow(medicinal_result), 1)
for (j in 1:nrow(medicinal_result)){
  pvalue[j, 1] = p.value = wilcox.test(medicinal_result[j, high_group],
                                       medicinal_result[j, low_group])$p.value
  log2FoldChange[j, 1] = mean(medicinal_result[j, high_group]) - 
    mean(medicinal_result[j, low_group])
}
padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
rTable <- data.frame(log2FoldChange, 
                     pvalue, 
                     padj,
                     row.names = rownames(medicinal_result))
high_group_res <- signif(apply(medicinal_result[rownames(rTable), high_group], 
                               1,
                               median), 4)
low_group_res <- signif(apply(medicinal_result[rownames(rTable), low_group], 
                              1, 
                              median), 4)
rTable <- data.frame(high_group_res, 
                     low_group_res,
                     rTable[, c("padj", "pvalue", "log2FoldChange")])
rTable$drugs <- rownames(rTable)
rTable$sig <- ifelse(rTable$pvalue < 0.05,
                     ifelse(rTable$pvalue < 0.01, 
                            ifelse(rTable$pvalue < 0.001,
                                   ifelse(rTable$pvalue < 0.0001,
                                          paste(rTable$drugs, "****",  sep = ""),
                                          paste(rTable$drugs, "***", sep = "")),
                                   paste(rTable$drugs, "**", sep = "")),
                            paste(rTable$drugs, "*",  sep = "")), 
                     rTable$drugs)
write.csv(rTable, file = "02.ic50_wilcox_test.csv")
DE.drug<-rTable[which(rTable$pvalue<0.05),]
DE.drug_up <- DE.drug[which(DE.drug$log2FoldChange>0),]
DE.drug_down <- DE.drug[which(DE.drug$log2FoldChange<0),]
write.csv(DE.drug_up, file = "03.ic50_up.csv")
write.csv(DE.drug_down, file = "04.ic50_down.csv")

all(rownames(rTable) == rownames(medicinal_result))
##Down-regulated
if (length(rownames(DE.drug_down)) != 0) {
  drugs_res <- data.frame(drugs=rownames(medicinal_result), medicinal_result, pvalue=rTable$pvalue)
  drugs_res <- drugs_res[rownames(drugs_res)%in%rownames(DE.drug_down),] 
  drugs_res <- drugs_res[which(rownames(drugs_res) %in% rownames(DE.drug)),]
  violin_dat <- gather(drugs_res, key=indivs, value=score, -c("drugs","pvalue"))
  violin_dat$indivs <- ifelse(gsub("\\.","-",violin_dat$indivs) %in% high_group,
                              "High Risk", "Low Risk") 
  violin_dat$indivs <- factor(violin_dat$indivs, levels = c("High Risk", "Low Risk"))
  ###head(violin_dat)
  drugs_hub_boxplot1 <- ggboxplot(violin_dat, x = "indivs", y = "score",
                                  color = "indivs", palette = c("#A73030FF", "#0073C2FF"),
                                  add = "jitter",
                                  short.panel.labs = T,
                                  ggtheme = theme_bw()) +
    stat_compare_means(label = "p.signif", label.x = 1.4, vjust = 0.5)
  drugs_hub_boxplot <- facet(drugs_hub_boxplot1,
                             facet.by = "drugs",
                             short.panel.labs = T,
                             panel.labs.background = list(fill = "white"),
                             ncol = 7,
                             scales = "free_y") + xlab("") + ylab("IC(50)") +
    # geom_text(data=data_text,
    #           mapping=aes(x=x,y=y,label=label),nudge_x=0.1,nudge_y=0.1)+
    theme(panel.grid = element_blank(),
          legend.position = "none",
          # strip.background = element_blank(),
          axis.title.y = element_text(size = 18, face = "bold"),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(size = 10, face = "bold",angle = 45,hjust = 1),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          text = element_text(size = 13, face = "bold"))
  drugs_hub_boxplot
  ggsave(filename = "ic50_drugs.plot(down).pdf", height = 10, width = 16)
  ggsave(filename = "ic50_drugs.plot(down).png", height = 10, width = 16)
}


###Up-regulated
if (length(rownames(DE.drug_up)) != 0){
  drugs_res <- data.frame(drugs=rownames(medicinal_result), medicinal_result, pvalue=rTable$pvalue)
  drugs_res <- drugs_res[rownames(drugs_res)%in%rownames(DE.drug_up),]
  # drugs_res <- drugs_res[rownames(rTable2),]
  drugs_res <- drugs_res[which(rownames(drugs_res) %in% rownames(DE.drug)),]
  violin_dat <- gather(drugs_res, key=indivs, value=score, -c("drugs","pvalue"))
  violin_dat$indivs <- ifelse(gsub("\\.","-",violin_dat$indivs) %in% high_group,
                              "High Risk", "Low Risk") 
  violin_dat$indivs <- factor(violin_dat$indivs, levels = c("High Risk", "Low Risk"))
  head(violin_dat)
  drugs_hub_boxplot1 <- ggboxplot(violin_dat, x = "indivs", y = "score",
                                  color = "indivs", palette = c("#A73030FF", "#0073C2FF"),
                                  add = "jitter",
                                  short.panel.labs = T,
                                  ggtheme = theme_bw()) +
    stat_compare_means(label = "p.signif", label.x = 1.4, vjust = 0.5)
  drugs_hub_boxplot <- facet(drugs_hub_boxplot1,
                             facet.by = "drugs",
                             short.panel.labs = T,
                             panel.labs.background = list(fill = "white"),
                             ncol = 7,
                             scales = "free_y") + xlab("") + ylab("IC(50)") +
    # geom_text(data=data_text,
    #           mapping=aes(x=x,y=y,label=label),nudge_x=0.1,nudge_y=0.1)+
    theme(panel.grid = element_blank(),
          legend.position = "none",
          # strip.background = element_blank(),
          axis.title.y = element_text(size = 18, face = "bold"),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(size = 10, face = "bold",angle = 45,hjust = 1),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          text = element_text(size = 13, face = "bold"))
  drugs_hub_boxplot
  ggsave(filename = "ic50_drugs.plot(up).pdf", height = 10, width = 16)
  ggsave(filename = "ic50_drugs.plot(up).png", height = 10, width = 16)
}


