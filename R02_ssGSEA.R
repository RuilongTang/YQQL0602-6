rm(list = ls())

setwd("")
if (! dir.exists("./02_ssGSEA")){
  dir.create("./02_ssGSEA")
}
setwd("./02_ssGSEA")

## 01 Read data----------
data <- read.csv("../00_Rowdata/fpkm.tumor.csv",header=TRUE, row.names=1, check.names=FALSE)%>% lc.tableToNum
data <- log2(data+1)

##Survival data
survival_coad <- read_tsv("../00_Rowdata/TCGA-COAD.survival.tsv")
survival_read <- read_tsv("../00_Rowdata/TCGA-READ.survival.tsv")
survival <- rbind(survival_coad,survival_read)
survival <- survival[,c(1,2,4)]

saveRDS(survival,"survival.rds")
write.csv(survival,file = "tcga_survival.csv")


###ICD score---------------------
library(GSVA)

icd <- read.csv("../00_Rowdata/ICD.csv",header=TRUE)
dat.final2 <- as.matrix(data[,colnames(data)%in%survival$sample])

gene_list <- split(as.matrix(icd)[,1],
                   icd[,2])

ssgsea_score = gsva(dat.final2, gene_list, 
                    method = "ssgsea", 
                    ssgsea.norm = TRUE, 
                    verbose = TRUE)
write.csv(ssgsea_score,file = "ssgsea_result.csv")

score <- t(ssgsea_score) %>% data.frame()
score$sample <- rownames(score)



####Add survival curve
icd_survival <- merge(survival,score,by="sample")
icd_survival$group <- ifelse(icd_survival$immunogenic.cell.death>median(icd_survival$immunogenic.cell.death),'High','Low')

###Draw the survival curve
library("survival")
library("survminer")


###Cut-off value calculation
res.cut <- surv_cutpoint(icd_survival,time = 'OS.time',event = 'OS',variables = 'immunogenic.cell.death')
res.cut
cutpoint <- res.cut$cutpoint[1]
icd_survival$group <- ifelse(icd_survival$immunogenic.cell.death>cutpoint$cutpoint,'High','Low')
write.csv(icd_survival,file = "icd.group.csv")


###Calculate the survival value
kmfit <- survfit(Surv(OS.time,OS) ~ group, data = icd_survival)
print(kmfit)
##Drawing
pdf(file = 'KM.pdf',w=8,h=7,onefile = F)
ggsurvplot(kmfit,
           pval = TRUE,
           conf.int = F,
           legend.labs=c('High','Low'),
           legend.title="group",
           title="",
           font.main = c(15,"bold"),font.legend = 15,
           font.x = 15,font.y = 15,
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           surv.median.line = "hv",
           ggtheme = theme_bw(),
           palette = c("#A73030FF", "#0073C2FF"))
dev.off()

png(file = 'KM.png',w=700,h=600)
ggsurvplot(kmfit,
           pval = TRUE,
           conf.int = F,
           legend.labs=c('High','Low'),
           legend.title="group",
           title="",
           font.main = c(15,"bold"),font.legend = 15,
           font.x = 15,font.y = 15,
           risk.table = TRUE,
           risk.table.col = "strata",
           linetype = "strata",
           surv.median.line = "hv",
           ggtheme = theme_bw(),
           palette = c("#A73030FF", "#0073C2FF"))
dev.off()

