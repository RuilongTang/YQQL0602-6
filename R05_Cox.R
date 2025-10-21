rm(list = ls())

setwd("")
if (! dir.exists("./05_Cox")){
  dir.create("./05_Cox")
}
setwd("./05_Cox")


# 01 Obtain the dataset----------
data <- read.csv("../00_Rowdata/fpkm.tumor.csv",header=TRUE, row.names=1, check.names=FALSE)
data <- log2(data+1)

##
gene <- read.csv("../04_Veen/candi.csv",header = T,row.names = 1)
survival <- read.csv("../02_ssGSEA/icd.group.csv",header = T,row.names = 1)
survival <- survival[,c(1,2,3)]


## Merge survival data
survival_dat<-t(data[gene$symbol,survival$sample])
train_dat<-t(scale(t(survival_dat))) %>% data.frame()
train_dat$sample<-rownames(train_dat)
train_dat<-merge(survival,train_dat,by='sample')
rownames(train_dat)<-train_dat$sample
train_dat<-train_dat[,-1]
colnames(train_dat)


### Univariate cox
library(survival)
library(survminer)
colnames_sum <- colnames(train_dat)
colnames_sum <- gsub("-","_",colnames_sum)
colnames_sum <- gsub(" ","_",colnames_sum)
colnames(train_dat) <- colnames_sum
covariates <- colnames_sum[-which(colnames_sum %in% c("OS", "OS.time"))]

#The Surv() function generates an impact of the survival time of a surviving object on survival and constructs a survival analysis formula for each variable
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste("Surv(OS.time, OS)~", x)))  #as.formula(). Converts a string into a formula. Construct the formula object
# The coxph function is used to calculate the cox model cycle for cox regression analysis on each feature
univ_models <- lapply(univ_formulas,
                      function(x) {coxph(x, data = train_dat)})

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #Obtain p value
                         p.value<-signif(x$wald["pvalue"], digits=3)
                         #Obtain HR
                         HR <-signif(x$coef[2], digits=3);
                         #Obtain the 95% confidence interval
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", 
                                      HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })

## coef is the regression coefficient b in the formula (sometimes also called the beta value). exp(coef) is the hazard ratio (HR) in the cox model.
## z represents the wald statistic, which is coef divided by its standard error se(coef). ower.95 and upper.95 represent the 95% confidence interval of exp(coef). The narrower the confidence interval, the higher the confidence, and the more precise and truthful your experiment is.
res_mod <- t(as.data.frame(univ_results, check.names = FALSE))
res_mod <- as.data.frame(res_mod)
##0.05
res_results_0.05 <- res_mod[which(as.numeric(res_mod$p.value) < 0.05),]
res_results_0.05 <- na.omit(res_results_0.05)
write.csv(res_results_0.05, file = "univariate_cox_result_0.05.csv")
dim(res_results_0.05)
library(tidyr)
res_results_0.05_2 <- separate(res_results_0.05, "HR (95% CI for HR)",
                               into = c("HR", "HR.95L"),
                               sep = " ")
res_results_0.05_2 <- separate(res_results_0.05_2, "HR.95L",
                               into = c("HR.95L", "HR.95H"),
                               sep = "\\-")
res_results_0.05_2$HR.95L <- gsub("\\(", "", res_results_0.05_2$HR.95L)
res_results_0.05_2$HR.95H <- gsub("\\)", "", res_results_0.05_2$HR.95H)

res_results_0.05_2[,1:ncol(res_results_0.05_2)] <- as.numeric(unlist(res_results_0.05_2[,1:ncol(res_results_0.05_2)]))


hz <- paste(round(res_results_0.05_2$HR,3),
            "(",round(res_results_0.05_2$HR.95L,3),
            "-",round(res_results_0.05_2$HR.95H,3),")",sep = "")


tabletext <- cbind(c(NA,"Gene",rownames(res_results_0.05_2)),
                   c(NA,"P value",ifelse(res_results_0.05_2$p.value<0.001,
                                         "< 0.001",
                                         round(res_results_0.05_2$p.value,3))),
                   c(NA,"Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
library(forestplot)
pdf(file = "univariate_cox_forest.pdf", height = 6, width = 10, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #It is the position where the Pvalue box plot is located
           is.summary = c(TRUE, TRUE,rep(FALSE, 57)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,NA,res_results_0.05_2$HR),
           lower=c(NA,NA,res_results_0.05_2$HR.95L), #Lower limit of the 95% confidence interval
           upper=c(NA,NA,res_results_0.05_2$HR.95H), #Upper limit of the 95% confidence interval
           boxsize=0.2,lwd.ci=3,   #The size of the box and the width of the line
           ci.vertices.height = 0.08,ci.vertices=TRUE, #Confidence intervals are defined by line width, height and shape
           zero=1,lwd.zero=0.5,    #The position of the zero line width reference line
           colgap=unit(5,"mm"),    #Column gap
           xticks = c(0, 1,2,3,4,5,6,7), #Horizontal coordinate scale
           lwd.xaxis=2,            #Width of the X-axis
           lineheight = unit(0.8,"cm"), #Fixed row height
           graphwidth = unit(.5,"npc"), #The width ratio of the figure in the table
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #Error bar display mode
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #A black line is added at the top of the third line, and the numbers within the quotation marks indicate the position of the line
           #                 "59" = gpar(lwd=2, col="black")),#At the bottom of the last line, add a black line. The number in "16" is nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#Graphic margins
           #The cex parameter in the fpTxtGp function sets the size of each component
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T) # The grid lines perpendicular to the X-axis correspond to each scale
dev.off()
png(filename = "univariate_cox_forest.png", height = 600, width = 800)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #It is the position where the Pvalue box plot is located
           is.summary = c(TRUE, TRUE,rep(FALSE, 57)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,NA,res_results_0.05_2$HR),
           lower=c(NA,NA,res_results_0.05_2$HR.95L), #Lower limit of the 95% confidence interval
           upper=c(NA,NA,res_results_0.05_2$HR.95H), #Upper limit of the 95% confidence interval
           boxsize=0.2,lwd.ci=3,   #The size of the box and the width of the line
           ci.vertices.height = 0.08,ci.vertices=TRUE, #Confidence intervals are defined by line width, height and shape
           zero=1,lwd.zero=0.5,    #The position of the zero line width reference line
           colgap=unit(5,"mm"),    #Column gap
           xticks = c(0, 1,2,3,4,5,6,7), #Horizontal coordinate scale
           lwd.xaxis=2,            #Width of the X-axis
           lineheight = unit(0.8,"cm"), #Fixed row height
           graphwidth = unit(.5,"npc"), #The width ratio of the figure in the table
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #Error bar display mode
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #A black line is added at the top of the third line, and the numbers within the quotation marks indicate the position of the line
           #                 "59" = gpar(lwd=2, col="black")),#At the bottom of the last line, add a black line. The number in "16" is nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#Graphic margins
           #The cex parameter in the fpTxtGp function sets the size of each component
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T) # The grid lines perpendicular to the X-axis correspond to each scale

dev.off()

