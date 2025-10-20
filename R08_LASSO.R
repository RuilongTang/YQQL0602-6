rm(list = ls())

setwd("/data/nas1/yuanyt/project/05.YQQL-0602-6/")
if (! dir.exists("./08_LASSO")){
  dir.create("./08_LASSO")
}
setwd("./08_LASSO")


library(lance)
library(readxl)
library(readr)
library(tidyverse)
###Read data ——————
# 01 Get the dataset ----------
data <- read.csv("../00_Rowdata/fpkm.tumor.csv",header=TRUE, row.names=1, check.names=FALSE)
data <- log2(data+1)

##
gene <- read.csv("../05_Cox/univariate_cox_result_0.05.csv",header = T,row.names = 1)
survival <- read.csv("../02_ssGSEA/icd.group.csv",header = T,row.names = 1)
survival <- survival[,c(1,2,3)]


## Merge survival data
survival_dat<-t(data[rownames(gene),survival$sample])
train_dat<-survival_dat %>% data.frame()
train_dat$sample<-rownames(train_dat)
train_dat<-merge(survival,train_dat,by='sample')
rownames(train_dat)<-train_dat$sample
train_dat<-train_dat[,-1]
colnames(train_dat)



### LASSO
library(survival)
library(survminer)
library(glmnet)
train_data<-train_dat
x_all <- subset(train_data, select = -c(OS, OS.time))
y_all <- subset(train_data, select = c(OS, OS.time))



# Fitting model
fit <- glmnet(as.matrix(x_all), Surv(y_all$OS.time,y_all$OS), 
              family = "cox") 
#dev.new()
png(filename = "01.lasso_model.png", height = 400, width = 500)
plot(fit, xvar = "lambda",label = TRUE, las=1)
dev.off()
pdf(file = "01.lasso_model.pdf", height = 5)
plot(fit, xvar = "lambda",label = TRUE, las=1)
dev.off()


# Cross-validate the fitting model
set.seed(2)
# 2 6gene 
{
  cvfit = cv.glmnet(as.matrix(x_all),
                    Surv(y_all$OS.time,y_all$OS),nfold=10,
                    family = "cox") 
  png(filename = "02.lasso_verify.png", height = 400, width = 500)
  plot(cvfit, las =1)
  dev.off()
  pdf(file = "02.lasso_verify.pdf", height = 5)
  plot(cvfit, las =1)
  dev.off()
  coef.min = coef(cvfit, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
  cvfit$lambda.min
  # [1] 0.03887252
  # Find out those regression coefficients that have not been penalized as 0
  active.min = which(coef.min@i != 0)
  
  coef.min
  
  # Extract gene names
  lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
  lasso_geneids
}

# Multivariate cox ------------------------------------------------------------------

library(survival)
library(survminer)
cox_data <- as.formula(paste0('Surv(OS.time,OS)~', paste(lasso_geneids, sep = '', collapse = '+')))#svmrfe_result;lasso_geneids;cox_symbol
cox_more <- coxph(cox_data, data = train_data)
cox_zph <- cox.zph(cox_more)#
cox_table <- cox_zph$table[-nrow(cox_zph$table),]
cox_formula <- as.formula(paste("Surv(OS.time,OS)~",
                                paste(rownames(cox_table)[cox_table[,3] > 0.05],
                                      collapse = "+")))
cox_more_2 <- coxph(cox_formula, data = train_data)

mul_cox_result <- summary(cox_more_2)$coefficients %>% as.data.frame()
mul_cox_gene <- rownames(mul_cox_result)
mul_cox_gene

train_HRforest <- ggforest(model = cox_more_2,
                           data = train_data,
                           main = "Hazard ratio of Candidate Genes",
                           fontsize = 1) +
  theme(plot.title = element_text(face = "bold", size = 10))
train_HRforest
ggsave(width=10,height=6,'03.multiCox.pdf',train_HRforest)
ggsave(width=10,height=6,'03.multiCox.png',train_HRforest)

# step --------------------------------------------------------------------
step_Cox <- step(cox_more_2 ,direction = "both") ##The parameters also include \forward\both\backward
cox_summary <- summary(step_Cox)
cox_result <- summary(step_Cox)$coefficients
cox_gene <- rownames(cox_result)
cox_gene

write.csv(cox_result, file = "multicox_coefficients_step.txt")

train_HRforest <- ggforest(model = step_Cox,
                           data = train_data,
                           main = "Hazard ratio of Candidate Genes",
                           fontsize = 1) +
  theme(plot.title = element_text(face = "bold", size = 10))
train_HRforest
ggsave(width=10,height=6,'04.multiCox_step.pdf',train_HRforest)
ggsave(width=10,height=6,'04.multiCox_step.png',train_HRforest)



# length(cox_gene)

# Calculate the risk score -------------
cox_gene_rt <- train_data[ , c('OS.time','OS', cox_gene)]  # The screening results cox_gene stepwise regression mul_cox_gene multivariate

riskScore <- predict(step_Cox, type="lp", newdata = cox_gene_rt) ###for stepwise regression

# riskScore <- predict(cox_more_2,type="lp",newdata=cox_gene_rt)  # for Multivariate cox
riskScore <- data.frame(riskScore)

write.csv(riskScore,file = "risk.csv")
saveRDS(cox_gene,"cox_gene.rds")
saveRDS(step_Cox,"step_Cox.rds")
saveRDS(riskScore, file="riskScore.rds")


identical(rownames(riskScore), rownames(cox_gene_rt))
cox_gene_rt$riskscore <- riskScore[,1]
cox_gene_rt$risk <- ifelse(cox_gene_rt$riskscore > median(cox_gene_rt$riskscore), "High Risk", "Low Risk")
dim(cox_gene_rt)

# K-M curve ----------
kmfit <- survfit(Surv(OS.time, OS) ~ risk, data = cox_gene_rt) 
surv_pvalue(kmfit)   

pdf('04.KM.pdf',w=7,h=7,onefile=F,family='Times')
s_surv <- ggsurvplot(kmfit,
                     pval = TRUE,   #!!
                     conf.int = F,
                     legend.labs=c("High Risk", "Low Risk"),
                     legend.title="Risk",
                     title="TCGA-CRC",
                     risk.table = TRUE,
                     risk.table.col = "strata",
                     #linetype = "strata",
                     # surv.median.line = "hv",
                     ggtheme = theme_bw(),
                     palette=c("#dc143c","#6495ed"))

s_surv$table <- s_surv$table +
  labs(x = "Time (days)") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 18, face = "bold", family = "Times"),
        axis.title.x = element_text(size = 18, face = "bold", family = "Times"),
        axis.text = element_text(size = 16, face = "bold", family = "Times"),
        plot.title = element_text(size = 20, face = "bold", family = "Times"),
        text = element_text(size = 18, face = "bold", family = "Times"))
s_surv$plot <- s_surv$plot +
  labs(x = "Time (days)") +
  labs(y = "CRC survival probability") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1,1),
        legend.position = c(1,1),
        legend.margin =  ggplot2::margin(c(-5,4,4,3)),
        legend.background = element_rect(fill = "white", color = "black", size = 0.2),
        axis.title.y = element_text(size = 16, face = "bold", family = "Times"),
        axis.title.x = element_text(size = 16, face = "bold", family = "Times"),
        axis.text = element_text(size = 14, face = "bold", family = "Times"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, family = "Times"),
        legend.text = element_text(size = 16, face = "bold", family = "Times"),
        text = element_text(size = 20, face = "bold", family = "Times"))
# s_surv[["plot"]][["layers"]][[4]][["aes_params"]][["label"]] <- as.expression(bquote(italic('p')==.(round(surv_pvalue(kmfit, method =  "survdiff")[,2],10))))
print(s_surv)
dev.off()

png('04.KM.png',w=7,h=7,family='Times',units='in',res=600)
s_surv <- ggsurvplot(kmfit,
                     pval = TRUE,   #!!
                     conf.int = F,
                     legend.labs=c("High Risk", "Low Risk"),
                     legend.title="Risk",
                     title="TCGA-CRC",
                     risk.table = TRUE,
                     risk.table.col = "strata",
                     #linetype = "strata",
                     # surv.median.line = "hv",
                     ggtheme = theme_bw(),
                     palette=c("#dc143c","#6495ed"))

s_surv$table <- s_surv$table +
  labs(x = "Time (days)") +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(size = 18, face = "bold", family = "Times"),
        axis.title.x = element_text(size = 18, face = "bold", family = "Times"),
        axis.text = element_text(size = 16, face = "bold", family = "Times"),
        plot.title = element_text(size = 20, face = "bold", family = "Times"),
        text = element_text(size = 18, face = "bold", family = "Times"))
s_surv$plot <- s_surv$plot +
  labs(x = "Time (days)") +
  labs(y = "CRC survival probability") +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.justification = c(1,1),
        legend.position = c(1,1),
        legend.margin =  ggplot2::margin(c(-5,4,4,3)),
        legend.background = element_rect(fill = "white", color = "black", size = 0.2),
        axis.title.y = element_text(size = 16, face = "bold", family = "Times"),
        axis.title.x = element_text(size = 16, face = "bold", family = "Times"),
        axis.text = element_text(size = 14, face = "bold", family = "Times"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5, family = "Times"),
        legend.text = element_text(size = 16, face = "bold", family = "Times"),
        text = element_text(size = 20, face = "bold", family = "Times"))
# s_surv[["plot"]][["layers"]][[4]][["aes_params"]][["label"]] <- as.expression(bquote(italic('p')==.(round(surv_pvalue(kmfit, method =  "survdiff")[,2],10))))
print(s_surv)
dev.off()


# 4 ROC ----------
library(survivalROC)

risk <- cox_gene_rt
ROC <- risk
cutoff_1<- 365*1
cutoff_3 <- 365*3
cutoff_5 <- 365*5
# cutoff_7 <- 365*7

year_1= survivalROC(Stime=ROC$OS.time,
                    status=ROC$OS,
                    marker = ROC$riskscore,
                    predict.time = cutoff_1,
                    method = 'KM')

year_3= survivalROC(Stime=ROC$OS.time,
                    status=ROC$OS,
                    marker = ROC$riskscore,
                    predict.time = cutoff_3,
                    method = 'KM')

year_5= survivalROC(Stime=ROC$OS.time,
                    status=ROC$OS,
                    marker = ROC$riskscore,
                    predict.time = cutoff_5,
                    method = 'KM')

# year_7 = survivalROC(Stime=ROC$OS.time,
#                      status=ROC$OS,
#                      marker = ROC$riskscore,
#                      predict.time = cutoff_7,
#                      method = 'KM')

pdf('05.ROC.pdf',w=10,h=9)
par(pin = c(4,4), mar = c(5,6,6,1),family="serif")      ##par(family="serif")  :all
# pin: Graphic dimensions (width and height) expressed in inches
# mai: Boundary size expressed as a numerical vector, in the order of "bottom, left, top, right", with the unit being inches
# mar: The size of the boundary expressed as a numerical vector, in the order of "bottom, left, top, right", with the unit being British minutes. Default: (5,4,4,2) + 0.1
plot(year_1$FP, year_1$TP,
     type="l",col="red",xlim=c(0,1), ylim=c(0,1),
     xlab="FP",
     ylab="TP",
     main="TCGA-CRC, Method=KM\n Year = 1,3,5",
     cex.axis=1.8,
     cex.lab=2.0,
     cex.main=2.0,
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub =2,
     lwd=2)

abline(0,1,col="gray",lty=2)

lines(year_3$FP, year_3$TP, type="l",col="green",xlim=c(0,1), ylim=c(0,1,lwd=2))

lines(year_5$FP, year_5$TP, type="l",col="blue",xlim=c(0,1), ylim=c(0,1),lwd=2)

legend(0.55,0.2,c(paste("AUC of 1 year =",round(year_1$AUC,2)),
                  paste("AUC of 3 year =",round(year_3$AUC,2)),
                  paste("AUC of 5 year =",round(year_5$AUC,2))),
       x.intersp=1, y.intersp=0.8,
       lty= 1 ,lwd= 2,col=c("red","green",'blue'),
       bty = "n",#
       seg.len=1,cex=1.5,
       text.font = 2)#
dev.off()


png('05.ROC.png',w=10,h=9,units='in',res=600)
par(pin = c(4,4), mar = c(5,6,6,1),family="serif")      ##par(family="serif")  :all
# pin: Graphic dimensions (width and height) expressed in inches
# mai: Boundary size expressed as a numerical vector, in the order of "bottom, left, top, right", with the unit being inches
# mar: The size of the boundary expressed as a numerical vector, in the order of "bottom, left, top, right", with the unit being British minutes. Default: (5,4,4,2) + 0.1
plot(year_1$FP, year_1$TP,
     type="l",col="red",xlim=c(0,1), ylim=c(0,1),
     xlab="FP",
     ylab="TP",
     main="TCGA-CRC, Method=KM\n Year = 1,3,5",
     cex.axis=1.8,
     cex.lab=2.0,
     cex.main=2.0,
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub =2,
     lwd=2)

abline(0,1,col="gray",lty=2)

lines(year_3$FP, year_3$TP, type="l",col="green",xlim=c(0,1), ylim=c(0,1,lwd=2))

lines(year_5$FP, year_5$TP, type="l",col="blue",xlim=c(0,1), ylim=c(0,1),lwd=2)

legend(0.55,0.2,c(paste("AUC of 1 year =",round(year_1$AUC,2)),
                  paste("AUC of 3 year =",round(year_3$AUC,2)),
                  paste("AUC of 5 year =",round(year_5$AUC,2))),
       x.intersp=1, y.intersp=0.8,
       lty= 1 ,lwd= 2,col=c("red","green",'blue'),
       bty = "n",#
       seg.len=1,cex=1.5,
       text.font = 2)#
dev.off()


# 5 Survival curve ----------
rt<-risk
rt=rt[order(rt$riskscore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="Low Risk"])
highLength=length(riskClass[riskClass=="High Risk"])
line=rt[,"riskscore"]
line[line>10]=10
pdf(file="06.riskScore.pdf",width = 10,height = 6)
par(pin = c(4,4), mar = c(6,6,6,1),family='serif')
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("#6495ed",lowLength),
           rep("#dc143c",highLength)),
     cex.axis=1.8,  ##The scaling factor of the scale text on the coordinate axes. Similar to cex.
     cex.lab=2.0, the scaling factor of the scale text on the coordinate axes. Similar to cex.
     cex.main=2.0, ## The scaling factor of the title.
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub =2)
legend("topleft",c(paste("High Risk"),
                   paste("Low Risk")),
       x.intersp=1, y.intersp=0.8,
       pch =20,col=c("#dc143c","#6495ed"),
       bty = 1,# The type of bty box
       seg.len=1,cex=1.5,
       text.font = 2)#
abline(h=median(rt$riskScore),v=lowLength,lty=2)
# title(main = "TCGA-AML",cex.main=2)
dev.off()

png(file="06.riskScore.png",width = 10,height = 6,units='in',res=600 )
par(pin = c(4,4), mar = c(6,6,6,1),family='serif')
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("#6495ed",lowLength),
           rep("#dc143c",highLength)),
     cex.axis=1.8,  ##The scaling factor of the scale text on the coordinate axes. Similar to cex.
     cex.lab=2.0, the scaling factor of the scale text on the coordinate axes. Similar to cex.
     cex.main=2.0, ## The scaling factor of the title.
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub =2)
legend("topleft",c(paste("High Risk"),
                   paste("Low Risk")),
       x.intersp=1, y.intersp=0.8,
       pch =20,col=c("#dc143c","#6495ed"),
       bty = 1,# The type of bty box
       seg.len=1,cex=1.5,
       text.font = 2)#
abline(h=median(rt$riskScore),v=lowLength,lty=2)
# title(main = "TCGA-AML",cex.main=2)
dev.off()


# 6 Survival state ----------
rt<-risk
rt <- rt[order(rt$riskscore),]
riskClass <- rt[,"risk"]
lowLength <- length(riskClass[riskClass=="Low Risk"])
highLength <- length(riskClass[riskClass=="High Risk"])
color=as.vector(rt$OS)
color[color==1]="#dc143c"     ##1:dead, 0:alive
  color[color==0]="#6495ed"
    pdf(file="07.survStat.pdf",width = 12,height = 6)
    par(pin = c(4,4), mar = c(6,6,6,1),family='serif')
    # pin: Graphic dimensions (width and height) expressed in inches
    # mai: Boundary size expressed as a numerical vector, in the order of "bottom, left, top, right", with the unit being inches
    # mar: The size of the boundary expressed as a numerical vector, in the order of "bottom, left, top, right", with the unit being British minutes. Default: (5,4,4,2) + 0.1
    plot(rt$OS.time /365,
         pch=19,
         xlab="Patients (increasing risk socre)",
         ylab="Time (years)",
         col=color,
         cex.axis=1.8,  ##The scaling factor of the scale text on the coordinate axes. Similar to cex.
         cex.lab=2.0,   ##The scaling factor of the scale text on the coordinate axes. Similar to cex.
         cex.main=2.0,   ##The scaling factor of the title.
         font.axis = 2,
         font.lab = 2,
         font.main = 2,
         font.sub =2,
         #main="\n Risk Score Distribution"
    )
    abline(v=lowLength,lty=2)
    legend("topright",c(paste("Dead"),
                        paste("Alive")),
           x.intersp=1, y.intersp=0.8,
           pch =20,col=c("#dc143c","#6495ed"),
           bty = "1",# The types of bty frames: n: Borderless 1: With borders
           seg.len=1,cex=1.5,
           text.font = 2)
    # title(main = "TCGA-AML",cex.main=2)
    dev.off()
    
    png(file="07.survStat.png",width = 10,height = 6,units='in',res=600 )
    par(pin = c(4,4), mar = c(6,6,6,1),family='serif')
    # pin: Graphic dimensions expressed in inches (width and height)
# mai: The size of the boundary expressed as a numerical vector, in the order of "bottom, left, top, right", with the unit being inches
# mar: The size of the boundary expressed as a numerical vector, in the order of "bottom, left, top, right", with the unit being British minutes. Default: (5,4,4,2) + 0.1
    plot(rt$OS.time/365,
         pch=19,
         xlab="Patients (increasing risk socre)",
         ylab="Time (years)",
         col=color,
         cex.axis=1.8,  ##The scaling factor of the scale text on the coordinate axes. Similar to cex.
         cex.lab=2.0,   ##The scaling factor of the scale text on the coordinate axes. Similar to cex.
         cex.main=2.0,   ##The scaling factor of the title.
         font.axis = 2,
         font.lab = 2,
         font.main = 2,
         font.sub =2,
         #main="\n Risk Score Distribution"
    )
    abline(v=lowLength,lty=2)
    legend("topright",c(paste("Dead"),
                        paste("Alive")),
           x.intersp=1, y.intersp=0.8,
           pch =20,col=c("#dc143c","#6495ed"),
           bty = "1",# The types of bty frames: n: Borderless 1: With borders
           seg.len=1,cex=1.5,
           text.font = 2)
    # title(main = "TCGA-AML",cex.main=2)
    dev.off()
    
    
    
    
    
    



