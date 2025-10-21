rm(list = ls())

setwd("/data/nas1/yuanyt/project/05.YQQL-0602-6/")
if (! dir.exists("./09_Verif")){
  dir.create("./09_Verif")
}
setwd("./09_Verif")


library(lance)
library(readxl)
library(readr)
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(Ipaper)
library(survival)
library(survminer)

###Read data ——————
###Risk coefficient
cox_gene <- readRDS("../08_LASSO/cox_gene.rds")
step_Cox <- readRDS("../08_LASSO/step_Cox.rds")


# 01 Obtain the dataset----------
###GSE17536  GSE17537 GSE103479 GSE28722  GSE29621  GSE72970
survival<-read.csv(file = '../00_Rowdata/survival(GSE17536).csv',header=TRUE, row.names=1,)
data <- read.csv("../00_Rowdata/dat(GSE17536).csv",header=TRUE, row.names=1, check.names=FALSE)
survival_dat<-t(data)
survival_dat<-survival_dat[rownames(survival_dat)%in%survival$sample,]

## Merge survival data
train_dat<-survival_dat[,colnames(survival_dat)%in%cox_gene]
train_dat<-as.data.frame(train_dat)
train_dat$sample<-rownames(train_dat)
train_dat<-merge(survival,train_dat,by='sample')
rownames(train_dat)<-train_dat$sample
train_data<-train_dat[,-1]
train_data <- na.omit(train_data)
colnames(train_data)



#Construction and validation of risk model ------
# Calculate the risk score-------------
cox_gene_rt <- train_data[ , c('OS.time','OS', cox_gene)]  # Screening results

riskScore <- predict(step_Cox, type="lp", newdata = cox_gene_rt) ### for stepwise regression
riskScore <- data.frame(riskScore)

identical(rownames(riskScore), rownames(cox_gene_rt))
cox_gene_rt$riskscore <- riskScore[,1]

####Median
cox_gene_rt$risk <- ifelse(cox_gene_rt$riskscore > median(cox_gene_rt$riskscore), "High Risk", "Low Risk")
dim(cox_gene_rt)

# K-M curve ----------
kmfit <- survfit(Surv(OS.time, OS) ~ risk, data = cox_gene_rt) 
surv_pvalue(kmfit)   

pdf('01.KM.pdf',w=7,h=7,onefile=F,family='Times')
s_surv <- ggsurvplot(kmfit,
                     pval = TRUE,   #!!
                     conf.int = F,
                     legend.labs=c("High Risk", "Low Risk"),
                     legend.title="Risk",
                     title="GEO-CRC",
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


png('01.KM.png',w=7,h=7,family='Times',units='in',res=600)
s_surv <- ggsurvplot(kmfit,
                     pval = TRUE,   #!!
                     conf.int = F,
                     legend.labs=c("High Risk", "Low Risk"),
                     legend.title="Risk",
                     title="GEO-CRC",
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
cutoff_1<- 12*1
cutoff_3 <- 12*3
cutoff_5 <- 12*5
# cutoff_7 <- 12*7

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
#                      method = 'KM')#


pdf('02.ROC.pdf',w=10,h=9)
par(pin = c(4,4), mar = c(5,6,6,1),family="serif")      ##par(family="serif")  :all
# pin: Graphic dimensions expressed in inches (width and height)
# mai: The size of the boundary expressed as a numerical vector, in the order of "bottom, left, top, right", with the unit being inches
# mar: The size of the boundary expressed as a numerical vector, in the order of "bottom, left, top, right", with the unit being British minutes. Default: (5,4,4,2) + 0.1
plot(year_1$FP, year_1$TP,
     type="l",col="red",xlim=c(0,1), ylim=c(0,1),
     xlab="FP",
     ylab="TP",
     main="GEO-CRC, Method=KM\n Year = 1,3,5",
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


png('02.ROC.png',w=10,h=9,units='in',res=600)
par(pin = c(4,4), mar = c(5,6,6,1),family="serif")      ##par(family="serif")  :all
# pin: Graphic dimensions expressed in inches (width and height)
# mai: The size of the boundary expressed as a numerical vector, in the order of "bottom, left, top, right", with the unit being inches
# mar: The size of the boundary expressed as a numerical vector, in the order of "bottom, left, top, right", with the unit being British minutes. Default: (5,4,4,2) + 0.1
plot(year_1$FP, year_1$TP,
     type="l",col="red",xlim=c(0,1), ylim=c(0,1),
     xlab="FP",
     ylab="TP",
     main="GEO-CRC, Method=KM\n Year = 1,3,5",
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


# 5 ----------
rt<-risk
rt=rt[order(rt$riskscore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="Low Risk"])
highLength=length(riskClass[riskClass=="High Risk"])
line=rt[,"riskscore"]
line[line>10]=10
pdf(file="03.riskScore.pdf",width = 10,height = 6)
par(pin = c(4,4), mar = c(6,6,6,1),family='serif')
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("#6495ed",lowLength),
           rep("#dc143c",highLength)),
     cex.axis=1.8,  ##The scaling factor of the scale text on the coordinate axis. Similar to cex.
     cex.lab=2.0,   ##The scaling factor of the scale text on the coordinate axis. Similar to cex.
     cex.main=2.0,   ##The zoom ratio of the title.
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub =2)
legend("topleft",c(paste("High Risk"),
                   paste("Low Risk")),
       x.intersp=1, y.intersp=0.8,
       pch =20,col=c("#dc143c","#6495ed"),
       bty = 1,# The types of bty frames
       seg.len=1,cex=1.5,
       text.font = 2)#
abline(h=median(rt$riskScore),v=lowLength,lty=2)
# title(main = "TCGA-AML",cex.main=2)
dev.off()

png(file="03.riskScore.png",width = 10,height = 6,units='in',res=600 )
par(pin = c(4,4), mar = c(6,6,6,1),family='serif')
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("#6495ed",lowLength),
           rep("#dc143c",highLength)),
     cex.axis=1.8,  ##The scaling factor of the scale text on the coordinate axis. Similar to cex.
     cex.lab=2.0,   ##The scaling factor of the scale text on the coordinate axis. Similar to cex.
     cex.main=2.0,   ##The zoom ratio of the title.
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub =2)
legend("topleft",c(paste("High Risk"),
                   paste("Low Risk")),
       x.intersp=1, y.intersp=0.8,
       pch =20,col=c("#dc143c","#6495ed"),
       bty = 1,# The types of bty frames
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
    pdf(file="04.survStat.pdf",width = 12,height = 6)
    par(pin = c(4,4), mar = c(6,6,6,1),family='serif')
    # pin: Graphic dimensions expressed in inches (width and height)
    # mai: The size of the boundary expressed as a numerical vector, in the order of "bottom, left, top, right", with the unit being inches
    # mar: The size of the boundary expressed as a numerical vector, in the order of "bottom, left, top, right", with the unit being British minutes. Default: (5,4,4,2) + 0.1
    plot(rt$OS.time /12,
         pch=19,
         xlab="Patients (increasing risk socre)",
         ylab="Time (years)",
         col=color,
         cex.axis=1.8,  ##The scaling factor of the scale text on the coordinate axis. Similar to cex.
         cex.lab=2.0,   ##The scaling factor of the scale text on the coordinate axis. Similar to cex.
         cex.main=2.0,   ##The zoom ratio of the title.
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
    
    png(file="04.survStat.png",width = 10,height = 6,units='in',res=600 )
    par(pin = c(4,4), mar = c(6,6,6,1),family='serif')
    # pin: Graphic dimensions expressed in inches (width and height)
    # mai: The size of the boundary expressed as a numerical vector, in the order of "bottom, left, top, right", with the unit being inches
    # mar: The size of the boundary expressed as a numerical vector, in the order of "bottom, left, top, right", with the unit being British minutes. Default: (5,4,4,2) + 0.1
    plot(rt$OS.time/12,
         pch=19,
         xlab="Patients (increasing risk socre)",
         ylab="Time (years)",
         col=color,
         cex.axis=1.8,  ##The scaling factor of the scale text on the coordinate axis. Similar to cex.
         cex.lab=2.0,   ##The scaling factor of the scale text on the coordinate axis. Similar to cex.
         cex.main=2.0,   ##The zoom ratio of the title.
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
    
    


