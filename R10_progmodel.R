rm(list = ls())

setwd("/data/nas1/yuanyt/project/05.YQQL-0602-6/")
if (! dir.exists("./10_progmodel")){
  dir.create("./10_progmodel")
}
setwd("./10_progmodel")

## Univariate Cox----------
# Select the required clinical data
survival <- read.csv("../02_ssGSEA/icd.group.csv",header = T,row.names = 1)
survival <- survival[,c(1,2,3)]
riskScore <- readRDS("../08_LASSO/riskScore.rds")
riskScore$sample <- rownames(riskScore)
survival <- merge(riskScore,survival,by= "sample")


phenotype_coad <- read_tsv(file = 'TCGA-COAD.GDC_phenotype.tsv.gz')
phenotype_read <- read_tsv(file = 'TCGA-READ.GDC_phenotype.tsv.gz')
coad <- phenotype_coad[,colnames(phenotype_coad) %in% colnames(phenotype_read)]
read <- phenotype_read[,colnames(phenotype_read) %in% colnames(phenotype_coad)]


phenotype <- rbind(coad,read)

train_phenotype <- data.frame(sample = phenotype$submitter_id.samples,
                              age = phenotype$age_at_initial_pathologic_diagnosis,
                              gender = phenotype$gender.demographic,
                              race = phenotype$race.demographic,
                              # TNM.stage = phenotype$tumor_stage.diagnoses,
                              T.stage = phenotype$pathologic_T,
                              N.stage = phenotype$pathologic_N,
                              M.stage = phenotype$pathologic_M)


#Retain the training set samples, list them as samples, behavioral survival time, status and various clinical manifestations
train_phenotype <- merge(train_phenotype, survival, by='sample')
write.csv(train_phenotype, file = '01.phenotype.csv', row.names = F, quote = F)

#Convert survival and various clinical manifestations into numerical values
train_phenotype2 <- train_phenotype
train_phenotype2$OS <- as.numeric(train_phenotype2$OS)
train_phenotype2$OS.time <- as.numeric(train_phenotype2$OS.time)

table(train_phenotype2$age)
train_phenotype2$age <- ifelse(train_phenotype2$age>60,1,0)

table(train_phenotype2$gender)
train_phenotype2$gender <- ifelse(train_phenotype2$gender == 'female',1,0)

table(train_phenotype2$race)
train_phenotype2$race <- gsub('american indian or alaska native',NA,train_phenotype2$race)
train_phenotype2$race <- gsub('asian','1',train_phenotype2$race)
train_phenotype2$race <- gsub('black or african american','2',train_phenotype2$race)
train_phenotype2$race <- gsub('white','3',train_phenotype2$race)
train_phenotype2$race <- gsub('not reported',NA,train_phenotype2$race)


table(train_phenotype2$T.stage)
train_phenotype2$T.stage <- gsub('a','',train_phenotype2$T.stage)
train_phenotype2$T.stage <- gsub('b','',train_phenotype2$T.stage)
#train_phenotype2$T.stage<-gsub('c','',train_phenotype2$T.stage)
#train_phenotype2$T.stage<-gsub('d','',train_phenotype2$T.stage)
train_phenotype2$T.stage <- gsub('Tis',NA,train_phenotype2$T.stage)
train_phenotype2$T.stage <- gsub('T3','3',train_phenotype2$T.stage)
train_phenotype2$T.stage <- gsub('T4','4',train_phenotype2$T.stage)
train_phenotype2$T.stage <- gsub('T2','2',train_phenotype2$T.stage)
train_phenotype2$T.stage <- gsub('T1','1',train_phenotype2$T.stage)

table(train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('a','',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('b','',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('c','',train_phenotype2$N.stage)
# train_phenotype2$N.stage<-gsub(' (i-)','',train_phenotype2$N.stage,fixed = T)
# train_phenotype2$N.stage<-gsub(' (i+)','',train_phenotype2$N.stage,fixed = T)
# train_phenotype2$N.stage<-gsub(' (mol+)','',train_phenotype2$N.stage,fixed = T)
# train_phenotype2$N.stage<-gsub('mi','',train_phenotype2$N.stage,fixed = T)
train_phenotype2$N.stage<-gsub('NX',NA,train_phenotype2$N.stage,fixed = T)
train_phenotype2$N.stage<-gsub('N2','2',train_phenotype2$N.stage)
# train_phenotype2$N.stage<-gsub('N3','3',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('N1','1',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('N0','0',train_phenotype2$N.stage)

table(train_phenotype2$M.stage)
# train_phenotype2$M.stage<-gsub('cM0 (i+)','0',train_phenotype2$M.stage,fixed = T)
train_phenotype2$M.stage<-gsub('a','',train_phenotype2$M.stage,fixed = T)
train_phenotype2$M.stage<-gsub('b','',train_phenotype2$M.stage,fixed = T)
train_phenotype2$M.stage<-gsub('MX',NA,train_phenotype2$M.stage,fixed = T)
train_phenotype2$M.stage<-gsub('M0','0',train_phenotype2$M.stage,fixed = T)
train_phenotype2$M.stage<-gsub('M1','1',train_phenotype2$M.stage,fixed = T)

# table(train_phenotype2$grade)
# train_phenotype2$grade <- gsub('G3','3',train_phenotype2$grade)
# train_phenotype2$grade <- gsub('G4','4',train_phenotype2$grade)
# train_phenotype2$grade <- gsub('G2','2',train_phenotype2$grade)
# train_phenotype2$grade <- gsub('G1','1',train_phenotype2$grade)



os_risk_clinical <- train_phenotype2 %>% 
  column_to_rownames(var = 'sample')
dim(os_risk_clinical)

#Extract covariates
colnames_os_risk_clinical <- colnames(os_risk_clinical)
covariates <- colnames_os_risk_clinical[-which(colnames_os_risk_clinical %in% c("OS", "OS.time"))]

os_risk_clinical$age <- factor(os_risk_clinical$age)
os_risk_clinical$gender <- factor(os_risk_clinical$gender)
os_risk_clinical$race <- factor(os_risk_clinical$race)
# os_risk_clinical$TNM.stage <- factor(os_risk_clinical$TNM.stage)
os_risk_clinical$T.stage <- factor(os_risk_clinical$T.stage)
os_risk_clinical$N.stage <- factor(os_risk_clinical$N.stage)
os_risk_clinical$M.stage <- factor(os_risk_clinical$M.stage)
# os_risk_clinical$grade <- factor(os_risk_clinical$grade)

#Univariate cox
library(survival)
res.risk = coxph(Surv(time = OS.time, event = OS) ~ riskScore, data = os_risk_clinical) %>% summary
res.risk = c(res.risk$conf.int[-2], res.risk$coefficients[5])

res.age = coxph(Surv(time = OS.time, event = OS) ~ age, data = os_risk_clinical) %>% summary
res.age = c(res.age$conf.int[-2], res.age$coefficients[5])

res.gender = coxph(Surv(time = OS.time, event = OS) ~ gender, data = os_risk_clinical) %>% summary
res.gender = c(res.gender$conf.int[-2], res.gender$coefficients[5])

res.race = coxph(Surv(time = OS.time, event = OS) ~ race, data = os_risk_clinical) %>% summary
res.race = cbind(res.race$conf.int[,-2], res.race$coefficients[,5])

# res.TNM.stage = coxph(Surv(time = OS.time, event = OS) ~ TNM.stage, data = os_risk_clinical) %>% summary
# res.TNM.stage = cbind(res.TNM.stage$conf.int[,-2], res.TNM.stage$coefficients[,5])

res.T.stage = coxph(Surv(time = OS.time, event = OS) ~ T.stage, data = os_risk_clinical) %>% summary
res.T.stage = cbind(res.T.stage$conf.int[,-2], res.T.stage$coefficients[,5])

res.N.stage = coxph(Surv(time = OS.time, event = OS) ~ N.stage, data = os_risk_clinical) %>% summary
res.N.stage = cbind(res.N.stage$conf.int[,-2], res.N.stage$coefficients[,5])

res.M.stage = coxph(Surv(time = OS.time, event = OS) ~ M.stage, data = os_risk_clinical) %>% summary
res.M.stage = c(res.M.stage$conf.int[-2], res.M.stage$coefficients[5])

# res.grade = coxph(Surv(time = OS.time, event = OS) ~ grade, data = os_risk_clinical) %>% summary
# res.grade = cbind(res.grade$conf.int[,-2], res.grade$coefficients[,5])

#Organize the single-factor table
res.ref = c(1,1,1,NA)

res = rbind(res.risk, res.age, res.gender, res.ref, res.race, res.ref, res.T.stage, res.ref, res.N.stage, res.M.stage) %>% as.data.frame()
rownames(res)
res$Indicators = c('riskScore',
                   'age',
                   'gender', 
                   'asian(Reference)', 'black or african american', 'white',
                   'T1(Reference)','T2','T3','T4',
                   'N0(Reference)','N1','N2',
                   'M1 vs M0')

colnames(res) = c("hr","low","up","pv","Indicator")
# res$p = signif(res$pv, 2) %>% paste0("p = ", .) #Keep two significant figures and add "p =" before the number
# res$p[is.na(res$pv)] = NA #is.na is used to detect A problem with a value of #N/A and returns TRUE or FALSE
# res$Indicator = factor(res$Indicator, levels = rev(res$Indicator)) #rev reverses the order of the elements in a vector or matrix
# rownames(res) <- res$Indicator
res2 <- data.frame(p.value=res$pv,
                   HR=res$hr,
                   HR.95L=res$low,
                   HR.95H=res$up,
                   Indicator=res$Indicator)
rownames(res2) <- res2$Indicator
res2 <- subset(res2, select = -c(Indicator))
write.csv(res2, file = "02.univariate_cox_prog_forest.csv", quote = F, row.names = T)

#Organize drawing data
library(tidyr)
#Organize it into HR(HR.95L-HR.95H)
hz <- paste(round(res2$HR,3), #Keep to three decimal places
            "(",
            round(res2$HR.95L,3),
            "-",
            round(res2$HR.95H,3),
            ")",
            sep = "") 
hz

hz[c(4,7,11)] <- "" #The reference value is empty

tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.05,"< 0.05",round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1

#draw
library(forestplot)
pdf(file = '01.uncoxforest.pdf',height = 6, width = 12, onefile = F,family='Times')
forestplot(labeltext=tabletext, 
           graph.pos=4,  #It is the position where the Pvalue box plot is located
           is.summary = c(TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,rep(FALSE, 70)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,res2$HR),
           lower=c(NA,res2$HR.95L), #Lower limit of the 95% confidence interval
           upper=c(NA,res2$HR.95H), #Higher limit of the 95% confidence interval
           boxsize=0.1,lwd.ci=3,   #The size of the box and the width of the line
           ci.vertices.height = 0.08,ci.vertices=TRUE, #Confidence intervals are defined by line width, height and shape
           zero=1,lwd.zero=0.5,    #The position of the zero line width reference line
           colgap=unit(5,"mm"),    #Column gap
           xticks = c(0,1,5,10,30), #Horizontal coordinate scale
           lwd.xaxis=2,            #Width of the X-axis
           lineheight = unit(1.2,"cm"), #Fixed row height
           graphwidth = unit(.6,"npc"), #The width ratio of the figure in the table
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #Error bar display mode
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #A black line is added at the top of the third line, and the numbers within the quotation marks indicate the position of the line
           #                 "19" = gpar(lwd=2, col="black")),#At the bottom of the last line, add a black line. The number in "16" is nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#Graphic margins
           #The cex parameter in the fpTxtGp function sets the size of each component
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Univariate Cox") # The grid lines perpendicular to the X-axis correspond to each scale
dev.off()

png(file = '01.uncoxforest.png',height = 6, width = 12, units = 'in',res = 600,family='Times')
forestplot(labeltext=tabletext, 
           graph.pos=4,  #It is the position where the Pvalue box plot is located
           is.summary = c(TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,rep(FALSE, 7)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,res2$HR),
           lower=c(NA,res2$HR.95L), #Lower limit of the 95% confidence interval
           upper=c(NA,res2$HR.95H), #Higher limit of the 95% confidence interval
           boxsize=0.1,lwd.ci=3,   #The size of the box and the width of the line
           ci.vertices.height = 0.08,ci.vertices=TRUE, #Confidence intervals are defined by line width, height and shape
           zero=1,lwd.zero=0.5,    #The position of the zero line width reference line
           colgap=unit(5,"mm"),    #Column gap
           xticks = c(0,1,5,10,30), #Horizontal coordinate scale
           lwd.xaxis=2,            #Width of the X-axis
           lineheight = unit(1.2,"cm"), #Fixed row height
           graphwidth = unit(.6,"npc"), #The width ratio of the figure in the table
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #Error bar display mode
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
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
           grid = T,
           title = "Univariate Cox") # The grid lines perpendicular to the X-axis correspond to each scale
dev.off()

## 07-2 Multivariate Cox----------
cox_more <- coxph(Surv(time = OS.time, event = OS) ~ riskScore + age + T.stage + N.stage + M.stage, data = os_risk_clinical)
cox_zph <- cox.zph(cox_more)
cox_table <- cox_zph$table[-nrow(cox_zph$table),]      #PH hypothesis test

res.mul <- coxph(Surv(time = OS.time, event = OS) ~ riskScore + age + T.stage + N.stage + M.stage, data = os_risk_clinical) %>% summary
res.mul  <-  cbind(res.mul$conf.int[,-2], res.mul$coefficients[,5]) %>% as.data.frame()
rownames(res.mul)

res.mul = rbind(res.mul[c(1:2),], res.ref, res.mul[c(3:5),], res.ref, res.mul[c(6:8),])
res.mul$Indicators = c("riskScore","age","T1(Reference)",'T2','T3','T4','N0(Reference)','N1','N2',"M1 vs.M0")
colnames(res.mul) = c("hr","low","up","pv","Indicator")
# res.mul$p = signif(res.mul$pv, 2) %>% paste0("p = ", .)
# res.mul$p[is.na(res.mul$pv)] = NA
# res.mul$Indicator = factor(res.mul$Indicator, levels = rev(res.mul$Indicator))
# rownames(res.mul) <- res.mul$Indicator

multi_res <- data.frame(p.value=res.mul$pv,
                        HR=res.mul$hr,
                        HR.95L=res.mul$low,
                        HR.95H=res.mul$up,
                        Indicator=res.mul$Indicator)
rownames(multi_res) <- multi_res$Indicator
multi_res
multi_res <- subset(multi_res, select = -c(Indicator))
write.csv(multi_res,
          file = "03.multivariate_cox_prog_result.csv",
          quote = F,
          row.names = T)

library(tidyr)
hz <- paste(round(multi_res$HR,3),
            "(",round(multi_res$HR.95L,3),
            "-",round(multi_res$HR.95H,3),")",sep = "")
hz

hz[c(3,7)] <- ""
tabletext <- cbind(c(NA,rownames(multi_res)),
                   c("P value",ifelse(multi_res$p.value<0.05,
                                      "< 0.05",
                                      round(multi_res$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
pdf(file = '02.mulcoxforest.pdf',height = 6, width = 12, onefile = F,family='Times')
forestplot(labeltext=tabletext, 
           graph.pos=4,  #It is the position where the Pvalue box plot is located
           is.summary = c(TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,multi_res$HR),
           lower=c(NA,multi_res$HR.95L), #Lower limit of the 95% confidence interval
           upper=c(NA,multi_res$HR.95H), #Higher limit of the 95% confidence interval
           boxsize=0.1,lwd.ci=3,   #The size of the box and the width of the line
           ci.vertices.height = 0.08,ci.vertices=TRUE, #Confidence intervals are defined by line width, height and shape
           zero=1,lwd.zero=0.5,    #The position of the zero line width reference line
           colgap=unit(5,"mm"),    #Column gap
           xticks = c(0,1,5,10,15,20,25), #Horizontal coordinate scale
           lwd.xaxis=2,            #Width of the X-axis
           lineheight = unit(1.2,"cm"), #Fixed row height
           graphwidth = unit(.6,"npc"), #The width ratio of the figure in the table
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #Error bar display mode
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
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
           grid = T,
           title = "Multivariate Cox") # The grid lines perpendicular to the X-axis correspond to each scale
dev.off()

png(file = '02.mulcoxforest.png',height = 6, width = 12, units = 'in',res = 600,family='Times')
forestplot(labeltext=tabletext, 
           graph.pos=4,  #It is the position where the Pvalue box plot is located
           is.summary = c(TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,multi_res$HR),
           lower=c(NA,multi_res$HR.95L), #Lower limit of the 95% confidence interval
           upper=c(NA,multi_res$HR.95H), #Higher limit of the 95% confidence interval
           boxsize=0.1,lwd.ci=3,   #The size of the box and the width of the line
           ci.vertices.height = 0.08,ci.vertices=TRUE, #Confidence intervals are defined by line width, height and shape
           zero=1,lwd.zero=0.5,    #The position of the zero line width reference line
           colgap=unit(5,"mm"),    #Column gap
           xticks = c(0,1,5,10,15,20,25), #Horizontal coordinate scale
           lwd.xaxis=2,            #Width of the X-axis
           lineheight = unit(1.2,"cm"), #Fixed row height
           graphwidth = unit(.6,"npc"), #The width ratio of the figure in the table
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #Error bar display mode
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
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
           grid = T,
           title = "Multivariate Cox") # The grid lines perpendicular to the X-axis correspond to each scale
dev.off()


##07-3 Construct the COX model and draw the nomogram---------
multi_cov<-c('riskScore',"age","N.stage","M.stage")
cox_data_prog <- as.formula(paste0('Surv(OS.time, OS)~',
                                   paste(multi_cov,
                                         sep = '',
                                         collapse = '+')))
cox_more_prog <- coxph(cox_data_prog,
                       data = as.data.frame(os_risk_clinical))

# Nomogram
library(rms)
ddist <- datadist(os_risk_clinical) #The function for generating data summary
options(datadist='ddist')

# Draw a nomogram
# os_risk_clinical$Age <- ifelse(os_risk_clinical$Age==1,'>60','<=60')
res.cox <- psm(cox_data_prog, data = os_risk_clinical, dist = 'lognormal')
surv <- Survival(res.cox) # Construct a survival probability function
function(x) surv(365, x) # The probability of an event occurring within one year
#function(x) surv(730, x) #The probability of the event occurring within two years
function(x) surv(1095, x) # The probability of the event occurring within three years
function(x) surv(1825, x) # The probability of the event occurring within five years
# function(x) surv(2555, x) # The probability of a 7-year event occurring
nom.cox <- nomogram(res.cox,
                    fun = list(function(x) surv(365, x),
                               function(x) surv(1095, x),
                               function(x) surv(1825, x)),
                    funlabel=c("1-year Survival Probability", "3-year Survival Probability", "5-year Survival Probability"),
                    maxscale = 10,
                    fun.at = c(0.01,seq(0.1,0.9,by=0.2),0.95,0.99),
                    lp=F)

library(regplot)

regplot(res.cox, plots=c("violin","bars"), observation=T, title="Survival Nomogram", failtime=c(365*1,365*3,365*5), prfail=T,  points=TRUE)

png("03.nomogram_classical.png",width = 12,height = 6,family='Times',units='in',res=600)
plot(nom.cox, lplabel="Linear Predictor",
     xfrac=.35,varname.label=TRUE, 
     varname.label.sep="=", ia.space=.2, 
     tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE, 
     cap.labels=FALSE,cex.var = 1.35,cex.axis = 1.05,lwd=5,
     label.every = 1,col.grid = gray(c(0.8, 0.95)))

dev.off()

pdf("03.nomogram_classical.pdf",width = 12,height = 6,family='Times')
plot(nom.cox, lplabel="Linear Predictor",
     xfrac=.35,varname.label=TRUE, varname.label.sep="=", ia.space=.2, 
     tck=NA, tcl=-0.20, lmgp=0.3,
     points.label='Points', total.points.label='Total Points',
     total.sep.page=FALSE, 
     cap.labels=FALSE,cex.var = 1.35,cex.axis = 1.05,lwd=5,
     label.every = 1,col.grid = gray(c(0.8, 0.95)))

dev.off()



##07-4 Construct the calibration curve---------
coxm_1 <- cph(cox_data_prog,
              data=os_risk_clinical,
              surv=T,
              x=T,
              y=T,
              time.inc = 365)
cal_1<-calibrate(coxm_1,u=365,cmethod='KM',m=150,B=1000)

coxm_3 <- cph(cox_data_prog,
              data=os_risk_clinical,
              surv=T,
              x=T,
              y=T,
              time.inc = 365*3)
cal_3<-calibrate(coxm_3,u=365*3,cmethod='KM',m=150,B=1000)

##Draw the 3-year survival calibration curve
##If time.in and u are the same, they both refer to the time nodes to be evaluated
coxm_5 <- cph(cox_data_prog,
              data=os_risk_clinical,
              surv=T,
              x=T,y=T,
              time.inc = 5*365)
cal_5 <-calibrate(coxm_5,u=5*365,cmethod='KM',m=150,B=1000)

# coxm_7 <- cph(cox_data_prog,
#               data=os_risk_clinical,
#               surv=T,
#               x=T,y=T,
#               time.inc = 7*365)
# cal_7 <-calibrate(coxm_7,u=7*365,cmethod='KM',m=150,B=1000)


pdf(file = '04.calibrate.pdf',w=9,h=8,family='Times')
par(mar=c(7,4,4,3),cex=1.5)

plot(cal_1,
     subtitles = F,
     lwd=2,lty=1, ##Set the shape and size of the lines
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##Set a color
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability',#Sticky note
     ylab='Actual 1-5 year Survival Probability',#Tag
     col="#00468b",#Set a color
     xlim = c(0,1),ylim = c(0,1)) ##The X-axis and Y-axis range
plot(cal_3,
     add = T,
     subtitles = F,
     lwd=2,lty=1, ##Set the shape and size of the lines
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##Set a color
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability',#Sticky note
     ylab='Actual 1-5 year Survival Probability',#Tag
     col="#ed0000",#Set a color
     xlim = c(0,1),ylim = c(0,1)) ##The X-axis and Y-axis range
plot(cal_5,
     add = T,
     subtitles = F,
     lwd=2,lty=1,  ##Set the line width and line type
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##Set a color
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability',#Sticky note
     ylab='Actual 1-5 year Survival Probability',#Tag
     col="#42b540",#Set a color
     xlim = c(0,1),ylim = c(0,1)) ##The X-axis and Y-axis range
# plot(cal_7,
#      add = T,
#      subtitles = F,
#      lwd=2,lty=1, ##Set the shape and size of the lines
#      errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##Set a color
#      xlab='Nomogram-Predicted Probability of 3-7 year Survival Probability',#Sticky note
#      ylab='Actual 3-7 year Survival Probability',#Tag
#      col="#42b540",#Set a color
#      xlim = c(0,1),ylim = c(0,1)) ##The X-axis and Y-axis range

#Add a legend
legend("bottomright", legend=c("1-year", "3-year", "5-year"), 
       col=c("#00468b", "#ed0000", "#42b540"), 
       lwd=2)
#Adjust the diagonal
abline(0,1,lty=5,lwd=2,col="grey")
dev.off()

png(file = '04.calibrate.png',w=9,h=8, units = 'in',res = 600,family='Times')
par(mar=c(7,4,4,3),cex=1.5)

plot(cal_1,
     subtitles = F,
     lwd=2,lty=1, ##Set the shape and size of the lines
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##Set a color
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability',#Sticky note
     ylab='Actual 1-5 year Survival Probability',#Tag
     col="#00468b",#Set a color
     xlim = c(0,1),ylim = c(0,1)) ##The X-axis and Y-axis range
plot(cal_3,
     add = T,
     subtitles = F,
     lwd=2,lty=1, ##Set the shape and size of the lines
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##Set a color
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability',#Sticky note
     ylab='Actual 1-5 year Survival Probability',#Tag
     col="#ed0000",#Set a color
     xlim = c(0,1),ylim = c(0,1)) ##The X-axis and Y-axis range
plot(cal_5,
     add = T,
     subtitles = F,
     lwd=2,lty=1,  ##Set the line width and line type
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##Set a color
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability',#Sticky note
     ylab='Actual 1-5 year Survival Probability',#Tag
     col="#42b540",#Set a color
     xlim = c(0,1),ylim = c(0,1)) ##The X-axis and Y-axis range
# plot(cal_7,
#      add = T,
#      subtitles = F,
#      lwd=2,lty=1, ##Set the shape and size of the lines
#      errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##Set a color
#      xlab='Nomogram-Predicted Probability of 3-7 year Survival Probability',#Sticky note
#      ylab='Actual 3-7 year Survival Probability',#Tag
#      col="#42b540",#Set a color
#      xlim = c(0,1),ylim = c(0,1)) ##The X-axis and Y-axis range

#Add a legend
legend("bottomright", legend=c("1-year", "3-year", "5-year"), 
       col=c("#00468b", "#ed0000", "#42b540"), 
       lwd=2)
#Adjust the diagonal
abline(0,1,lty=5,lwd=2,col="grey")
dev.off()



