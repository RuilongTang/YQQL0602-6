rm(list = ls())

setwd("/data/nas1/yuanyt/project/05.YQQL-0602-6/")
if (! dir.exists("./11_Clinical")){
  dir.create("./11_Clinical")
}
setwd("./11_Clinical")


###Read data-------------------------
survival <- read.csv("../02_ssGSEA/icd.group.csv",header = T,row.names = 1)
survival <- survival[,c(1,2,3)]
riskScore <- readRDS("../08_LASSO/riskScore.rds")
riskScore$sample <- rownames(riskScore)
survival <- merge(riskScore,survival,by= "sample")


phenotype_coad <- read_tsv(file = '../10_progmodel/TCGA-COAD.GDC_phenotype.tsv.gz')
phenotype_read <- read_tsv(file = '../10_progmodel/TCGA-READ.GDC_phenotype.tsv.gz')
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


train_phenotype2 <- train_phenotype
train_phenotype2$risk <- ifelse(train_phenotype2$riskScore > median(train_phenotype2$riskScore),"High Risk","Low Risk")
colnames(train_phenotype2)

###age---------------------------------------------
age <- train_phenotype2[,c("age","risk")]
age$age <- ifelse(age$age<60,">=60","<60")

age1 <- xtabs(~age+risk,data=age)
p <- chisq.test(age1)[[3]]

###draw
age2 <- table(age)%>%as.data.frame()
## Stacked graph
library(reshape2)
library(plyr)
library(RColorBrewer)

age3 <- ggplot(age2,aes(x=risk,y=Freq,fill=age))+
  geom_bar(position = "fill",stat="identity",alpha=0.7)+
  theme_bw()+
  scale_fill_manual(values=brewer.pal(3, "Set2"))+ 
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, face = "bold",family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  face = "bold",family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18, color='black', face = "bold",family='Times'),
        title=element_text(size=30, color='black', face = "bold",family='Times'),
        strip.text = element_text(size = 14,family = "Times", face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x="",y="Percentage",fill="Age",title = "Age")
age3
# ggsave(age3,filename = '01.age.pdf',w=7,h=6)
# ggsave(age3,filename = '01.age.png',w=7,h=6)






###gender------------------------------
gender <- train_phenotype2[,c("gender","risk")]

gender1 <- xtabs(~gender+risk,data=gender)
p <- chisq.test(gender1)[[3]]

###draw
gender2 <- table(gender)%>%as.data.frame()
## Stacked graph
library(reshape2)
library(plyr)
library(RColorBrewer)

gender3 <- ggplot(gender2,aes(x=risk,y=Freq,fill=gender))+
  geom_bar(position = "fill",stat="identity",alpha=0.7)+
  theme_bw()+
  scale_fill_manual(values=brewer.pal(3, "Set2"))+ 
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, face = "bold",family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  face = "bold",family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18, color='black', face = "bold",family='Times'),
        title=element_text(size=30, color='black', face = "bold",family='Times'),
        strip.text = element_text(size = 14,family = "Times", face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x="",y="Percentage",fill="Gender",title = "Gender")
gender3
# ggsave(gender3,filename = '02.gender.pdf',w=7,h=6)
# ggsave(gender3,filename = '02.gender.png',w=7,h=6)


###race-------------------
race <- train_phenotype2[,c("race","risk")]
race$race <- gsub('american indian or alaska native',NA,race$race)
race$race <- gsub('not reported',NA,race$race)

race1 <- xtabs(~race+risk,data=race)
p <- chisq.test(race1)[[3]]

###draw
race2 <- table(race)%>%as.data.frame()
## Stacked graph

race3 <- ggplot(race2,aes(x=risk,y=Freq,fill=race))+
  geom_bar(position = "fill",stat="identity",alpha=0.7)+
  theme_bw()+
  scale_fill_manual(values=brewer.pal(3, "Set2"))+ 
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, face = "bold",family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  face = "bold",family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18, color='black', face = "bold",family='Times'),
        title=element_text(size=30, color='black', face = "bold",family='Times'),
        strip.text = element_text(size = 14,family = "Times", face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x="",y="Percentage",fill="race",title = "race")
race3
# ggsave(race3,filename = '03.race.pdf',w=7,h=6)
# ggsave(race3,filename = '03.race.png',w=7,h=6)

###T.stage-------------------
T.stage <- train_phenotype2[,c("T.stage","risk")]
T.stage$T.stage <- gsub('a','',T.stage$T.stage)
T.stage$T.stage <- gsub('b','',T.stage$T.stage)
T.stage$T.stage <- gsub('Tis',NA,T.stage$T.stage)

T.stage1 <- xtabs(~T.stage+risk,data=T.stage)
p <- chisq.test(T.stage1)[[3]]

###draw
T.stage2 <- table(T.stage)%>%as.data.frame()
## Stacked graph

T.stage3 <- ggplot(T.stage2,aes(x=risk,y=Freq,fill=T.stage))+
  geom_bar(position = "fill",stat="identity",alpha=0.7)+
  theme_bw()+
  scale_fill_manual(values=brewer.pal(4, "Set2"))+ 
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, face = "bold",family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  face = "bold",family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18, color='black', face = "bold",family='Times'),
        title=element_text(size=30, color='black', face = "bold",family='Times'),
        strip.text = element_text(size = 14,family = "Times", face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x="",y="Percentage",fill="T.stage",title = "T.stage****")
T.stage3
# ggsave(T.stage3,filename = '04.T.stage.pdf',w=7,h=6)
# ggsave(T.stage3,filename = '04.T.stage.png',w=7,h=6)

###N.stage-------------------
N.stage <- train_phenotype2[,c("N.stage","risk")]
table(N.stage$N.stage)
N.stage$N.stage <- gsub('a','',N.stage$N.stage)
N.stage$N.stage <- gsub('b','',N.stage$N.stage)
N.stage$N.stage <- gsub('c','',N.stage$N.stage)
N.stage$N.stage <- gsub('NX',NA,N.stage$N.stage)

N.stage1 <- xtabs(~N.stage+risk,data=N.stage)
p <- chisq.test(N.stage1)[[3]]

###draw
N.stage2 <- table(N.stage)%>%as.data.frame()
## Stacked graph

N.stage3 <- ggplot(N.stage2,aes(x=risk,y=Freq,fill=N.stage))+
  geom_bar(position = "fill",stat="identity",alpha=0.7)+
  theme_bw()+
  scale_fill_manual(values=brewer.pal(4, "Set2"))+ 
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, face = "bold",family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  face = "bold",family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18, color='black', face = "bold",family='Times'),
        title=element_text(size=30, color='black', face = "bold",family='Times'),
        strip.text = element_text(size = 14,family = "Times", face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x="",y="Percentage",fill="N.stage",title = "N.stage***")
N.stage3
# ggsave(N.stage3,filename = '05.N.stage.pdf',w=7,h=6)
# ggsave(N.stage3,filename = '05.N.stage.png',w=7,h=6)

###M.stage-------------------
M.stage <- train_phenotype2[,c("M.stage","risk")]
table(M.stage$M.stage)
M.stage$M.stage <- gsub('a','',M.stage$M.stage)
M.stage$M.stage <- gsub('b','',M.stage$M.stage)
M.stage$M.stage <- gsub('MX',NA,M.stage$M.stage)

M.stage1 <- xtabs(~M.stage+risk,data=M.stage)
p <- chisq.test(M.stage1)[[3]]

###draw
M.stage2 <- table(M.stage)%>%as.data.frame()
## Stacked graph

M.stage3 <- ggplot(M.stage2,aes(x=risk,y=Freq,fill=M.stage))+
  geom_bar(position = "fill",stat="identity",alpha=0.7)+
  theme_bw()+
  scale_fill_manual(values=brewer.pal(4, "Set2"))+ 
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, face = "bold",family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  face = "bold",family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18, color='black', face = "bold",family='Times'),
        title=element_text(size=30, color='black', face = "bold",family='Times'),
        strip.text = element_text(size = 14,family = "Times", face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x="",y="Percentage",fill="M.stage",title = "M.stage****")
M.stage3
# ggsave(M.stage3,filename = '06.M.stage.pdf',w=7,h=6)
# ggsave(M.stage3,filename = '06.M.stage.png',w=7,h=6)



##Merge output -----------------------------------
library(patchwork)
all_clinical_index <- age3 + gender3 + race3 + T.stage3 + N.stage3 + M.stage3 +
  plot_layout(ncol = 3) &
  theme_bw() &
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold",size = 12))
all_clinical_index
ggsave(filename = '01.clinical.pdf',all_clinical_index,w=12,h=8)
ggsave(filename = '01.clinical.png',all_clinical_index,w=12,h=8)





##clinical----------
### Draw heat maps (combined with clinical information)

###Express data and grouping information
data <- read.csv("../00_Rowdata/tumor.tcga.csv",header=TRUE, row.names=1, check.names=FALSE)
gene <- readRDS("../08_LASSO/cox_gene.rds")
dat_gene <- data[gene,train_phenotype$sample]
###Clinical information collation

train_phenotype2 <- train_phenotype
train_phenotype2$risk <- ifelse(train_phenotype2$riskScore > median(train_phenotype2$riskScore),"High Risk","Low Risk")
colnames(train_phenotype2)

table(train_phenotype2$age)
train_phenotype2$age<- ifelse(train_phenotype2$age > 60,">60","<=60")

table(train_phenotype2$gender)

table(train_phenotype2$race)
train_phenotype2$race <- gsub('american indian or alaska native',"Unknown",train_phenotype2$race)
train_phenotype2$race <- gsub('not reported',"Unknown",train_phenotype2$race)

table(train_phenotype2$T.stage)
train_phenotype2$T.stage <- gsub('a','',train_phenotype2$T.stage)
train_phenotype2$T.stage <- gsub('b','',train_phenotype2$T.stage)
train_phenotype2$T.stage <- gsub('Tis',"Unknown",train_phenotype2$T.stage)

table(train_phenotype2$N.stage)
train_phenotype2$N.stage <- gsub('a','',train_phenotype2$N.stage)
train_phenotype2$N.stage <- gsub('b','',train_phenotype2$N.stage)
train_phenotype2$N.stage <- gsub('c','',train_phenotype2$N.stage)
train_phenotype2$N.stage <- gsub('NX',"Unknown",train_phenotype2$N.stage)

table(train_phenotype2$M.stage)
train_phenotype2$M.stage <- gsub('a','',train_phenotype2$M.stage)
train_phenotype2$M.stage <- gsub('b','',train_phenotype2$M.stage)
train_phenotype2$M.stage <- gsub('MX',"Unknown",train_phenotype2$M.stage)


table(train_phenotype2$OS)
train_phenotype2$OS <- ifelse(train_phenotype2$OS ==1,"Dead","Alive")

train_phenotype2 <- train_phenotype2[order(train_phenotype2$risk,train_phenotype2$OS,train_phenotype2$age,train_phenotype2$gender),]
heat.group <- train_phenotype2



# heat.dat <- t(dat_gene) %>% data.frame()
heat.dat <- dat_gene

colnames(heat.group)
heat.group<-dplyr::select(heat.group,c('risk','OS','age','race','gender','T.stage','N.stage','M.stage','sample'))
rownames(heat.group) <- heat.group$sample
heat.group <- heat.group[,-9]
heat.group$risk <- as.factor(heat.group$risk)
heat.group$age<-as.factor(heat.group$age)
heat.group$gender <- as.factor(heat.group$gender)
heat.group$race <- as.factor(heat.group$race)
heat.group$OS <- as.factor(heat.group$OS)
heat.group$T.stage <- as.factor(heat.group$T.stage)
heat.group$N.stage <- as.factor(heat.group$N.stage)
heat.group$M.stage <- as.factor(heat.group$M.stage)

rt_dat<-log2(heat.dat+1)
rt_dat <- rt_dat[gene,rownames(heat.group)]
rt_dat <- as.matrix(rt_dat)
# rt_dat <- data.frame(t(scale(t(rt_dat))))
# colnames(rt_dat)<-gsub('.','-',colnames(rt_dat),fixed = T)
ann_colors<-list(
  risk=c('Low Risk'='#00FF7F','High Risk'='#FF5336'),
  gender=c('male'='#FF8C00','female'='#20B2AA'),
  T.stage=c('T1'='#33CCCC','T2'='#CCCCFF','T3'='#FF9966','T4'='#FFCCFF','Unknown'='#40E0D0'),
  race=c('black or african american'='#FF8C00','white'='blue','Unknown'='green','asian' = '#FF1493'),
  age=c('<=60'='#00CED1','>60'='#EE82EE'),
  N.stage=c('N0'='#FFE4E1','N1'='#FF8C00','N2'='#FFB6C1','Unknown'='#40E0D0'),
  M.stage=c('M0'='#D2691E','M1'='#20B2AA','Unknown'='#40E0D0'),
  OS=c('Dead'='red','Alive'='green')
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
png("heatmap.png",w=8,h=8,units = "in",res = 300)
print(hp)
dev.off()
pdf("heatmap.pdf",w=8,h=8)
print(hp)
dev.off()


