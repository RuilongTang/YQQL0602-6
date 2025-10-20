rm(list = ls())

setwd("/data/nas1/yuanyt/project/05.YQQL-0602-6/")
if (! dir.exists("./16_Wilcoxon")){
  dir.create("./16_Wilcoxon")
}
setwd("./16_Wilcoxon")


###Read data -------------------------
survival <- read.csv("../02_ssGSEA/icd.group.csv",header = T,row.names = 1)
survival <- survival[,c(1,2,4)]

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
colnames(train_phenotype2)




## 01 Age-----
age <- train_phenotype2[,c("age","immunogenic.cell.death")]
my_comparisons <- list(c(">=60","<60"))
age$age <- ifelse(age$age<60,">=60","<60")
age <- na.omit(age)
age<-ggboxplot(age, x = "age", y = "immunogenic.cell.death",
               color = "age", palette = c("#F2F921", "#1E90FF"),
               add = "jitter",
               short.panel.labs = T,
               ggtheme = theme_bw())+
  scale_y_continuous(name = "",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("Age") +
  geom_signif(comparisons = my_comparisons,
              test = t.test,
              map_signif_level = T,
              y_position = 1.2)+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        # strip.background = element_blank(),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 20, face = "bold"),
        text = element_text(size = 30, face = "bold"))+
  guides(fill='none')
age
# ggsave(filename = '01.age.pdf',age,w=7,h=6)
# ggsave(filename = '01.age.png',age,w=7,h=6)

## 02 gender-----
gender <- train_phenotype2[,c("gender","immunogenic.cell.death")]
my_comparisons <- list(c("female","male"))
table(gender$gender)
gender <- na.omit(gender)
gender<-ggboxplot(gender, x = "gender", y = "immunogenic.cell.death",
                  color = "gender", palette = c("#F2F921", "#1E90FF"),
                  add = "jitter",
                  short.panel.labs = T,
                  ggtheme = theme_bw())+
  scale_y_continuous(name = "",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("gender") +
  geom_signif(comparisons = my_comparisons,
              test = t.test,
              map_signif_level = T,
              y_position = 1.2)+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        # strip.background = element_blank(),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 20, face = "bold"),
        text = element_text(size = 30, face = "bold"))+
  guides(fill='none')
gender
# ggsave(filename = '02.gender.pdf',gender,w=7,h=6)
# ggsave(filename = '02.gender.png',gender,w=7,h=6)


## 03 race-----
race <- train_phenotype2[,c("race","immunogenic.cell.death")]
race <- race[which(race$race %in%c("white","black or african american")),]
table(race$race)
race <- na.omit(race)
my_comparisons <- list(c("white","black or african american"))
race<-ggboxplot(race, x = "race", y = "immunogenic.cell.death",
                color = "race", palette = c("#F2F921", "#1E90FF"),
                add = "jitter",
                short.panel.labs = T,
                ggtheme = theme_bw())+
  scale_y_continuous(name = "",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("race") +
  geom_signif(comparisons = my_comparisons,
              map_signif_level = TRUE,
              test = t.test,
              y_position = 1.2)+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        # strip.background = element_blank(),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 20, face = "bold"),
        text = element_text(size = 30, face = "bold"))+
  guides(fill='none')
race
# ggsave(filename = '03.race.pdf',race,w=7,h=6)
# ggsave(filename = '03.race.png',race,w=7,h=6)


## 04 T.stage-----
T.stage <- train_phenotype2[,c("T.stage","immunogenic.cell.death")]
table(T.stage$T.stage)
T.stage$T.stage <- gsub('a','',T.stage$T.stage)
T.stage$T.stage <- gsub('b','',T.stage$T.stage)
T.stage$T.stage <- gsub('Tis',NA,T.stage$T.stage)
T.stage$T.stage <- gsub('T1','T1/2',T.stage$T.stage)
T.stage$T.stage <- gsub('T2','T1/2',T.stage$T.stage)
T.stage$T.stage <- gsub('T3','T3/4',T.stage$T.stage)
T.stage$T.stage <- gsub('T4','T3/4',T.stage$T.stage)

my_comparisons <- list(c("T1/2","T3/4"))
T.stage <- na.omit(T.stage)
T.stage<-ggboxplot(T.stage, x = "T.stage", y = "immunogenic.cell.death",
                    color = "T.stage", palette = c("#F2F921", "#1E90FF"),
                    add = "jitter",
                    short.panel.labs = T,
                    ggtheme = theme_bw())+
  scale_y_continuous(name = "",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("T.stage") +
  geom_signif(comparisons = my_comparisons,
              map_signif_level = TRUE,
              test = t.test,
              y_position = 1.2)+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        # strip.background = element_blank(),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 20, face = "bold"),
        text = element_text(size = 30, face = "bold"))+
  guides(fill='none')
T.stage
# ggsave(filename = '04.T.stage.pdf',T.stage,w=7,h=6)
# ggsave(filename = '04.T.stage.png',T.stage,w=7,h=6)


## 05 N.stage-----
N.stage <- train_phenotype2[,c("N.stage","immunogenic.cell.death")]
table(N.stage$N.stage)
N.stage$N.stage <- gsub('a','',N.stage$N.stage)
N.stage$N.stage <- gsub('b','',N.stage$N.stage)
N.stage$N.stage <- gsub('c','',N.stage$N.stage)
N.stage$N.stage <- gsub('NX',NA,N.stage$N.stage)
N.stage$N.stage <- gsub('N1','N1/2',N.stage$N.stage)
N.stage$N.stage <- gsub('N2','N1/2',N.stage$N.stage)


N.stage <- na.omit(N.stage)
my_comparisons <- list(c("N0","N1/2"))
N.stage<-ggboxplot(N.stage, x = "N.stage", y = "immunogenic.cell.death",
                     color = "N.stage", palette = c("#F2F921", "#1E90FF"),
                     add = "jitter",
                     short.panel.labs = T,
                     ggtheme = theme_bw())+
  scale_y_continuous(name = "",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("N.stage") +
  geom_signif(comparisons = my_comparisons,
              map_signif_level = TRUE,
              test = t.test,
              y_position = 1.2)+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        # strip.background = element_blank(),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 20, face = "bold"),
        text = element_text(size = 30, face = "bold"))+
  guides(fill='none')
N.stage
# ggsave(filename = '05.N.stage.pdf',N.stage,w=7,h=6)
# ggsave(filename = '05.N.stage.png',N.stage,w=7,h=6)



## 06 M.stage-----
M.stage <- train_phenotype2[,c("M.stage","immunogenic.cell.death")]
table(M.stage$M.stage)
M.stage$M.stage <- gsub('a','',M.stage$M.stage)
M.stage$M.stage <- gsub('b','',M.stage$M.stage)
M.stage$M.stage <- gsub('MX',NA,M.stage$M.stage)
M.stage <- na.omit(M.stage)

my_comparisons <- list(c("M0","M1"))
M.stage<-ggboxplot(M.stage, x = "M.stage", y = "immunogenic.cell.death",
               color = "M.stage", palette = c("#F2F921", "#1E90FF"),
               add = "jitter",
               short.panel.labs = T,
               ggtheme = theme_bw())+
  scale_y_continuous(name = "",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("M.stage") +
  geom_signif(comparisons = my_comparisons,
              map_signif_level = TRUE,
              test = t.test,
              y_position = 1.2)+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        # strip.background = element_blank(),
        axis.title.y = element_text(size = 24, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 16, face = "bold"),
        axis.text.x = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 20, face = "bold"),
        text = element_text(size = 30, face = "bold"))+
  guides(fill='none')
M.stage
# ggsave(filename = '06.M.stage.pdf',M.stage,w=9,h=8)
# ggsave(filename = '06.M.stage.png',M.stage,w=9,h=8)





##Merge output-----------------
library(patchwork)
all_clinical_index <- age + gender + race + T.stage + N.stage + M.stage+
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


