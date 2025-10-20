rm(list = ls())

setwd("/data/nas1/yuanyt/project/05.YQQL-0602-6/")
if (! dir.exists("./00_Rowdata")){
  dir.create("./00_Rowdata")
}
setwd("./00_Rowdata")

library(TCGAbiolinks)
library(readr)
library(readxl)
library(tidyverse)


protein <- read.table("/data/nas1/xieyunpeng/project/pipeline/PCG/protein_coding.txt",header = T,sep = "\t")

### 01 Data processing
###Read the data downloaded from xena count
tcga.COAD.expr <- read_tsv("../00_Rowdata/TCGA-COAD.htseq_counts.tsv.gz")
tcga.READ.expr <- read_tsv("../00_Rowdata/TCGA-READ.htseq_counts.tsv.gz")
tcga.expr <- merge(tcga.COAD.expr,tcga.READ.expr,by="Ensembl_ID")
tcga.expr<-as.data.frame(tcga.expr)
rownames(tcga.expr)<-tcga.expr[,1]
tcga.expr<-tcga.expr[,-1]
tcga.expr <- 2^tcga.expr-1


## Perform id conversion on the data
genecode<-read.table(file = '/data/nas1/luchunlin/pipeline/GENEANNO/gencode.v22.annotation.gene.probeMap')
probe2symbol<-genecode[,(1:2)]
colnames(probe2symbol)<-c('ID','symbol')
probe2symbol<-probe2symbol[-1,]
dat.tcga<-tcga.expr
dat.tcga$ID <- rownames(dat.tcga)
dat.tcga$ID<-as.character(dat.tcga$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)

####Table data processing
dat.tcga<-dat.tcga %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     ## Remove redundant information
  dplyr::select(symbol,everything())%>%     ## Rearrange
  mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    ## Calculate the average
  arrange(desc(rowMean))%>%       ## Sort the average values of the expressions from largest to smallest
  distinct(symbol,.keep_all = T)%>%      ##symbol left the first one
  dplyr::select(-rowMean)%>%     ## Reverse selection removes the column "rowMean"
  tibble::column_to_rownames(colnames(.)[1])   ## Change the first column to the row name and delete it
dim(dat.tcga)


### get mRNA
dat.tcga <- dat.tcga[rownames(dat.tcga)%in%protein$gene_name,]


### Screen cancerous tissues and remove adjacent tissues. 01 to 09 are tumors, and 10 to 19 are normal controls
group <- data.frame(colnames(dat.tcga))  # Take the sample id of the first row
for (i in 1:length(group[,1])) {
  num=as.numeric(as.character(substring(group[i,1],14,15)))
  if(num %in% seq(1,9)){group[i,2]="CRC"}
  if(num %in% seq(10,29)){group[i,2]="Control"}
}


###tumor group
names(group) <- c("id","group")
group$group <- as.factor(group$group)
exp_group <- group$id[which(group$group == "CRC")]
exp_tumor<-dat.tcga[,which(colnames(dat.tcga) %in% exp_group)]
exp_tumor<-as.data.frame(exp_tumor)  

###control group
exp_control<-dat.tcga[,which(!colnames(dat.tcga)%in%exp_group)]
exp_control<-as.data.frame(exp_control)
control.sample<-colnames(exp_control)
tumor.sample<-colnames(exp_tumor)

dat.final<-cbind(exp_control,exp_tumor)
write.csv(dat.final,file = 'dat.tcga.csv')
write.csv(exp_tumor,file = 'tumor.tcga.csv')
write.csv(group,file = 'tcga.group.csv')


###fpkm
fpkm.COAD.expr <- read_tsv("../00_Rowdata/TCGA-COAD.htseq_fpkm.tsv.gz")
fpkm.READ.expr <- read_tsv("../00_Rowdata/TCGA-READ.htseq_fpkm.tsv.gz")
fpkm.expr <- merge(fpkm.COAD.expr,fpkm.READ.expr,by="Ensembl_ID")
fpkm.expr<-as.data.frame(fpkm.expr)
rownames(fpkm.expr)<-fpkm.expr[,1]
fpkm.expr<-fpkm.expr[,-1]
fpkm.expr <- 2^fpkm.expr-1

## Perform id conversion on the data
dat_fpkm<-fpkm.expr
dat_fpkm$ID <- rownames(dat_fpkm)
dat_fpkm$ID<-as.character(dat_fpkm$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)

dat_fpkm<-dat_fpkm %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     ##  Remove redundant information
  dplyr::select(symbol,everything())%>%     ##  Rearrange
  mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    ##  Calculate the average
  arrange(desc(rowMean))%>%       ##  Sort the average values of the expressions from largest to smallest
  distinct(symbol,.keep_all = T)%>%      ## symbol left the first one
  dplyr::select(-rowMean)%>%     ##  Reverse selection removes the column "rowMean"
  tibble::column_to_rownames(colnames(.)[1])   ##  Change the first column to the row name and delete it
dim(dat_fpkm)

####mRNA
dat_fpkm <- dat_fpkm[rownames(dat_fpkm)%in%protein$gene_name,]



###tumor group
exp_tumor<-dat_fpkm[,which(colnames(dat_fpkm)%in%exp_group)]
exp_tumor<-as.data.frame(exp_tumor)  

###control group
exp_control<-dat_fpkm[,which(!colnames(dat_fpkm)%in%exp_group)]
exp_control<-as.data.frame(exp_control)
control.sample<-colnames(exp_control)
tumor.sample<-colnames(exp_tumor)


dat.final<-cbind(exp_control,exp_tumor)
write.csv(dat.final,file = 'fpkm.tcga.csv')
write.csv(exp_tumor,file = 'fpkm.tumor.csv')



###Read the data downloaded from GEO
library(GEOquery)
library(tidyverse)
library(lance)

###GSE103479--------------------------------
gset<-getGEO("GSE103479",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)

expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL23985",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','Gene.Symbol')%>%
  filter('Gene.Symbol'!='')%>%
  separate('Gene.Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
library(tidyverse)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>% 
  dplyr::select(-ID)%>%     ##  Remove redundant information
  dplyr::select(symbol,everything())%>%     ##  Rearrange
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## Calculate the average
  arrange(desc(rowMean))%>%       ##  Sort the average values of the expressions from largest to smallest
  distinct(symbol,.keep_all = T)%>%      ## symbol left the first one
  dplyr::select(-rowMean)%>%     ## Reverse selection removes the column "rowMean"
  tibble::column_to_rownames(colnames(.)[1])   ## Change the first column to the row name and delete it



###Group information
pd<-pData(a)
survival<-data.frame(sample=pd$geo_accession,
                  OS = pd$`status alive.dead:ch1`,
                  OS.time = pd$`overall survival time:ch1`)
survival <- survival[which(survival$OS.time != "NA"),]

table(survival$OS)
library(stringr)
survival$OS<-ifelse(survival$OS=='Alive','0','1')
survival$OS.time<- as.numeric(survival$OS.time) * 30

dat<-dat[,survival$sample]
write.csv(dat,file = 'dat(GSE103479).csv')
write.csv(survival,file = 'survival(GSE103479).csv')




###GSE17536----------------------
gset<-getGEO("GSE17536",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)

expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL570",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
library(tidyverse)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>% 
  dplyr::select(-ID)%>%     ##  Remove redundant information
  dplyr::select(symbol,everything())%>%     ##  Rearrange
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## Calculate the average
  arrange(desc(rowMean))%>%       ##  Sort the average values of the expressions from largest to smallest
  distinct(symbol,.keep_all = T)%>%      ## symbol left the first one
  dplyr::select(-rowMean)%>%     ## Reverse selection removes the column "rowMean"
  tibble::column_to_rownames(colnames(.)[1])   ## Change the first column to the row name and delete it


###Group information
pd<-pData(a)
survival<-data.frame(sample=pd$geo_accession,
                     OS = pd$`overall_event (death from any cause):ch1`,
                     OS.time = pd$`overall survival follow-up time:ch1`)
survival <- survival[which(survival$OS.time != "NA"),]

table(survival$OS)
library(stringr)
survival$OS<-ifelse(survival$OS=='no death','0','1')
survival$OS.time<- as.numeric(survival$OS.time) * 30

dat<-dat[,survival$sample]
write.csv(dat,file = 'dat(GSE17536).csv')
write.csv(survival,file = 'survival(GSE17536).csv')


###GSE17537----------------------
gset<-getGEO("GSE17537",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)

expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL570",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
library(tidyverse)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>% 
  dplyr::select(-ID)%>%     ##  Remove redundant information
  dplyr::select(symbol,everything())%>%     ##  Rearrange
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## Calculate the average
  arrange(desc(rowMean))%>%       ##  Sort the average values of the expressions from largest to smallest
  distinct(symbol,.keep_all = T)%>%      ## symbol left the first one
  dplyr::select(-rowMean)%>%     ## Reverse selection removes the column "rowMean"
  tibble::column_to_rownames(colnames(.)[1])   ## Change the first column to the row name and delete it


###Group information
pd<-pData(a)
survival<-data.frame(sample=pd$geo_accession,
                     OS = pd$`overall_event (death from any cause):ch1`,
                     OS.time = pd$`overall survival follow-up time:ch1`)
survival <- survival[which(survival$OS.time != "NA"),]

table(survival$OS)
library(stringr)
survival$OS<-ifelse(survival$OS=='no death','0','1')
survival$OS.time<- as.numeric(survival$OS.time) * 30

dat<-dat[,survival$sample]
write.csv(dat,file = 'dat(GSE17537).csv')
write.csv(survival,file = 'survival(GSE17537).csv')


###GSE17538----------------------
gset<-getGEO("GSE17538",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)

expr<-as.data.frame(exprs(gset[[2]]))
gpl<-getGEO("GPL570",destdir = '.')
a=gset[[2]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
library(tidyverse)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>% 
  dplyr::select(-ID)%>%     ##  Remove redundant information
  dplyr::select(symbol,everything())%>%     ##  Rearrange
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## Calculate the average
  arrange(desc(rowMean))%>%       ##  Sort the average values of the expressions from largest to smallest
  distinct(symbol,.keep_all = T)%>%      ## symbol left the first one
  dplyr::select(-rowMean)%>%     ## Reverse selection removes the column "rowMean"
  tibble::column_to_rownames(colnames(.)[1])   ## Change the first column to the row name and delete it


###Group information
pd<-pData(a)
survival<-data.frame(sample=pd$geo_accession,
                     OS = pd$`overall_event (death from any cause):ch1`,
                     OS.time = pd$`overall survival follow-up time:ch1`)
survival <- survival[which(survival$OS.time != "NA"),]

table(survival$OS)
library(stringr)
survival$OS<-ifelse(survival$OS=='no death','0','1')
survival$OS.time<- as.numeric(survival$OS.time) * 30

dat<-dat[,survival$sample]
write.csv(dat,file = 'dat(GSE17538).csv')
write.csv(survival,file = 'survival(GSE17538).csv')


###GSE28722----------------------
gset<-getGEO("GSE28722",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)

expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL13425",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','ORF')%>%
  filter('ORF'!='')%>%
  separate('ORF',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
library(tidyverse)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>% 
  dplyr::select(-ID)%>%     ##  Remove redundant information
  dplyr::select(symbol,everything())%>%     ##  Rearrange
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## Calculate the average
  arrange(desc(rowMean))%>%       ##  Sort the average values of the expressions from largest to smallest
  distinct(symbol,.keep_all = T)%>%      ## symbol left the first one
  dplyr::select(-rowMean)%>%     ## Reverse selection removes the column "rowMean"
  tibble::column_to_rownames(colnames(.)[1])   ## Change the first column to the row name and delete it


###Group information
pd<-pData(a)
survival<-data.frame(sample=pd$geo_accession,
                     OS = pd$`overall survival censor (1-censored,0-non-censored):ch2`,
                     OS.time = pd$`overall survival (years):ch2`)
survival <- survival[which(survival$OS.time != "NA"),]

table(survival$OS)
library(stringr)
survival$OS.time<- as.numeric(survival$OS.time) * 365

dat<-dat[,survival$sample]
write.csv(dat,file = 'dat(GSE28722).csv')
write.csv(survival,file = 'survival(GSE28722).csv')


###GSE29621----------------------
gset<-getGEO("GSE29621",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)

expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL570",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
library(tidyverse)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>% 
  dplyr::select(-ID)%>%     ##  Remove redundant information
  dplyr::select(symbol,everything())%>%     ##  Rearrange
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## Calculate the average
  arrange(desc(rowMean))%>%       ##  Sort the average values of the expressions from largest to smallest
  distinct(symbol,.keep_all = T)%>%      ## symbol left the first one
  dplyr::select(-rowMean)%>%     ## Reverse selection removes the column "rowMean"
  tibble::column_to_rownames(colnames(.)[1])   ## Change the first column to the row name and delete it


###Group information
pd<-pData(a)
survival<-data.frame(sample=pd$geo_accession,
                     OS = pd$`os event:ch1`,
                     OS.time = pd$`overall survival (os):ch1`)
survival <- survival[which(survival$OS.time != "NA"),]

table(survival$OS)
library(stringr)
survival$OS<-ifelse(survival$OS=='alive','0','1')
survival$OS.time<- as.numeric(survival$OS.time) * 30

dat<-dat[,survival$sample]
write.csv(dat,file = 'dat(GSE29621).csv')
write.csv(survival,file = 'survival(GSE29621).csv')


###GSE72970----------------------
gset<-getGEO("GSE72970",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)

expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL570",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
library(tidyverse)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>% 
  dplyr::select(-ID)%>%     ##  Remove redundant information
  dplyr::select(symbol,everything())%>%     ##  Rearrange
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## Calculate the average
  arrange(desc(rowMean))%>%       ##  Sort the average values of the expressions from largest to smallest
  distinct(symbol,.keep_all = T)%>%      ## symbol left the first one
  dplyr::select(-rowMean)%>%     ## Reverse selection removes the column "rowMean"
  tibble::column_to_rownames(colnames(.)[1])   ## Change the first column to the row name and delete it


###Group information
pd<-pData(a)
survival<-data.frame(sample=pd$geo_accession,
                     OS = pd$`os censored:ch1`,
                     OS.time = pd$`os:ch1`)
survival <- survival[which(survival$OS.time != "NA"),]

table(survival$OS)
library(stringr)
survival$OS.time<- as.numeric(survival$OS.time) * 30

dat<-dat[,survival$sample]
write.csv(dat,file = 'dat(GSE72970).csv')
write.csv(survival,file = 'survival(GSE72970).csv')


###GSE39582-----------------
gset<-getGEO("GSE39582",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)

expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL570",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
library(tidyverse)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>% 
  dplyr::select(-ID)%>%     ##  Remove redundant information
  dplyr::select(symbol,everything())%>%     ##  Rearrange
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## Calculate the average
  arrange(desc(rowMean))%>%       ##  Sort the average values of the expressions from largest to smallest
  distinct(symbol,.keep_all = T)%>%      ## symbol left the first one
  dplyr::select(-rowMean)%>%     ## Reverse selection removes the column "rowMean"
  tibble::column_to_rownames(colnames(.)[1])   ## Change the first column to the row name and delete it


###Group information
pd<-pData(a)
survival<-data.frame(sample=pd$geo_accession,
                     OS = pd$`os.event:ch1`,
                     OS.time = pd$`os.delay (months):ch1`)
survival <- survival[which(survival$OS != "N/A"),]

table(survival$OS)
library(stringr)
survival$OS.time<- as.numeric(survival$OS.time) * 30

dat<-dat[,survival$sample]
write.csv(dat,file = 'dat(GSE39582).csv')
write.csv(survival,file = 'survival(GSE39582).csv')






