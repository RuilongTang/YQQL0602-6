rm(list = ls())

setwd("/data/nas1/yuanyt/project/05.YQQL-0602-6/")
if (! dir.exists("./03_WGCNA")){
  dir.create("./03_WGCNA")
}
setwd("./03_WGCNA")


data <- read.csv("../00_Rowdata/fpkm.tumor.csv",header=TRUE, row.names=1, check.names=FALSE)%>% lc.tableToNum
data <- log2(data+1)

group<-read.csv('../02_ssGSEA/icd.group.csv',row.names = 1)
group <- group[,c("sample","group")]
table(group$group)


###Loading package
library(WGCNA)
library(reshape2)
library(stringr)
options(stringsAsFactors = F)
enableWGCNAThreads()

exprMat<-data[,group$sample]
dim(exprMat)
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
# When calculating the correlation of binary variables, such as sample trait information, or when gene expression is highly dependent on the disease state, the following parameters need to be set
maxPOutliers = ifelse(corType=="pearson",1,0.05)
# When associating the binary variables of sample properties, set
robustY = ifelse(corType=="pearson",T,F)
dataExpr <- exprMat[rownames(data),]



## 01 Data filtering-----
##Screen the genes with the top 75% median absolute deviation (MAD), with at least a MAD greater than 0.01
## After screening, the computational load will be reduced and some information will also be lost
## Alternatively, no filtering is required; just make the MAD greater than 0
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad >
                                max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
# dataExprVar <- dataExpr
## Convert it into a matrix with samples in rows and genes in columns
dataExpr <- as.data.frame(t(dataExprVar))
## Detect missing values
gsg = goodSamplesGenes(dataExpr, verbose = 3)
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
dim(dataExpr)
# [1] 635 13603

## 02 Soft threshold screening ----
## Sample clustering to check if there are any outlier samples

tree=hclust(dist(dataExpr),method ='complete')   ##Observe whether there are any outlier samples  c("single", "complete", "average", "mcquitty", "ward.D", "centroid", "median", "ward.D2")
#The hclust() function is a function in the stats package that can implement hierarchical clustering based on the distance matrix.
par(pin = c(4,4), mar = c(6,6,6,1),family='serif')
plot(tree,xlab="", sub="", main="Sample Clustering",
     labels=F,
     cex=1.0,  ##label size
     font=2,
     cex.axis=1.6,  ##The scaling factor of the scale text on the coordinate axis. Similar to cex.
     cex.lab=1.8,   ##The scaling factor of the scale text on the coordinate axis. Similar to cex.
     cex.main=1.8,   ##The zoom ratio of the title
     font.axis = 2,
     font.lab = 2,
     font.main = 2,
     font.sub =2)
abline(h =500, col = 'red')
dev.off()

# ## As can be seen from the above picture, it seems abnormal and should be considered for removal.
# clust = cutreeStatic(tree, cutHeight =210, minSize = 10)
# table(clust)  ###There are no outlier samples
# ## The expression matrix after removing abnormal samples
# datExpr = datExpr1[clust == 1, ]
# nGenes = ncol(datExpr)      #Number of genes
# nSample =nrow(datExpr) #Sample name
# nGenes ;nSample
SampleName<-rownames(dataExpr)

## Phenotypic data
condition <- group
condition$High<-ifelse(condition$group=='High',1,0)
condition$Low<-ifelse(condition$group=='Low',1,0)
rownames(condition) <- condition$sample
condition<-condition[,-c(1:2)]
# library(tidyverse)
datTraits<-condition

dataExpr[1:4,1:4]

traitColors = numbers2colors(datTraits, signed = FALSE)
tree2=hclust(dist(dataExpr),method ='complete')   ##Visualize the connection between phenotypic data and gene expression level data, and reconstruct the sample clustering tree  c("single", "complete", "average", "mcquitty", "ward.D", "centroid", "median", "ward.D2")

pdf(file='01.sampleClustering.pdf',w=10,h=7)
par(pin = c(4,4), mar = c(6,6,6,1),family='serif')
plotDendroAndColors(tree2,
                    traitColors,
                    dendroLabels = FALSE,
                    hang = 0.03,
                    cex.dendroLabels = 1,
                    # addGuide = TRUE,   ##Add grid lines
                    font=2,
                    guideHang = 0.05,
                    cex.axis=1.4,  ##The scaling factor of the scale text on the coordinate axis. Similar to cex.
                    cex.lab=1.6,   ##The scaling factor of the scale text on the coordinate axis. Similar to cex.
                    cex.main=1.6,   ##The zoom ratio of the title.
                    font.axis = 2,
                    font.lab = 2,
                    font.main = 2,
                    font.sub =2,
                    cex.colorLabels=1,
                    groupLabels = colnames(datTraits),
                    main = "Sample Clustering and trait heatmap")
dev.off()

png(file='01.sampleClustering.png',w=10,h=7,units='in',res=600,bg='white')
par(pin = c(4,4), mar = c(6,6,6,1),family='serif')
plotDendroAndColors(tree2,
                    traitColors,
                    dendroLabels = FALSE,
                    hang = 0.03,
                    cex.dendroLabels = 1,
                    # addGuide = TRUE,   ##Add grid lines
                    font=2,
                    guideHang = 0.05,
                    cex.axis=1.4,  ##The scaling factor of the scale text on the coordinate axis. Similar to cex.
                    cex.lab=1.6,   ##The scaling factor of the scale text on the coordinate axis. Similar to cex.
                    cex.main=1.6,   ##The zoom ratio of the title.
                    font.axis = 2,
                    font.lab = 2,
                    font.main = 2,
                    font.sub =2,
                    cex.colorLabels=1,
                    groupLabels = colnames(datTraits),
                    main = "Sample Clustering and trait heatmap")
dev.off()

# Determine cluster under the line
# Pruning algorithm: cutHeight - the height of pruned branches; minSize - the minimum number of the cluster
# clust = cutreeStatic(sampleTree, cutHeight = 70, minSize = 10)
# table(clust)
# # clust 1 contains the samples we want to keep.
# keepSamples = (clust==1)
## Data that meets the requirements
# dataExpr = dataExpr[keepSamples,]
## Extract column
#nGenes = ncol(dataExpr)
## Extract line
#nSamples = nrow(dataExpr)


# Set the selection range of network construction parameters and calculate the scale-free distribution topology matrixpowers <-  c(c(1:20))
sft <- pickSoftThreshold(dataExpr, powerVector=powers,
                         networkType=type, verbose=5)

# The horizontal axis represents Soft threshold (power), and the vertical axis is the evaluation parameter of the scale-free network. The higher the value, the more the network conforms to the scale-free feature (non-scale).
pdf('02.softThreshold.pdf', w=12, h=8)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fitsigned R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# Screening criteria. R - square = 0.85
abline(h=0.85,col="red")
# Soft threshold and average connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
     cex=cex1, col="red")
dev.off()

png('02.softThreshold.png', w=800, h=500)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fitsigned R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# Screening criteria. R - square = 0.85
abline(h=0.85,col="red")
# Soft threshold and average connectivity
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
     cex=cex1, col="red")
dev.off()


power = sft$powerEstimate
power
## 15


## 03 One-step network construction---------
## One-step network construction and module detection##
# power: The soft threshold calculated in the previous step
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 50,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType,
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = "DiffGene_TOM",
                       verbose = 3)
# According to the number of genes in the module, they are arranged in descending order and numbered successively as' 1- maximum number of modules'.
# **0 (grey)** indicates genes that have not been assigned to any module.
table(net$colors)  
#     0     1     2     3     4     5     6 
# 10178  1617   780   452   402    90    84 
## 04-4 Hierarchical clustering number display each module
## The gray ones are genes that have not been classified into modules.
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
table(moduleLabels)
table(moduleColors)
# Plot the dendrogram and the module colors underneath
# If you are not satisfied with the result, you can also recutBlockwiseTrees to save computing time
# The plotdendroandcolors function accepts a clustered object and the colors corresponding to all the individuals contained within that object.
png(filename = "03.cluster_dendrogram.png", height = 600, width = 800)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
pdf(file = "03.cluster_dendrogram.pdf", height = 6, width = 8)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


## 04-5 Draw a heat map of the correlations between modules-------
# module eigengene can draw a line graph as a display of the gene expression trend of each module
MEs = net$MEs
### There's no need to recalculate. Just change the following names
### The official tutorial is recalculated
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
# The correlation graph among each module obtained by clustering based on the expression levels between genes
# Set the margins at the bottom, left, top and right of marDendro/marHeatmap
sizeGrWindow(6,6)
pdf('04.Eigengene_dendrogram_and_heatmap.pdf',w=5,h=5)
plotEigengeneNetworks(MEs_col, 
                      setLabels = "Eigengene dendrogram and heatmap", 
                      marDendro = c(1,3,4,4),
                      marHeatmap = c(3,3,0,2), 
                      plotDendrograms = T, 
                      xLabelsAngle = 90)

dev.off()
png('04.Eigengene_dendrogram_and_heatmap.png',w=500,h=500)
plotEigengeneNetworks(MEs_col, 
                      setLabels = "Eigengene dendrogram and heatmap", 
                      marDendro = c(1,3,4,4),
                      marHeatmap = c(3,3,0,2), 
                      plotDendrograms = T, 
                      xLabelsAngle = 90)

dev.off()
# Plot the dendrogram
plotEigengeneNetworks(MEs_col, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


## 04-6 Associated phenotypic data------
group_traits<-group
rownames(group_traits)<-group_traits$sample
group_traits<-group_traits[rownames(dataExpr),]
group_traits<-group_traits[,-1]
group_traits<-as.data.frame(group_traits)
colnames(group_traits)<-"Group"
rownames(group_traits)<-rownames(dataExpr)
datTraits=data.frame(samples=rownames(dataExpr),subtype=group_traits)
design_traits<-model.matrix(~0+datTraits$Group)
design_traits<-as.data.frame(design_traits)
colnames(design_traits)=levels(factor(datTraits$Group))
moduleColors <- labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(dataExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0); ##The ME value matrix of modules of different colors (sample vs module)
moduleTraitCor = cor(MEs, design_traits , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(10,10)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

pdf("05.wgcna.Module-trait.heatmap.pdf", width = 8, height =12)
par(mar = c(6, 8.5, 6, 3))
labeledHeatmap(Matrix = moduleTraitCor[-16,],
               xLabels = names(design_traits),
               yLabels = names(MEs[,-16]),
               ySymbols = names(MEs[,-16]),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix[-16,],
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

png("05.wgcna.Module-trait.heatmap.png", width = 8, height =12,unit='in',res=600,bg='white')
par(mar = c(6, 8.5, 6, 3))
labeledHeatmap(Matrix = moduleTraitCor[-16,],
               xLabels = names(design_traits),
               yLabels = names(MEs[,-16]),
               ySymbols = names(MEs[,-16]),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix[-16,],
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


save.image(file = "WGCNA.RData")
load(file = "WGCNA.RData")
##MEbrown
module='brown'
probes=colnames(dataExpr)
inModule=(moduleColors==module)
modProbes=probes[inModule]
modGenes<-as.data.frame(modProbes)
colnames(modGenes)<-'modgene'
write.csv(modGenes,file = 'brown.modgene.csv')

##MEblue
module='blue'
probes=colnames(dataExpr)
inModule=(moduleColors==module)
modProbes=probes[inModule]
modGenes<-as.data.frame(modProbes)
colnames(modGenes)<-'modgene'
write.csv(modGenes,file = 'blue.modgene.csv')
#
### The correlation between modules and genes
if (corType=="pearson") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
  as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}
# Calculate the correlation matrix between traits and genes

## Only continuous traits can be calculated. If it is a discrete variable, it should be converted to a 0-1 matrix when constructing the sample table.

if (corType=="pearson") {
  geneTraitCor = as.data.frame(cor(dataExpr, design_traits, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, design_traits, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}
### brown---------
module = "brown"
pheno = "T"
modNames = substring(colnames(MEs_col), 3)
# Get the columns of attention
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(design_traits))
# Obtain the genes within the module
moduleGenes = moduleColors == module
MM<-abs(geneModuleMembership[moduleGenes,module_column])
GS<-abs(geneTraitCor[moduleGenes, 1])
c<-as.data.frame(cbind(MM,GS))
rownames(c)=modGenes$modgene
brown_hub <- subset(c,c$MM > 0.5 & c$GS > 0.2)
sizeGrWindow(7, 7)
par(mfrow = c(1,1),
    mai = c(1,1,1,0.5))
pdf('06.lightgreen_gene_cor.pdf',w=6,h=6)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
abline(h=0.2,v=0.5,col="red",lwd=1.5)
dev.off()
png('06.brown_gene_cor.png',w=600,h=600)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
abline(h=0.2,v=0.5,col="red",lwd=1.5)
dev.off()




