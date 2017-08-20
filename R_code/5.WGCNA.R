#Contant:
#5.1 Network analysis of expression data from XX:
#       5.1.1 Data input and cleaning 
#       5.1.2 Network construction and module detection
#               5.1.2a. Automatic, one-step network construction and module detection
#               5.1.2b Step-by-step network construction and module detection
#               5.1.2c Dealing with large datasets: block-wise network construction and module detection
#       5.1.3 Relating modules to external clinical traits and identifying important genes
#       5.1.4 Interfacing network analysis with other data such as functional annotation and gene ontology  (required)=============================
#       5.1.5 Network visualization using WGCNA functions
#       5.1.6 Export of networks to external software
#5.2 Consensus analysis of naive state Stem Cell and Danwei's expression data

########################################################################
#
#  0 preparation and parameter adjustion
# 
# ######################################################################

list.of.cran.packages<- c("easypackages")
new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(easypackages)

source("https://bioconductor.org/biocLite.R")
list.of.bio.packages <- c("WGCNA")
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)

libraries(list.of.bio.packages,list.of.cran.packages)

# Display the current working directory=======================
getwd()

if (Sys.info()[['sysname']]=="Darwin"){
        setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.

if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) disableWGCNAThreads()
if (is.na(Sys.getenv("RSTUDIO", unset = NA)))  enableWGCNAThreads()

#extract the top 5000 most variant genes for WGCNA studies.
#set up the extract
extract <- TRUE
#extract <- FALSE

########################################################################################
#
# 5.1 Network analysis of expression data from XX:
#  finding modules related to Cell.Stage
#
#  https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
#  labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html
########################################################################################
#====5.1.1 Data input and cleaning (required)=============================

#=====================================================================================
#
#  Code chunk5.1.1-1
#
#=====================================================================================

# Read in rnaseq_MicroArrayData
Danwei_Vassena = read.csv("rnaseq_MicroArrayData.csv",row.names = 1)
# Take a quick look at what is in the data set:
dim(Danwei_Vassena)
names(Danwei_Vassena)


#=====================================================================================
#
#  Code chunk5.1.1-2
#
#=====================================================================================


datExpr0 = as.data.frame(t(Danwei_Vassena[, -c(1:9)]))
head(names(datExpr0))
rownames(datExpr0)

ProbeNames = rownames(Danwei_Vassena)
Danwei.Data = data.frame(t(Danwei_Vassena[,1:9]))
Vassena.Data = data.frame(t(Danwei_Vassena[,10:33]))
Subsets = c(rep(2, times=9), rep(1, times=23))


#=====================================================================================
#
#  Code chunk5.1.1-3
#
#=====================================================================================


gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK


#-------------------------------------------------------------------------
#
#  Code chunk5.1.1-4
#
#-------------------------------------------------------------------------


if (!gsg$allOK)
{
        # Optionally, print the gene and sample names that were removed:
        if (sum(!gsg$goodGenes)>0) 
                printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
        if (sum(!gsg$goodSamples)>0) 
                printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
        # Remove the offending genes and samples from the data:
        datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


#=====================================================================================
#
#  Code chunk5.1.1-5
#
#=====================================================================================
#replace all datExpr0 with Vassena.Data

sampleTree = hclust(dist(Vassena.Data), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.

#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9)
par(cex = 0.6)
par(mfrow = c(1,1))
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


#=====================================================================================
#
#  Code chunk5.1.1-6
#
#=====================================================================================


# Plot a line to show the cut
abline(h = 160, col = "red")
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 160, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = Vassena.Data[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#extract the top 5000 most variant genes for WGCNA studies.
#transpose matrix to correlate genes in the following
if(extract==TRUE){datExpr = datExpr[,order(apply(datExpr,2,mad), decreasing = T)[1:5000]]}


#-------------------------------------------------------------------------
#
#  Code chunk5.1.1-7
#
#-------------------------------------------------------------------------


traitData = read.csv("ClinicalTraits.csv",row.names = 1)
dim(traitData)
names(traitData)

# remove columns that hold information we do not need.
# convert character to numeric in r
allTraits = apply(traitData,2,function(x) as.numeric(as.factor(x)))
allTraits <- as.data.frame(allTraits)
names(allTraits)
rownames(allTraits) <- traitData$Cell
head(allTraits,10)
# Form a data frame analogous to expression data that will hold the clinical traits.

#replace all femaleSamples with Vassena.Samples


Vassena.Samples = rownames(datExpr)
traitRows = match(Vassena.Samples, traitData$Cell)
datTraits = allTraits[traitRows, ]
#rownames(datTraits) = allTraits[traitRows, 1]

collectGarbage()

#-------------------------------------------------------------------------
#
#  Code chunk5.1.1-8
#
#-------------------------------------------------------------------------


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)
#Error in numbers2colors(datTraits[, "Cell.Stage"], signed = FALSE) :
#'x' must be numeric. For a factor, please use as.numeric(x) in the call.

# Plot the sample dendrogram and the colors underneath.
par(mfrow = c(1,1))
par(mar = c(0,5,0,0))
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")


#=====================================================================================
#
#  Code chunk5.1.1-9
#
#=====================================================================================


save(datExpr, datTraits, file = "Vassena-01-dataInput.RData")

#====5.1.2 Network construction and module detection (required)=============================

#====5.1.2a. Automatic, one-step network construction and module detection

#=====================================================================================
#  Code chunk5.1.2a-1
#
#=====================================================================================

# Allow multi-threading within WGCNA. This helps speed up certain calculations.
# At present this call is necessary for the code to work.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments. 
# See note above.

if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) disableWGCNAThreads()
if (is.na(Sys.getenv("RSTUDIO", unset = NA)))  enableWGCNAThreads()
# enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "Vassena-01-dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames
#extract the top 5000 most variant genes for WGCNA studies.
#transpose matrix to correlate genes in the following
if(extract==TRUE){datExpr = datExpr[,order(apply(datExpr,2,mad), decreasing = T)[1:5000]]}
#=====================================================================================
#
#  Code chunk5.1.2a-2
# labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html
#=====================================================================================


# Choose a set of soft-thresholding powers ?????????????????????
powers = c(c(1:10), seq(from = 12, to=20, by=2)) 
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5) #takes long time
# Plot the results:
par(mar=c(2,4,2,4))
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.60,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#chose power for Unsigned and signed hybrid networks networks 
#https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html
cuts <- c(0, 20, 30, 40,60, Inf)
chose.power <- c(10, 9,8,7,6)
power <- chose.power[findInterval(nrow(datExpr), cuts)]
power

#=====================================================================================
#
#  Code chunk5.1.2a-3
# 
#=====================================================================================

# Constructing the gene network and identifying modules is now a simple function call:
#?????????????????????
net = blockwiseModules(datExpr, 
                       power = power,#power = 6 soft-thresholding power for network construction.
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "Vassena.TOM", 
                       verbose = 3) #takes long time
# To see how many modules were identified and what the module sizes are
table(net$colors)

#=====================================================================================
#
#  Code chunk5.1.2a-4
#
#=====================================================================================


# ??????????????????

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#=====================================================================================
#
#  Code chunk5.1.2a-5
#
#=====================================================================================


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "Vassena-02-networkConstruction-auto.RData")

#====> # jump to chunk5.1.6-4


#=====5.1.2b Step-by-step network construction and module detection======
#=====================================================================================
#
#  Code chunk5.1.2b-1
#
#=====================================================================================

# Load the data saved in the first part
lnames = load(file = "Vassena-01-dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames
#extract the top 5000 most variant genes for WGCNA studies.
#transpose matrix to correlate genes in the following
if(extract==TRUE){datExpr = datExpr[,order(apply(datExpr,2,mad), decreasing = T)[1:5000]]}

#=====================================================================================
#
#  Code chunk5.1.2b-2
# labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:

par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.7,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#chose power for Unsigned and signed hybrid networks networks 
#https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html
cuts <- c(0, 20, 30, 40,60, Inf)
chose.power <- c(10, 9,8,7,6)
power <- chose.power[findInterval(nrow(datExpr), cuts)]
power
#=====================================================================================
#
#  Code chunk5.1.2b-3
#
#=====================================================================================
#generated similarity matrices S based on Pearson correlations between all gene pairs"
S <- 0.5 + 0.5*cor(datExpr, method="spearman") #Sij = 0.5 + 0.5 × cor(i,j)


#????????????:(1) Co-expression similarity and adjacency
#built adjacency matrices for each dataset using a soft power threshold of 60
adjacency = adjacency.fromSimilarity(S, power = power,type = "signed") #type = "signed"! #power =softPower
SubGeneNames <-rownames(adjacency)

#=====================================================================================
#
#  Code chunk5.1.2b-4f
# https://www.researchgate.net/post/What_do_adjacency_matrix_and_Topology_Overlap_Matrix_from_WGCNA_package_tell_about_the_data
#=====================================================================================
# To minimize effects of noise and spurious associations
# ????????????????????????????????????,Turn adjacency into topological overlap
# calculate the corresponding dissimilarity
TOM = TOMsimilarity(adjacency)
rownames(TOM) <- SubGeneNames
colnames(TOM) <- SubGeneNames #same like abov
dissTOM = 1-TOM

#=====================================================================================
#
#  Code chunk5.1.2b-5
#
#=====================================================================================

#??????????????????
# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
# Plot the resulting clustering tree (dendrogram)
par(mfrow = c(1,1))
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)


#=====================================================================================
#
#  Code chunk5.1.2b-6
#
#=====================================================================================

#?????????????????????dynamicTreeCut
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)


#=====================================================================================
#
#  Code chunk5.1.2b-7
#
#=====================================================================================

#??????????????????
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


#=====================================================================================
#
#  Code chunk5.1.2b-8
#
#=====================================================================================

#?????????????????????????????????,Merging of modules whose expression profiles are very similar
#?????????????????????leaf???????????????,??????????????????,
#?????????????????????????????????????????????????????????,???????????????????????????modules????????????
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
#moduleEigengenes represents the module expressions of the q-th module by the module eigengene E
#The eigengene E can be thought of as a weighted average expression profile.

MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result

plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")


#=====================================================================================
#
#  Code chunk5.1.2b-9
#
#=====================================================================================

#?????????75%????????????????????????
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs


#=====================================================================================
#
#  Code chunk5.1.2b-10
#
#=====================================================================================


#???????????????(Dynamic Tree Cut)????????????(Merged dynamic)????????????
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()
#??????????????????????????????
#plotDendroAndColors(geneTree,mergedColors,"Merged dynamic",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05)
#=====================================================================================
#
#  Code chunk5.1.2b-11
#
#=====================================================================================


# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "Vassena-02-networkConstruction-stepByStep.RData")

#Extract modules
module_colors <- setdiff(unique(mergedColors), "grey")
for (color in module_colors){
        module=SubGeneNames[which(mergedColors==color)]
        write.table(module, paste("Vassena.module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
        
}
#Look at expression patterns of these genes, as they are clustered

module.order <- unlist(tapply(1:ncol(datExpr),as.factor(mergedColors),I)) #dynamicColors
m<-t(t(datExpr[,module.order])/apply(datExpr[,module.order],2,max))
par(mfrow = c(1,1))
par(mar = c(7,2,1,1))
cex1 = 0.9
heatmap(t(m),zlim=c(0,1),col=gray.colors(100),Rowv=NA,Colv=NA,labRow=NA,scale="none",RowSideColors=mergedColors[module.order])
#===> jump to chunk5.1.6-2

#===5.1.2c Dealing with large datasets: block-wise network construction and module detection====

#=====================================================================================
#
#  Code chunk5.1.2c-1
#
#=====================================================================================

# Load the data saved in the first part
lnames = load(file = "Vassena-01-dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames
#extract the top 5000 most variant genes for WGCNA studies.
#transpose matrix to correlate genes in the following
#datExpr = datExpr[,order(apply(datExpr,2,mad), decreasing = T)[1:5000]]

#=====================================================================================
#
#  Code chunk5.1.2c-2
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:

par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#chose power for Unsigned and signed hybrid networks networks 
#https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html
cuts <- c(0, 20, 30, 40,60, Inf)
chose.power <- c(10, 9,8,7,6)
power <- chose.power[findInterval(nrow(datExpr), cuts)]
power

#=====================================================================================
#
#  Code chunk5.1.2c-3
#
#=====================================================================================


bwnet = blockwiseModules(datExpr, maxBlockSize = 2000,
                         power = power, TOMType = "signed", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "Vassena.TOM-blockwise",
                         verbose = 3)


#=====================================================================================
#
#  Code chunk5.1.2c-4
#
#=====================================================================================


# Load the results of single-block analysis
load(file = "Vassena-02-networkConstruction-auto.RData") #number doesn't match
# Relabel blockwise modules
bwLabels = matchLabels(bwnet$colors, moduleLabels)
# Convert labels to colors for plotting
bwModuleColors = labels2colors(bwLabels)
#To see how many modules were identified and what the module sizes are
table(bwLabels)
#=====================================================================================
#
#  Code chunk5.1.2c-5
#
#=====================================================================================


# open a graphics window

# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Plot the dendrogram and the module colors underneath for block 2
plotDendroAndColors(bwnet$dendrograms[[2]], bwModuleColors[bwnet$blockGenes[[2]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 2", 
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#=====================================================================================
#
#  Code chunk5.1.2c-6
#
#=====================================================================================



plotDendroAndColors(geneTree,
                    cbind(moduleColors, bwModuleColors),
                    c("Single block", "2 blocks"),
                    main = "Single block gene dendrogram and module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#=====================================================================================
#
#  Code chunk5.1.2c-7
#
#=====================================================================================


singleBlockMEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
blockwiseMEs = moduleEigengenes(datExpr, bwModuleColors)$eigengenes


#=====================================================================================
#
#  Code chunk5.1.2c-8
#
#=====================================================================================


single2blockwise = match(names(singleBlockMEs), names(blockwiseMEs))
signif(diag(cor(blockwiseMEs[, single2blockwise], singleBlockMEs)), 3)


#====5.1.3 Relating modules to external clinical traits and identifying important genes (required)=============================
#=====================================================================================
#
#  Code chunk5.1.3-1
#
#=====================================================================================

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Load the expression and trait data saved in the first part
lnames = load(file = "Vassena-01-dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Vassena-02-networkConstruction-stepByStep.RData")
lnames


#=====================================================================================
#
#  Code chunk5.1.3-2
#
#=====================================================================================


# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p") #method = "spearman"
#Remove columns from dataframe where some of values are NA
moduleTraitCor <- moduleTraitCor[ , apply(moduleTraitCor, 2, function(x) !any(is.na(x)))]

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


#=====================================================================================
#
#  Code chunk5.1.3-3
#
#=====================================================================================



# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mfrow =c(1,1))
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
#WGCNA::greenWhiteRed: this palette is not suitable for people
#with green-red color blindness (the most common kind of color blindness).
#Consider using the function blueWhiteRed instead.

#=====================================================================================
#
#  Code chunk5.1.3-4
#
#=====================================================================================


# Define variable Cell.Stage containing the Cell.Stage column of datTrait
Cell.Stage = as.data.frame(datTraits$Cell.Stage)
names(Cell.Stage) = "Cell.Stage"
#names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(datExpr, Cell.Stage, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(Cell.Stage), sep="")
names(GSPvalue) = paste("p.GS.", names(Cell.Stage), sep="")

#=====================================================================================
#
#  Code chunk5.1.3-5
#
#=====================================================================================


module = "turquoise" # based on Module-trait relationships
column = match(module, modNames)
moduleGenes = moduleColors==module


par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for cell stage",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "turquoise") #yellow

#=====================================================================================
#
#  Code chunk5.1.3-6
#
#=====================================================================================
names(datExpr)

#=====================================================================================
#
#  Code chunk5.1.3-7
#
#=====================================================================================
names(datExpr)[moduleColors=="turquoise"]

#---------------------------------------------------------------#
#
# Code chunk5.1.3-8
#
#---------------------------------------------------------------#
#Skip this 

#annot = read.csv(file = "GeneAnnotation.csv")
#dim(annot)
#names(annot)
#probes = names(datExpr)
#probes2annot = match(probes, annot$substanceBXH)
# The following is the number or probes without annotation:
#sum(is.na(probes2annot))
# Should return 0.


#=====================================================================================
#
#  Code chunk5.1.3-9
#
#=====================================================================================
# We now create a data frame holding the following information for all probes
# Create the starting data frame
geneInfo0 = data.frame(
        #     substanceBXH = probes,
        geneSymbol = colnames(datExpr),
        #     LocusLinkID = annot$LocusLinkID[probes2annot],
        moduleColor = moduleColors,
        geneTraitSignificance,
        GSPvalue)
# Order modules by their significance for Cell.Stage
modOrder = order(-abs(cor(MEs, Cell.Stage, use = "p")))
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
        oldNames = names(geneInfo0)
        geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                               MMPvalue[, modOrder[mod]])
        names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                             paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
head(geneInfo0[1:3,])
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Cell.Stage))
geneInfo = geneInfo0[geneOrder, ]


#=====================================================================================
#
#  Code chunk5.1.3-10
#
#=====================================================================================
write.csv(geneInfo, file = "geneInfo.csv")

#====5.1.4 Interfacing network analysis with other data such as functional annotation and gene ontology  (required)=============================

#=====================================================================================
#
#  Code chunk5.1.4-1
#
#=====================================================================================

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Load the expression and trait data saved in the first part
lnames = load(file = "Vassena-01-dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Vassena-02-networkConstruction-stepByStep.RData")
lnames


#=====================================================================================
#
#  Code chunk5.1.4-2
#  consider Code chunk5.1.3-8
#=====================================================================================
allgenes = names(datExpr)
#library(biomaRt)
# downlaod annotationfrom biomAT
#mart <- useMart("ensembl")
#ensembl <- useDataset("hsapiens_gene_ensembl", mart)
#annot <-getBM(attributes = c("ensembl_gene_id", #Gene stable ID
#                     "external_gene_name", #Gene name
#                     "go_id",   #GO term accession
#                     "affy_hugene_1_0_st_v1"), #AFFY HuGene 1 0 st v1
#      filters = "external_gene_name", values = allgenes,  #save time with filter
#      mart = ensembl)

# Read in the annotation
annot = read.csv(file = "Human_GeneAnnotation.csv")
# Match probes in the data set to the gene IDs in the annotation file 
gene2annot = match(allgenes, annot$Gene.name)
gene2annot <- gene2annot[!is.na(gene2annot)] #Remove all NA values from a vector
# Get the corresponding GO IDs
allLLIDs = annot$GO.term.accession[gene2annot] 
allLLIDs <- allLLIDs[allLLIDs !=""] #Remove the empty element in the vector
# $ Choose interesting modules
#intModules = c("brown", "red", "salmon")
intModules = substring(names(MEs), 3)
for (module in intModules)
{
        # Select module probes
        modGenes = (moduleColors==module)
        # Get their entrez ID codes
        modLLIDs = allLLIDs[modGenes]
        # Write them into a file
        fileName = paste("Vassena.ensembl-", module, ".txt", sep="")
        write.table(as.data.frame(modLLIDs), file = fileName,
                    row.names = FALSE, col.names = FALSE)
}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("Vassena.ensembl-all.txt", sep="")
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)


#=====================================================================================
#
#  Code chunk5.1.4-3
#
#=====================================================================================

#GOenrichmentAnalysis takes long time
GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "human", nBestP = 10)
# This function is deprecated and will be removed in the near future. 

#=====================================================================================
#
#  Code chunk5.1.4-4
#
#=====================================================================================
tab = GOenr$bestPTerms[[4]]$enrichment


#=====================================================================================
#
#  Code chunk5.1.4-5
#
#=====================================================================================
names(tab)


#=====================================================================================
#
#  Code chunk5.1.4-6
#
#=====================================================================================
write.table(tab, file = "GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)


#=====================================================================================
#
#  Code chunk5.1.4-7
#
#=====================================================================================


keepCols = c(1, 2, 5, 6, 7, 12, 13)
screenTab = tab[, keepCols]
# Round the numeric columns to 2 decimal places:
numCols = c(3, 4)
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTab[, 7] = substring(screenTab[, 7], 1, 40)
# Shorten the column names:
colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name")
rownames(screenTab) = NULL
# Set the width of R's output. The reader should play with this number to obtain satisfactory output.
options(width=95)
# Finally, display the enrichment table:
screenTab

#===5.1.5 Network visualization using WGCNA functions===========


#=====================================================================================
#
#  Code chunk5.1.5-1
#
#=====================================================================================
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Load the expression and trait data saved in the first part
lnames = load(file = "Vassena-01-dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Vassena-02-networkConstruction-auto.RData")
#or----------------
lnames = load(file = "Vassena-02-networkConstruction-stepByStep.RData")

lnames
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


#=====================================================================================
#
#  Code chunk5.1.5-2
#
#=====================================================================================

# ???????????????????????????
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = power)
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA

# Call the plot function

TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes") #takes long time
# WARNING: On some computers, this code can take a while to run (20 minutes??). 
# I suggest you skip it.

#=====================================================================================
#
#  Code chunk5.1.5-3
#
#=====================================================================================

#????????????400?????????????????????
nSelect = 400
# For reproducibility, we set the random seed
set.seed(10)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]


# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7
diag(plotDiss) = NA
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")


#=====================================================================================
#
#  Code chunk5.1.5-4
#
#=====================================================================================


# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
Cell.Stage = as.data.frame(datTraits$Cell.Stage)
names(Cell.Stage) = "Cell.Stage"
# Add the cell.stage to existing module eigengenes
MET = orderMEs(cbind(MEs, Cell.Stage))
MET = orderMEs(MEs)
# Plot the relationships among the eigengenes and the trait

par(cex = 0.9)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, 
                      xLabelsAngle = 90)


#=====================================================================================
#
#  Code chunk5.1.5-5
#
#=====================================================================================


# Plot the dendrogram

par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


#===5.1.6 Export of networks to external software=============

#=====================================================================================
#
#  Code chunk5.1.6-1
#
#=====================================================================================

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Load the expression and trait data saved in the first part
lnames = load(file = "Vassena-01-dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Vassena-02-networkConstruction-auto.RData")
lnames


#=====================================================================================
#
#  Code chunk5.1.6-2
#
#=====================================================================================
#chose soft power threshold
cuts <- c(0, 20, 30, 40,60, Inf)
chose.power <- c(10, 9,8,7,6)
power <- chose.power[findInterval(nrow(datExpr), cuts)]
power

# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power = softPower)
# Read in the annotation file
annot = read.csv(file = "GeneAnnotation.csv")
# Match probes in the data set to the gene IDs in the annotation file 
gene2annot = match(allgenes, annot$Gene.name)
gene2annot <- gene2annot[!is.na(gene2annot)] #Remove all NA values from a vector
# Get the corresponding GO IDs
allLLIDs = annot$GO.term.accession[gene2annot] 
allLLIDs <- allLLIDs[allLLIDs !=""] #Remove the empty element in the vector

# Select module
module = "turquoise"
# Select module probes
probes = names(datExpr)
inModule = (moduleColors==module)
modProbes = probes[inModule]
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0)
#                            probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )


#=====================================================================================
#
#  Code chunk5.1.6-3
#
#=====================================================================================


nTop = 30
IMConn = softConnectivity(datExpr[, modProbes])
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("VisANTInput-", module, "-top30.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )


#=====================================================================================
#
#  Code chunk5.1.6-4
#
#=====================================================================================
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 6)
# Read in the annotation file
#annot = read.csv(file = "GeneAnnotation.csv")
# Select modules
modules = c("brown", "red")
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules))
modProbes = probes[inModule]
#modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)]
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               #altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])
#===> jump to chunk5.1.5-1

########################################################################################
#
# 5.2 Consensus analysis of naive state Stem Cell and Danwei's expression data
#  https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
#
########################################################################################


#=====================================================================================
#
#  Code chunk5.2.1-1
#
#=====================================================================================
# Load the package

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
#Read in the female liver data set
Danwei_Vassena = read.csv("LiverFemale3600.csv")
# Read in the male liver data set
maleData = read.csv("LiverMale3600.csv")
# Take a quick look at what is in the data sets (caution, longish output):
dim(Danwei_Vassena)
names(Danwei_Vassena)
dim(maleData)
names(maleData)


#=====================================================================================
#
#  Code chunk5.2.1-2
#
#=====================================================================================


# We work with two sets:
nSets = 2
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("Female liver", "Male liver")
shortLabels = c("Female", "Male")
# Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = vector(mode = "list", length = nSets)

multiExpr[[1]] = list(data = as.data.frame(t(Danwei_Vassena[-c(1:8)])))
names(multiExpr[[1]]$data) = Danwei_Vassena$substanceBXH
rownames(multiExpr[[1]]$data) = names(Danwei_Vassena)[-c(1:8)]
multiExpr[[2]] = list(data = as.data.frame(t(maleData[-c(1:8)])))
names(multiExpr[[2]]$data) = maleData$substanceBXH
rownames(multiExpr[[2]]$data) = names(maleData)[-c(1:8)]
# Check that the data has the correct format for many functions operating on multiple sets:
exprSize = checkSets(multiExpr)


#=====================================================================================
#
#  Code chunk5.2.1-3
#
#=====================================================================================


# Check that all genes and samples have sufficiently low numbers of missing values.
gsg = goodSamplesGenesMS(multiExpr, verbose = 3)
gsg$allOK


#=====================================================================================
#
#  Code chunk5.2.1-4
#
#=====================================================================================


if (!gsg$allOK)
{
        # Print information about the removed genes:
        if (sum(!gsg$goodGenes) > 0)
                printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes], 
                                                          collapse = ", ")))
        for (set in 1:exprSize$nSets)
        {
                if (sum(!gsg$goodSamples[[set]]))
                        printFlush(paste("In set", setLabels[set], "removing samples",
                                         paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
                # Remove the offending genes and samples
                multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes]
        }
        # Update exprSize
        exprSize = checkSets(multiExpr)
}


#=====================================================================================
#
#  Code chunk5.2.1-5
#
#=====================================================================================


sampleTrees = list()
for (set in 1:nSets)
{
        sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}


#=====================================================================================
#
#  Code chunk5.2.1-6
#
#=====================================================================================


#skip pdf(file = "Plots/SampleClustering.pdf", width = 12, height = 12)
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
        plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
             xlab="", sub="", cex = 0.7)
#skip dev.off()


#=====================================================================================
#
#  Code chunk5.2.1-7
#
#=====================================================================================


# Choose the "base" cut height for the female data set
baseHeight = 16
# Adjust the cut height for the male data set for the number of samples
cutHeights = c(16, 16*exprSize$nSamples[2]/exprSize$nSamples[1])
# Re-plot the dendrograms including the cut lines
#skip pdf(file = "Plots/SampleClustering.pdf", width = 12, height = 12)
par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets)
{
        plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
             xlab="", sub="", cex = 0.7)
        abline(h=cutHeights[set], col = "red")
}
#skip dev.off()


#=====================================================================================
#
#  Code chunk5.2.1-8
#
#=====================================================================================


for (set in 1:nSets)
{
        # Find clusters cut by the line
        labels = cutreeStatic(sampleTrees[[set]], cutHeight = cutHeights[set])
        # Keep the largest one (labeled by the number 1)
        keep = (labels==1)
        multiExpr[[set]]$data = multiExpr[[set]]$data[keep, ]
}
collectGarbage()
# Check the size of the leftover data
exprSize = checkSets(multiExpr)
exprSize


#=====================================================================================
#
#  Code chunk5.2.1-9
#
#=====================================================================================


traitData = read.csv("ClinicalTraits.csv")
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
allTraits = traitData[, -c(31, 16)]
allTraits = allTraits[, c(2, 11:36) ]
# See how big the traits are and what are the trait and sample names
dim(allTraits)
names(allTraits)
allTraits$Mice
# Form a multi-set structure that will hold the clinical traits.
Traits = vector(mode="list", length = nSets)
for (set in 1:nSets)
{
        setSamples = rownames(multiExpr[[set]]$data)
        traitRows = match(setSamples, allTraits$Mice)
        Traits[[set]] = list(data = allTraits[traitRows, -1])
        rownames(Traits[[set]]$data) = allTraits[traitRows, 1]
}
collectGarbage()
# Define data set dimensions
nGenes = exprSize$nGenes
nSamples = exprSize$nSamples


#=====================================================================================
#
#  Code chunk5.2.1-10
#
#=====================================================================================


save(multiExpr, Traits, nGenes, nSamples, setLabels, shortLabels, exprSize, 
     file = "Consensus-dataInput.RData")

#5.2.2 Network construction and consensus module detection============
# a.Automatic, one-step network construction and consensus module detection:
#=====================================================================================
#
#  Code chunk5.2.2a-1
#
#=====================================================================================

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Allow multi-threading within WGCNA. 
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
# enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "Consensus-dataInput.RData")
# The variable lnames contains the names of loaded variables.
lnames
# Get the number of sets in the multiExpr structure.
nSets = checkSets(multiExpr)$nSets


#=====================================================================================
#
#  Code chunk5.2.2a-2
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,20, by=2))
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets)
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
        powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                           verbose = 2)[[2]])
collectGarbage()
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity")
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4)
for (set in 1:nSets)
{
        for (col in 1:length(plotCols))
        {
                ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
                ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
        }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
pdf(file = "Plots/scaleFreeAnalysis.pdf", wi = 8, he = 6)
par(mfcol = c(2,2))
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
        if (set==1)
        {
                plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
                     xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
                     main = colNames[col])
                addGrid()
        }
        if (col==1)
        {
                text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
                     labels=powers,cex=cex1,col=colors[set])
        } else
                text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
                     labels=powers,cex=cex1,col=colors[set])
        if (col==1)
        {
                legend("bottomright", legend = setLabels, col = colors, pch = 20) 
        } else
                legend("topright", legend = setLabels, col = colors, pch = 20) 
}
dev.off()


#=====================================================================================
#
#  Code chunk5.2.2a-3
#
#=====================================================================================


net = blockwiseConsensusModules(
        multiExpr, power = 6, minModuleSize = 30, deepSplit = 2,
        pamRespectsDendro = FALSE, 
        mergeCutHeight = 0.25, numericLabels = TRUE,
        minKMEtoStay = 0,
        saveTOMs = TRUE, verbose = 5)


#=====================================================================================
#
#  Code chunk5.2.2a-4
#
#=====================================================================================


consMEs = net$multiMEs
moduleLabels = net$colors
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
consTree = net$dendrograms[[1]] 


#=====================================================================================
#
#  Code chunk5.2.2a-5
#
#=====================================================================================



pdf(file = "Plots/ConsensusDendrogram-auto.pdf", wi = 8, he = 6)
plotDendroAndColors(consTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")

dev.off()


#=====================================================================================
#
#  Code chunk5.2.2a-6
#
#=====================================================================================


save(consMEs, moduleLabels, moduleColors, consTree, file = "Consensus-NetworkConstruction-auto.RData")


#5.2.2b Step-by-step network construction and module detection, including scaling of Topological Overlap Matrices

#=====================================================================================
#
#  Code chunk5.2.2b-1
#
#=====================================================================================
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Allow multi-threading within WGCNA. 
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.

# Load the data saved in the first part
lnames = load(file = "Consensus-dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames
# Get the number of sets in the multiExpr structure.
nSets = checkSets(multiExpr)$nSets


#=====================================================================================
#
#  Code chunk5.2.2b-2
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,20, by=2))
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets)
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
        powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                           verbose = 2)[[2]])
collectGarbage()
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity")
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4)
for (set in 1:nSets)
{
        for (col in 1:length(plotCols))
        {
                ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
                ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
        }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power

par(mfcol = c(2,2))
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
        if (set==1)
        {
                plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
                     xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
                     main = colNames[col])
                addGrid()
        }
        if (col==1)
        {
                text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
                     labels=powers,cex=cex1,col=colors[set])
        } else
                text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
                     labels=powers,cex=cex1,col=colors[set])
        if (col==1)
        {
                legend("bottomright", legend = setLabels, col = colors, pch = 20) 
        } else
                legend("topright", legend = setLabels, col = colors, pch = 20) 
}


#=====================================================================================
#
#  Code chunk5.2.2b-3
#
#=====================================================================================


softPower = 6
# Initialize an appropriate array to hold the adjacencies
adjacencies = array(0, dim = c(nSets, nGenes, nGenes))
# Calculate adjacencies in each individual data set
for (set in 1:nSets)
        adjacencies[set, , ] = abs(cor(multiExpr[[set]]$data, use = "p"))^softPower


#=====================================================================================
#
#  Code chunk5.2.2b-4
#
#=====================================================================================


# Initialize an appropriate array to hold the TOMs
TOM = array(0, dim = c(nSets, nGenes, nGenes))
# Calculate TOMs in each individual data set
for (set in 1:nSets)
        TOM[set, , ] = TOMsimilarity(adjacencies[set, , ])


#=====================================================================================
#
#  Code chunk5.2.2b-5
#
#=====================================================================================


# Define the reference percentile
scaleP = 0.95
# Set RNG seed for reproducibility of sampling
set.seed(12345)
# Sample sufficiently large number of TOM entries
nSamples = as.integer(1/(1-scaleP) * 1000)
# Choose the sampled TOM entries
scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
TOMScalingSamples = list()
# These are TOM values at reference percentile
scaleQuant = rep(1, nSets)
# Scaling powers to equalize reference TOM values
scalePowers = rep(1, nSets)
# Loop over sets
for (set in 1:nSets)
{
        # Select the sampled TOM entries
        TOMScalingSamples[[set]] = as.dist(TOM[set, , ])[scaleSample]
        # Calculate the 95th percentile
        scaleQuant[set] = quantile(TOMScalingSamples[[set]],
                                   probs = scaleP, type = 8)
        # Scale the male TOM
        if (set>1)
        {
                scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set])
                TOM[set, ,] = TOM[set, ,]^scalePowers[set]
        }
}


#=====================================================================================
#
#  Code chunk5.2.2b-6
#
#=====================================================================================


# For plotting, also scale the sampled TOM entries
scaledTOMSamples = list()
for (set in 1:nSets)
        scaledTOMSamples[[set]] = TOMScalingSamples[[set]]^scalePowers[set]
# Open a suitably sized graphics window

#skip pdf(file = "Plots/TOMScaling-QQPlot.pdf", wi = 6, he = 6)
# qq plot of the unscaled samples
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[2]], plot.it = TRUE, cex = 0.6,
                    xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[2]),
                    main = "Q-Q plot of TOM", pch = 20)
# qq plot of the scaled samples
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[2]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20)
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off()


#=====================================================================================
#
#  Code chunk5.2.2b-7
#
#=====================================================================================


consensusTOM = pmin(TOM[1, , ], TOM[2, , ])


#=====================================================================================
#
#  Code chunk5.2.2b-8
#
#=====================================================================================


# Clustering
consTree = hclust(as.dist(1-consensusTOM), method = "average")
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
# Module identification using dynamic tree cut:
unmergedLabels = cutreeDynamic(dendro = consTree, distM = 1-consensusTOM,
                               deepSplit = 2, cutHeight = 0.995,
                               minClusterSize = minModuleSize,
                               pamRespectsDendro = FALSE )
unmergedColors = labels2colors(unmergedLabels)


#=====================================================================================
#
#  Code chunk5.2.2b-9
#
#=====================================================================================



plotDendroAndColors(consTree, unmergedColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#=====================================================================================
#
#  Code chunk5.2.2b-10
#
#=====================================================================================


# Calculate module eigengenes
unmergedMEs = multiSetMEs(multiExpr, colors = NULL, universalColors = unmergedColors)
# Calculate consensus dissimilarity of consensus module eigengenes
consMEDiss = consensusMEDissimilarity(unmergedMEs)
# Cluster consensus modules
consMETree = hclust(as.dist(consMEDiss), method = "average")
# Plot the result

par(mfrow = c(1,1))
plot(consMETree, main = "Consensus clustering of consensus module eigengenes",
     xlab = "", sub = "")
abline(h=0.25, col = "red")


#=====================================================================================
#
#  Code chunk5.2.2b-11
#
#=====================================================================================


merge = mergeCloseModules(multiExpr, unmergedLabels, cutHeight = 0.25, verbose = 3)


#=====================================================================================
#
#  Code chunk5.2.2b-12
#
#=====================================================================================


# Numeric module labels
moduleLabels = merge$colors
# Convert labels to colors
moduleColors = labels2colors(moduleLabels)
# Eigengenes of the new merged modules:
consMEs = merge$newMEs


#=====================================================================================
#
#  Code chunk5.2.2b-13
#
#=====================================================================================

plotDendroAndColors(consTree, cbind(unmergedColors, moduleColors),
                    c("Unmerged", "Merged"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#=====================================================================================
#
#  Code chunk5.2.2b-14
#
#=====================================================================================


save(consMEs, moduleColors, moduleLabels, consTree, file = "Consensus-NetworkConstruction-man.RData")


#5.2.2c Dealing with large datasets: block-wise network construction and consensus module detection, including comparing the block-wise approach to the standard single-block method


#=====================================================================================
#
#  Code chunk5.2.2c-1
#
#=====================================================================================

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Allow multi-threading within WGCNA. 
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.

# Load the data saved in the first part
lnames = load(file = "Consensus-dataInput.RData")
# The variable lnames contains the names of loaded variables.
lnames
# Get the number of sets in the multiExpr structure.
nSets = checkSets(multiExpr)$nSets


#=====================================================================================
#
#  Code chunk5.2.2c-2
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,20, by=2))
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets)
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
        powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                           verbose = 2)[[2]])
collectGarbage()
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity")
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4)
for (set in 1:nSets)
{
        for (col in 1:length(plotCols))
        {
                ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
                ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE)
        }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power

#pdf(file = "Plots/scaleFreeAnalysis.pdf", wi = 8, he = 6)
par(mfcol = c(2,2))
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
        if (set==1)
        {
                plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
                     xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
                     main = colNames[col])
                addGrid()
        }
        if (col==1)
        {
                text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
                     labels=powers,cex=cex1,col=colors[set])
        } else
                text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
                     labels=powers,cex=cex1,col=colors[set])
        if (col==1)
        {
                legend("bottomright", legend = setLabels, col = colors, pch = 20) 
        } else
                legend("topright", legend = setLabels, col = colors, pch = 20) 
}
dev.off()


#=====================================================================================
#
#  Code chunk5.2.2c-3
#
#=====================================================================================


bnet = blockwiseConsensusModules(
        multiExpr, maxBlockSize = 2000, power = 6, minModuleSize = 30,
        deepSplit = 2, 
        pamRespectsDendro = FALSE, 
        mergeCutHeight = 0.25, numericLabels = TRUE,
        minKMEtoStay = 0,
        saveTOMs = TRUE, verbose = 5)


#=====================================================================================
#
#  Code chunk5.2.2c-4
#
#=====================================================================================


load(file = "Consensus-NetworkConstruction-auto.RData")
bwLabels = matchLabels(bnet$colors, moduleLabels, pThreshold = 1e-7)
bwColors = labels2colors(bwLabels)


#=====================================================================================
#
#  Code chunk5.2.2c-5
#
#=====================================================================================


# Here we show a more flexible way of plotting several trees and colors on one page

pdf(file = "Plots/BlockwiseGeneDendrosAndColors.pdf", wi = 12, he = 6)
# Use the layout function for more involved screen sectioning
layout(matrix(c(1:4), 2, 2), heights = c(0.8, 0.2), widths = c(1,1))
#layout.show(4)
nBlocks = length(bnet$dendrograms)
# Plot the dendrogram and the module colors underneath for each block
for (block in 1:nBlocks)
        plotDendroAndColors(bnet$dendrograms[[block]], moduleColors[bnet$blockGenes[[block]]],
                            "Module colors", 
                            main = paste("Gene dendrogram and module colors in block", block), 
                            dendroLabels = FALSE, hang = 0.03,
                            addGuide = TRUE, guideHang = 0.05,
                            setLayout = FALSE)
dev.off()


#=====================================================================================
#
#  Code chunk5.2.2c-6
#
#=====================================================================================



pdf(file="Plots/SingleDendro-BWColors.pdf", wi = 12, he = 9)
plotDendroAndColors(consTree,
                    cbind(moduleColors, bwColors),
                    c("Single block", "Blockwise"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Single block consensus gene dendrogram and module colors")
dev.off()

#===5.2.3 Relating the consensus modules to female set-specific modules (this section requires the results of Section 2.a of the female turorial)

#=====================================================================================
#
#  Code chunk5.2.3-1
#
#=====================================================================================

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Load the data saved in the first part
lnames = load(file = "Consensus-dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames
# Load the results of network analysis, tutorial part 2.a
lnames = load(file = "Consensus-NetworkConstruction-auto.RData")
lnames


#=====================================================================================
#
#  Code chunk5.2.3-2
#
#=====================================================================================


lnames = load("../Mouse-Female/Vassena-02-networkConstruction-auto.RData")
lnames
# Rename variables to avoid conflicts
femaleLabels = moduleLabels
femaleColors = moduleColors
femaleTree = geneTree
femaleMEs = orderMEs(MEs, greyName = "ME0")


#=====================================================================================
#
#  Code chunk5.2.3-3
#
#=====================================================================================


lnames = load("Consensus-NetworkConstruction-auto.RData")
lnames


#=====================================================================================
#
#  Code chunk5.2.3-4
#
#=====================================================================================


# Isolate the module labels in the order they appear in ordered module eigengenes
femModuleLabels = substring(names(femaleMEs), 3)
consModuleLabels = substring(names(consMEs[[1]]$data), 3)
# Convert the numeric module labels to color labels
femModules = labels2colors(as.numeric(femModuleLabels))
consModules = labels2colors(as.numeric(consModuleLabels))
# Numbers of female and consensus modules
nFemMods = length(femModules)
nConsMods = length(consModules)
# Initialize tables of p-values and of the corresponding counts
pTable = matrix(0, nrow = nFemMods, ncol = nConsMods)
CountTbl = matrix(0, nrow = nFemMods, ncol = nConsMods)
# Execute all pairwaise comparisons
for (fmod in 1:nFemMods)
        for (cmod in 1:nConsMods)
        {
                femMembers = (femaleColors == femModules[fmod])
                consMembers = (moduleColors == consModules[cmod])
                pTable[fmod, cmod] = -log10(fisher.test(femMembers, consMembers, alternative = "greater")$p.value)
                CountTbl[fmod, cmod] = sum(femaleColors == femModules[fmod] & moduleColors ==
                                                   consModules[cmod])
        }


#=====================================================================================
#
#  Code chunk5.2.3-5
#
#=====================================================================================


# Truncate p values smaller than 10^{-50} to 10^{-50} 
pTable[is.infinite(pTable)] = 1.3*max(pTable[is.finite(pTable)])
pTable[pTable>50 ] = 50 
# Marginal counts (really module sizes)
femModTotals = apply(CountTbl, 1, sum)
consModTotals = apply(CountTbl, 2, sum)
# Actual plotting
pdf(file = "Plots/ConsensusVsFemaleModules.pdf", wi = 10, he = 7)
par(mfrow=c(1,1))
par(cex = 1.0)
par(mar=c(8, 10.4, 2.7, 1)+0.3)
# Use function labeledHeatmap to produce the color-coded table with all the trimmings
labeledHeatmap(Matrix = pTable,
               xLabels = paste(" ", consModules),
               yLabels = paste(" ", femModules),
               colorLabels = TRUE,
               xSymbols = paste("Cons ", consModules, ": ", consModTotals, sep=""),
               ySymbols = paste("Fem ", femModules, ": ", femModTotals, sep=""),
               textMatrix = CountTbl,
               colors = greenWhiteRed(100)[50:100],
               main = "Correspondence of Female set-specific and Female-Male consensus modules",
               cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE)
dev.off()

#====5.2.4 Relating consensus module to external microarray sample traits and exporting the results of network analysis

#=====================================================================================
#
#  Code chunk5.2.4-1
#
#=====================================================================================


# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Load the data saved in the first part
lnames = load(file = "Consensus-dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames
# Also load results of network analysis
lnames = load(file = "Consensus-NetworkConstruction-auto.RData")
lnames
exprSize = checkSets(multiExpr)
nSets = exprSize$nSets


#=====================================================================================
#
#  Code chunk5.2.4-2
#
#=====================================================================================


# Set up variables to contain the module-trait correlations
moduleTraitCor = list()
moduleTraitPvalue = list()
# Calculate the correlations
for (set in 1:nSets)
{
        moduleTraitCor[[set]] = cor(consMEs[[set]]$data, Traits[[set]]$data, use = "p")
        moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set])
}


#=====================================================================================
#
#  Code chunk5.2.4-3
#
#=====================================================================================


# Convert numerical lables to colors for labeling of modules in the plot
MEColors = labels2colors(as.numeric(substring(names(consMEs[[1]]$data), 3)))
MEColorNames = paste("ME", MEColors, sep="")
# Open a suitably sized window (the user should change the window size if necessary)

pdf(file = "Plots/ModuleTraitRelationships-female.pdf", wi = 10, he = 7)
# Plot the module-trait relationship table for set number 1
set = 1
textMatrix =  paste(signif(moduleTraitCor[[set]], 2), "\n(",
                    signif(moduleTraitPvalue[[set]], 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor[[set]])
par(mar = c(6, 8.8, 3, 2.2))
labeledHeatmap(Matrix = moduleTraitCor[[set]],
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module--trait relationships in", setLabels[set]))
dev.off()
# Plot the module-trait relationship table for set number 2
set = 2
textMatrix =  paste(signif(moduleTraitCor[[set]], 2), "\n(",
                    signif(moduleTraitPvalue[[set]], 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor[[set]])

pdf(file = "Plots/ModuleTraitRelationships-male.pdf", wi = 10, he = 7)
par(mar = c(6, 8.8, 3, 2.2))
labeledHeatmap(Matrix = moduleTraitCor[[set]],
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module--trait relationships in", setLabels[set]))
dev.off()


#=====================================================================================
#
#  Code chunk5.2.4-4
#
#=====================================================================================


# Initialize matrices to hold the consensus correlation and p-value
consensusCor = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]))
consensusPvalue = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]))
# Find consensus negative correlations
negative = moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0
consensusCor[negative] = pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative])
consensusPvalue[negative] = pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative])
# Find consensus positive correlations
positive = moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0
consensusCor[positive] = pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive])
consensusPvalue[positive] = pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive])


#=====================================================================================
#
#  Code chunk5.2.4-5
#
#=====================================================================================


textMatrix =  paste(signif(consensusCor, 2), "\n(",
                    signif(consensusPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor[[set]])

pdf(file = "Plots/ModuleTraitRelationships-consensus.pdf", wi = 10, he = 7)
par(mar = c(6, 8.8, 3, 2.2))
labeledHeatmap(Matrix = consensusCor,
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Consensus module--trait relationships across\n",
                            paste(setLabels, collapse = " and ")))


#=====================================================================================
#
#  Code chunk5.2.4-6
#
#=====================================================================================


file = "GeneAnnotation.csv"
annot = read.csv(file = file)
# Match probes in the data set to the probe IDs in the annotation file 
probes = names(multiExpr[[1]]$data)
probes2annot = match(probes, annot$substanceBXH)


#=====================================================================================
#
#  Code chunk5.2.4-7
#
#=====================================================================================


consMEs.unord = multiSetMEs(multiExpr, universalColors = moduleLabels, excludeGrey = TRUE)
GS = list()
kME = list()
for (set in 1:nSets)
{
        GS[[set]] = corAndPvalue(multiExpr[[set]]$data, Traits[[set]]$data)
        kME[[set]] = corAndPvalue(multiExpr[[set]]$data, consMEs.unord[[set]]$data)
}


#=====================================================================================
#
#  Code chunk5.2.4-8
#
#=====================================================================================


GS.metaZ = (GS[[1]]$Z + GS[[2]]$Z)/sqrt(2)
kME.metaZ = (kME[[1]]$Z + kME[[2]]$Z)/sqrt(2)
GS.metaP = 2*pnorm(abs(GS.metaZ), lower.tail = FALSE)
kME.metaP = 2*pnorm(abs(kME.metaZ), lower.tail = FALSE)


#=====================================================================================
#
#  Code chunk5.2.4-9
#
#=====================================================================================


GSmat = rbind(GS[[1]]$cor, GS[[2]]$cor, GS[[1]]$p, GS[[2]]$p, GS.metaZ, GS.metaP)
nTraits = checkSets(Traits)$nGenes
traitNames = colnames(Traits[[1]]$data)
dim(GSmat) = c(nGenes, 6*nTraits)
rownames(GSmat) = probes
colnames(GSmat) = spaste(
        c("GS.set1.", "GS.set2.", "p.GS.set1.", "p.GS.set2.", "Z.GS.meta.", "p.GS.meta"),
        rep(traitNames, rep(6, nTraits)))
# Same code for kME:
kMEmat = rbind(kME[[1]]$cor, kME[[2]]$cor, kME[[1]]$p, kME[[2]]$p, kME.metaZ, kME.metaP)
MEnames = colnames(consMEs.unord[[1]]$data)
nMEs = checkSets(consMEs.unord)$nGenes
dim(kMEmat) = c(nGenes, 6*nMEs)
rownames(kMEmat) = probes
colnames(kMEmat) = spaste(
        c("kME.set1.", "kME.set2.", "p.kME.set1.", "p.kME.set2.", "Z.kME.meta.", "p.kME.meta"),
        rep(MEnames, rep(6, nMEs)))


#=====================================================================================
#
#  Code chunk5.2.4-10
#
#=====================================================================================


info = data.frame(Probe = probes, GeneSymbol = annot$gene_symbol[probes2annot],
                  EntrezID = annot$LocusLinkID[probes2annot],
                  ModuleLabel = moduleLabels,
                  ModuleColor = labels2colors(moduleLabels),
                  GSmat,
                  kMEmat)
write.csv(info, file = "consensusAnalysis-CombinedNetworkResults.csv",
          row.names = FALSE, quote = FALSE)

#==5.2.5 Studying and comparing the relationships among modules and traits between the two data sets, including the visualization of consensus eigengene networks and the results of the differential analysis


#=====================================================================================
#
#  Code chunk5.2.5-1
#
#=====================================================================================
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Basic settings: we work with two data sets
nSets = 2
# For easier labeling of plots, create a vector holding descriptive names of the two sets.
setLabels = c("Female liver", "Male liver")
shortLabels = c("Female", "Male")
# Load the data saved in the first part
lnames = load(file = "Consensus-dataInput.RData")
#The variable lnames contains the names of loaded variables.
lnames
# Load the results of network analysis, tutorial part 2.a
lnames = load(file = "Consensus-NetworkConstruction-auto.RData")
lnames


#=====================================================================================
#
#  Code chunk5.2.5-2
#
#=====================================================================================


# Create a variable weight that will hold just the body weight of mice in both sets
weight = vector(mode = "list", length = nSets)
for (set in 1:nSets)
{
        weight[[set]] = list(data = as.data.frame(Traits[[set]]$data$weight_g))
        names(weight[[set]]$data) = "weight"
}
# Recalculate consMEs to give them color names
consMEsC = multiSetMEs(multiExpr, universalColors = moduleColors)
# We add the weight trait to the eigengenes and order them by consesus hierarchical clustering:
MET = consensusOrderMEs(addTraitToMEs(consMEsC, weight))


#=====================================================================================
#
#  Code chunk5.2.5-3
#
#=====================================================================================



pdf(file = "Plots/EigengeneNetworks.pdf", width= 8, height = 10)
par(cex = 0.9)
plotEigengeneNetworks(MET, setLabels, marDendro = c(0,2,2,1), marHeatmap = c(3,3,2,1),
                      zlimPreservation = c(0.5, 1), xLabelsAngle = 90)
dev.off()


########################################################################################
#
# 5.3 Analysis of simulated data
#  https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
#
########################################################################################

#===5.3.1 Simulation of expression and trait data============

#=====================================================================================
#
#  Code chunk5.3.1-1
#
#=====================================================================================
options(stringsAsFactors = FALSE)


#=====================================================================================
#
#  Code chunk5.3.1-2
#
#=====================================================================================


# Here are input parameters of the simulation model
# number of samples or microarrays in the training data
no.obs=50
# now we specify the true measures of eigengene significance
# recall that ESturquoise=cor(y,MEturquoise)
ESturquoise=0   ESbrown= -.6
ESgreen=.6ESyellow=0
# Note that we donate specify the eigengene significance of the blue module
# since it is highly correlated with the turquoise module.
ESvector=c(ESturquoise,ESbrown,ESgreen,ESyellow)
# number of genes 
nGenes1=3000
# proportion of genes in the turquoise, blue, brown, green, and yellow module #respectively.
simulateProportions1=c(0.2,0.15, 0.08, 0.06, 0.04)
# Note that the proportions donate add up to 1. The remaining genes will be colored grey,
# ie the grey genes are non-module genes.
# set the seed of the random number generator. As a homework exercise change this seed.
set.seed(1)
#Step 1: simulate a module eigengene network.
# Training Data Set I
MEgreen=rnorm(no.obs)
scaledy=MEgreen*ESgreen+sqrt(1-ESgreen^2)*rnorm(no.obs)
y=ifelse( scaledy>median(scaledy),2,1)
MEturquoise= ESturquoise*scaledy+sqrt(1-ESturquoise^2)*rnorm(no.obs)
# we simulate a strong dependence between MEblue and MEturquoise
MEblue= .6*MEturquoise+ sqrt(1-.6^2) *rnorm(no.obs)
MEbrown= ESbrown*scaledy+sqrt(1-ESbrown^2)*rnorm(no.obs)
MEyellow= ESyellow*scaledy+sqrt(1-ESyellow^2)*rnorm(no.obs)
ModuleEigengeneNetwork1=data.frame(y,MEturquoise,MEblue,MEbrown,MEgreen, MEyellow)


#=====================================================================================
#
#  Code chunk5.3.1-3
#
#=====================================================================================


dat1=simulateDatExpr5Modules(MEturquoise=ModuleEigengeneNetwork1$MEturquoise,
                             MEblue=ModuleEigengeneNetwork1$MEblue,
                             MEbrown=ModuleEigengeneNetwork1$MEbrown,
                             MEyellow=ModuleEigengeneNetwork1$MEyellow,
                             MEgreen=ModuleEigengeneNetwork1$MEgreen, 
                             nGenes=nGenes1, 
                             simulateProportions=simulateProportions1)


#=====================================================================================
#
#  Code chunk5.3.1-4
#
#=====================================================================================


datExpr = dat1$datExpr
truemodule = dat1$truemodule
datME = dat1$datME
attach(ModuleEigengeneNetwork1)


#=====================================================================================
#
#  Code chunk5.3.1-5
#
#=====================================================================================


table(truemodule)
dim(datExpr)


#=====================================================================================
#
#  Code chunk5.3.1-6
#
#=====================================================================================


datExpr=data.frame(datExpr)
ArrayName=paste("Sample",1:dim(datExpr)[[1]], sep="" )   
# The following code is useful for outputting the simulated data 
GeneName=paste("Gene",1:dim(datExpr)[[2]], sep="" )   
dimnames(datExpr)[[1]]=ArrayName
dimnames(datExpr)[[2]]=GeneName


#=====================================================================================
#
#  Code chunk5.3.1-7
#
#=====================================================================================


rm(dat1) collectGarbage()
# The following command will save all variables defined in the current session.
save.image("Simulated-dataSimulation.RData")



#5.3.2 Loading of expression data, an alternative to data simulation, provided to illustrate data loading of real data


#=====================================================================================
#
#  Code chunk5.3.2-1
#
#=====================================================================================

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)


#=====================================================================================
#
#  Code chunk5.3.2-2
#
#=====================================================================================


datGeneSummary=read.csv("GeneSummaryTutorial.csv")
datTraits=read.csv("TraitsTutorial.csv")
datMicroarrays=read.csv("MicroarrayDataTutorial.csv")


#=====================================================================================
#
#  Code chunk5.3.2-3
#
#=====================================================================================


# This vector contains the microarray sample names
ArrayName= names(data.frame(datMicroarrays[,-1]))
# This vector contains the gene names
GeneName= datMicroarrays$GeneName
# We transpose the data so that the rows correspond to samples and the columns correspond to genes
# Since the first column contains the gene names, we exclude it.
datExpr=data.frame(t(datMicroarrays[,-1]))
names(datExpr)=datMicroarrays[,1]
dimnames(datExpr)[[1]]=names(data.frame(datMicroarrays[,-1]))
#Also, since we simulated the data, we know the true module color:
truemodule= datGeneSummary$truemodule
rm(datMicroarrays)
collectGarbage()


#=====================================================================================
#
#  Code chunk5.3.2-4
#
#=====================================================================================


# First, make sure that the array names in the file datTraits line up with those in the microarray data 
table( dimnames(datExpr)[[1]]==datTraits$ArrayName)
# Next, define the microarray sample trait 
y=datTraits$y


#====5.3.3 Basic data preprocessing illustrates rudimentary techniques for handling missing data and removing outliers

#=====================================================================================
#
#  Code chunk5.3.3-1
#
#=====================================================================================


# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Load the previously saved data
load("Simulated-dataSimulation.RData")
attach(ModuleEigengeneNetwork1)


#=====================================================================================
#
#  Code chunk5.3.3-2
#
#=====================================================================================


meanExpressionByArray=apply( datExpr,1,mean, na.rm=T)  
NumberMissingByArray=apply( is.na(data.frame(datExpr)),1, sum)  


#=====================================================================================
#
#  Code chunk5.3.3-3
#
#=====================================================================================



barplot(meanExpressionByArray,
        xlab = "Sample", ylab = "Mean expression",
        main ="Mean expression across samples",
        names.arg = c(1:50), cex.names = 0.7)


#=====================================================================================
#
#  Code chunk5.3.3-4
#
#=====================================================================================


# Keep only arrays containing less than 500 missing entries
KeepArray= NumberMissingByArray<500
table(KeepArray)
datExpr=datExpr[KeepArray,]
y=y[KeepArray]
ArrayName[KeepArray]


#=====================================================================================
#
#  Code chunk5.3.3-5
#
#=====================================================================================


NumberMissingByGene =apply( is.na(data.frame(datExpr)),2, sum)
# One could do a barplot(NumberMissingByGene), but the barplot is empty in this case.
# It may be better to look at the numbers of missing samples using the summary method:
summary(NumberMissingByGene)
# Calculate the variances of the probes and the number of present entries
variancedatExpr=as.vector(apply(as.matrix(datExpr),2,var, na.rm=T))
no.presentdatExpr=as.vector(apply(!is.na(as.matrix(datExpr)),2, sum) )
# Another way of summarizing the number of pressent entries
table(no.presentdatExpr)
# Keep only genes whose variance is non-zero and have at least 4 present entries
KeepGenes= variancedatExpr>0 & no.presentdatExpr>=4
table(KeepGenes)
datExpr=datExpr[, KeepGenes]
GeneName=GeneName[KeepGenes]


#=====================================================================================
#
#  Code chunk5.3.3-6
#
#=====================================================================================



plotClusterTreeSamples(datExpr=datExpr, y=y)

#==5.3.4 Standard gene screening illustrates gene selection based on Pearson correlation and shows that the results are not satisfactory

#=====================================================================================
#
#  Code chunk5.3.4-1
#
#=====================================================================================

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Load the previously saved data
load("Simulated-dataSimulation.RData")
attach(ModuleEigengeneNetwork1)


#=====================================================================================
#
#  Code chunk5.3.4-2
#
#=====================================================================================


GS1= as.numeric(cor(y, datExpr, use="p"))
# Network terminology: GS1 will be referred to as signed gene significance measure
p.Standard=corPvalueFisher(GS1, nSamples =length(y) )
# since the q-value function has problems with missing data, we use the following trick
p.Standard2=p.Standard
p.Standard2[is.na(p.Standard)]=1
q.Standard=qvalue(p.Standard2)$qvalues
# Form a data frame to hold the results
StandardGeneScreeningResults=data.frame(GeneName,PearsonCorrelation=GS1, p.Standard, q.Standard)


#=====================================================================================
#
#  Code chunk5.3.4-3
#
#=====================================================================================


NoiseGeneIndicator=is.element( truemodule, c("turquoise", "blue", "yellow", "grey"))+.0
SignalGeneIndicator=1-NoiseGeneIndicator


#=====================================================================================
#
#  Code chunk5.3.4-4
#
#=====================================================================================


table(q.Standard<.20)


#=====================================================================================
#
#  Code chunk5.3.4-5
#
#=====================================================================================


mean(NoiseGeneIndicator[q.Standard<=0.20]) 


#=====================================================================================
#
#  Code chunk5.3.4-6
#
#=====================================================================================


save.image(file = "Simulated-StandardScreening.RData")

#==5.3.5 Construction of a weighted gene co-expression network and network modules illustrated step-by-step includes a discussion of alternate clustering techniques

#=====================================================================================
#
#  Code chunk5.3.5-1
#
#=====================================================================================

# Load additional necessary packages
library(cluster)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Load the previously saved data
load("Simulated-StandardScreening.RData") 
attach(ModuleEigengeneNetwork1)


#=====================================================================================
#
#  Code chunk5.3.5-2
#
#=====================================================================================


# here we define the adjacency matrix using soft thresholding with beta=6
ADJ1=abs(cor(datExpr,use="p"))^6
# When you have relatively few genes (<5000) use the following code
k=as.vector(apply(ADJ1,2,sum, na.rm=T))
# When you have a lot of genes use the following code
k=softConnectivity(datE=Vassena.Data,power=6) 
# Plot a histogram of k and a scale free topology plot

par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")


#=====================================================================================
#
#  Code chunk5.3.5-3
#
#=====================================================================================


datExpr=datExpr[, rank(-k,ties.method="first" )<=3600]


#=====================================================================================
#
#  Code chunk5.3.5-4
#
#=====================================================================================


# Turn adjacency into a measure of dissimilarity
dissADJ=1-ADJ1


#=====================================================================================
#
#  Code chunk5.3.5-5
#
#=====================================================================================


dissTOM=TOMdist(ADJ1)
collectGarbage()


#=====================================================================================
#
#  Code chunk5.3.5-6
#
#=====================================================================================


pam4=pam(as.dist(dissADJ), 4)
pam5=pam(as.dist(dissADJ), 5)
pam6=pam(as.dist(dissADJ), 6)
# Cross-tabulte the detected and the true (simulated) module membership:
table(pam4$clustering, truemodule)
table(pam5$clustering, truemodule)
table(pam6$clustering, truemodule)


#=====================================================================================
#
#  Code chunk5.3.5-7
#
#=====================================================================================


pamTOM4=pam(as.dist(dissTOM), 4)
pamTOM5=pam(as.dist(dissTOM), 5)
pamTOM6=pam(as.dist(dissTOM), 6)
# Cross-tabulte the detected and the true (simulated) module membership:
table(pamTOM4$clustering, truemodule)
table(pamTOM5$clustering, truemodule)
table(pamTOM6$clustering, truemodule)


#=====================================================================================
#
#  Code chunk5.3.5-8
#
#=====================================================================================


hierADJ=hclust(as.dist(dissADJ), method="average" )
# Plot the resulting clustering tree together with the true color assignment

plotDendroAndColors(hierADJ, colors = data.frame(truemodule), dendroLabels = FALSE, hang = 0.03, 
                    main = "Gene hierarchical clustering dendrogram and simulated module colors" )


#=====================================================================================
#
#  Code chunk5.3.5-9
#
#=====================================================================================


colorStaticADJ=as.character(cutreeStaticColor(hierADJ, cutHeight=.99, minSize=20))
# Plot the dendrogram with module colors

plotDendroAndColors(hierADJ, colors = data.frame(truemodule, colorStaticADJ),
                    dendroLabels = FALSE, abHeight = 0.99,
                    main = "Gene dendrogram and module colors")


#=====================================================================================
#
#  Code chunk5.3.5-10
#
#=====================================================================================


branch.number=cutreeDynamic(hierADJ,method="tree")
# This function transforms the branch numbers into colors
colorDynamicADJ=labels2colors(branch.number )


#=====================================================================================
#
#  Code chunk5.3.5-11
#
#=====================================================================================


colorDynamicHybridADJ=labels2colors(cutreeDynamic(hierADJ,distM= dissADJ, 
                                                  cutHeight = 0.998, deepSplit=2, pamRespectsDendro = FALSE))

# Plot results of all module detection methods together:

plotDendroAndColors(dendro = hierADJ, 
                    colors=data.frame(truemodule, colorStaticADJ, 
                                      colorDynamicADJ, colorDynamicHybridADJ), 
                    dendroLabels = FALSE, marAll = c(0.2, 8, 2.7, 0.2),
                    main = "Gene dendrogram and module colors")


#=====================================================================================
#
#  Code chunk5.3.5-12
#
#=====================================================================================


# Calculate the dendrogram
hierTOM = hclust(as.dist(dissTOM),method="average")
# The reader should vary the height cut-off parameter h1 
# (related to the y-axis of dendrogram) in the following
colorStaticTOM = as.character(cutreeStaticColor(hierTOM, cutHeight=.99, minSize=20))
colorDynamicTOM = labels2colors (cutreeDynamic(hierTOM,method="tree"))
colorDynamicHybridTOM = labels2colors(cutreeDynamic(hierTOM, distM= dissTOM , cutHeight = 0.998,
                                                    deepSplit=2, pamRespectsDendro = FALSE))
# Now we plot the results

plotDendroAndColors(hierTOM, 
                    colors=data.frame(truemodule, colorStaticTOM, 
                                      colorDynamicTOM, colorDynamicHybridTOM), 
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                    main = "Gene dendrogram and module colors, TOM dissimilarity")


#=====================================================================================
#
#  Code chunk5.3.5-13
#
#=====================================================================================


tabStaticADJ=table(colorStaticADJ,truemodule)
tabStaticTOM=table(colorStaticTOM,truemodule)
tabDynamicADJ=table(colorDynamicADJ, truemodule)
tabDynamicTOM=table(colorDynamicTOM,truemodule)
tabDynamicHybridADJ =table(colorDynamicHybridADJ,truemodule)
tabDynamicHybridTOM =table(colorDynamicHybridTOM,truemodule)


#=====================================================================================
#
#  Code chunk5.3.5-14
#
#=====================================================================================


randIndex(tabStaticADJ,adjust=F)
randIndex(tabStaticTOM,adjust=F)
randIndex(tabDynamicADJ,adjust=F)
randIndex(tabDynamicTOM,adjust=F)
randIndex(tabDynamicHybridADJ ,adjust=F)
randIndex(tabDynamicHybridTOM ,adjust=F)


#=====================================================================================
#
#  Code chunk5.3.5-15
#
#=====================================================================================


colorh1= colorDynamicHybridTOM
# remove the dissimilarities, adjacency matrices etc to free up space
rm(ADJ1) rm(dissADJ)              
collectGarbage()
save.image("Simulated-NetworkConstruction.RData")

#==5.3.6 Relating modules and module eigengenes to external data illustrates methods for relating modules to external microarray sample traits


#=====================================================================================
#
#  Code chunk5.3.6-1
#
#=====================================================================================


# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Load the previously saved data
load("Simulated-NetworkConstruction.RData") 
attach(ModuleEigengeneNetwork1)


#=====================================================================================
#
#  Code chunk5.3.6-2
#
#=====================================================================================


datME=moduleEigengenes(datExpr,colorh1)$eigengenes
signif(cor(datME, use="p"), 2)


#=====================================================================================
#
#  Code chunk5.3.6-3
#
#=====================================================================================


dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
# Plot the eigengene dendrogram
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes")


#=====================================================================================
#
#  Code chunk5.3.6-4
#
#=====================================================================================



plotMEpairs(datME,y=y)


#=====================================================================================
#
#  Code chunk5.3.6-5
#
#=====================================================================================


signif(cor(datME, ModuleEigengeneNetwork1[,-1]),2)


#=====================================================================================
#
#  Code chunk5.3.6-6
#
#=====================================================================================



par(mfrow=c(3,1), mar=c(1, 2, 4, 1))
which.module="turquoise" 
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
# for the second (blue) module we use
which.module="blue"  
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
which.module="brown" 
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )


#=====================================================================================
#
#  Code chunk5.3.6-7
#
#=====================================================================================



which.module="green"
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")


#=====================================================================================
#
#  Code chunk5.3.6-8
#
#=====================================================================================


signif(cor(y,datME, use="p"),2)


#=====================================================================================
#
#  Code chunk5.3.6-9
#
#=====================================================================================


cor.test(y, datME$MEbrown)


#=====================================================================================
#
#  Code chunk5.3.6-10
#
#=====================================================================================


p.values = corPvalueStudent(cor(y,datME, use="p"), nSamples = length(y))


#=====================================================================================
#
#  Code chunk5.3.6-11
#
#=====================================================================================


GS1=as.numeric(cor(y,datExpr, use="p"))
GeneSignificance=abs(GS1)
# Next module significance is defined as average gene significance.
ModuleSignificance=tapply(GeneSignificance, colorh1, mean, na.rm=T)


#=====================================================================================
#
#  Code chunk5.3.6-12
#
#=====================================================================================



par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance,colorh1)


#=====================================================================================
#
#  Code chunk5.3.6-13
#
#=====================================================================================


collectGarbage()
save.image("Simulated-RelatingToExt.RData")

#5.3.7 Module membership, intramodular connectivity, and screening for 
#intramodular hub genes illustrates using the intramodular connectivity 
#to define measures of module membership and to screen for genes based on network information

#=====================================================================================
#
#  Code chunk5.3.7-1
#
#=====================================================================================
library(cluster)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Load the previously saved data
load("Simulated-RelatingToExt.RData") 


#=====================================================================================
#
#  Code chunk5.3.7-2
#
#=====================================================================================


ADJ1=abs(cor(datExpr,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ1, colorh1)
head(Alldegrees1)


#=====================================================================================
#
#  Code chunk5.3.7-3
#
#=====================================================================================


colorlevels=unique(colorh1)

par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))
for (i in c(1:length(colorlevels))) 
{
        whichmodule=colorlevels[[i]] 
        restrict1 = (colorh1==whichmodule)
        verboseScatterplot(Alldegrees1$kWithin[restrict1], 
                           GeneSignificance[restrict1], col=colorh1[restrict1],
                           main=whichmodule, 
                           xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}


#=====================================================================================
#
#  Code chunk5.3.7-4
#
#=====================================================================================


datKME=signedKME(datExpr, datME, outputColumnName="MM.")
# Display the first few rows of the data frame
head(datKME)


#=====================================================================================
#
#  Code chunk5.3.7-5
#
#=====================================================================================


FilterGenes= abs(GS1)> .2 & abs(datKME$MM.brown)>.8
table(FilterGenes)


#=====================================================================================
#
#  Code chunk5.3.7-6
#
#=====================================================================================


dimnames(data.frame(datExpr))[[2]][FilterGenes]


#=====================================================================================
#
#  Code chunk5.3.7-7
#
#=====================================================================================



par(mfrow=c(2,2))
# We choose 4 modules to plot: turquoise, blue, brown, green. 
# For simplicity we write the code out explicitly for each module.
which.color="turquoise" 
restrictGenes=colorh1==which.color 
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes], 
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                   col=which.color, 
                   xlab="Intramodular Connectivity", 
                   ylab="(Module Membership)^6")

which.color="blue" 
restrictGenes=colorh1==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                   col=which.color,
                   xlab="Intramodular Connectivity",
                   ylab="(Module Membership)^6")

which.color="brown" 
restrictGenes=colorh1==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                   col=which.color,
                   xlab="Intramodular Connectivity",
                   ylab="(Module Membership)^6")

which.color="green"
restrictGenes=colorh1==which.color 
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes], 
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                   col=which.color, 
                   xlab="Intramodular Connectivity", 
                   ylab="(Module Membership)^6")


#=====================================================================================
#
#  Code chunk5.3.7-8
#
#=====================================================================================


NS1=networkScreening(y=y, datME=datME, datExpr=datExpr,
                     oddPower=3, blockSize=1000, minimumSampleSize=4,
                     addMEy=TRUE, removeDiag=FALSE, weightESy=0.5)


#=====================================================================================
#
#  Code chunk5.3.7-9
#
#=====================================================================================


# network screening analysis
mean(NoiseGeneIndicator[rank(NS1$p.Weighted,ties.method="first")<=100])
# standard analysis based on the correlation p-values (or Student T test)
mean(NoiseGeneIndicator[rank(NS1$p.Standard,ties.method="first")<=100]) 


#=====================================================================================
#
#  Code chunk5.3.7-10
#
#=====================================================================================


topNumbers=c(10,20,50,100)
for (i in c(1:length(topNumbers)) ) 
{
        print(paste("Proportion of noise genes in the top", topNumbers[i], "list"))
        WGCNApropNoise=mean(NoiseGeneIndicator[rank(NS1$p.Weighted,ties.method="first")<=topNumbers[i]])
        StandardpropNoise=mean(NoiseGeneIndicator[rank(NS1$p.Standard,ties.method="first")<=topNumbers[i]])
        print(paste("WGCNA, proportion of noise=", WGCNApropNoise, 
                    ", Standard, prop. noise=", StandardpropNoise))
        if (WGCNApropNoise< StandardpropNoise) print("WGCNA wins")
        if (WGCNApropNoise==StandardpropNoise) print("both methods tie")
        if (WGCNApropNoise>StandardpropNoise) print("standard screening wins")
} 


#=====================================================================================
#
#  Code chunk5.3.7-11
#
#=====================================================================================


rm(dissTOM) collectGarbage()


#=====================================================================================
#
#  Code chunk5.3.7-12
#
#=====================================================================================


#Form a data frame containing standard and network screening results
CorPrediction1=data.frame(GS1,NS1$cor.Weighted)
cor.Weighted=NS1$cor.Weighted
# Plot the comparison

verboseScatterplot(cor.Weighted, GS1,
                   main="Network-based weighted correlation versus Pearson correlation\n",
                   col=truemodule, cex.main = 1.2)
abline(0,1)


#=====================================================================================
#
#  Code chunk5.3.7-13
#
#=====================================================================================


set.seed(2)
nSamples2=2000
MEgreen=rnorm(nSamples2)
scaledy2=MEgreen*ESgreen+sqrt(1-ESgreen^2)*rnorm(nSamples2)
y2=ifelse( scaledy2>median(scaledy2),2,1)
MEturquoise= ESturquoise*scaledy2+sqrt(1-ESturquoise^2)*rnorm(nSamples2)
# we simulate a strong dependence between MEblue and MEturquoise
MEblue= .6*MEturquoise+ sqrt(1-.6^2) *rnorm(nSamples2)
MEbrown= ESbrown*scaledy2+sqrt(1-ESbrown^2)*rnorm(nSamples2)
MEyellow= ESyellow*scaledy2+sqrt(1-ESyellow^2)*rnorm(nSamples2)
# Put together a data frame of eigengenes
ModuleEigengeneNetwork2=data.frame(y=y2,MEturquoise,MEblue,MEbrown,MEgreen, MEyellow)
# Simulate the expression data
dat2=simulateDatExpr5Modules(MEturquoise=ModuleEigengeneNetwork2$MEturquoise,
                             MEblue=ModuleEigengeneNetwork2$MEblue,MEbrown=ModuleEigengeneNetwork2$MEbrown,
                             MEyellow=ModuleEigengeneNetwork2$MEyellow,
                             MEgreen=ModuleEigengeneNetwork2$MEgreen,simulateProportions=simulateProportions1, 
                             nGenes=nGenes1)
# recall that this is the signed gene significance in the training data
GS1= as.numeric(cor(y, datExpr, use="p"))
# the following is the signed gene significance in the test data
GS2=as.numeric( cor(ModuleEigengeneNetwork2$y, dat2$datExpr, use="p"))


#=====================================================================================
#
#  Code chunk5.3.7-14
#
#=====================================================================================



par(mfrow=c(1,1))
verboseScatterplot(GS1,GS2,
                   main="Trait-based gene significance in test set vs. training set\n",
                   xlab = "Training set gene significance",
                   ylab = "Test set gene significance",
                   col=truemodule, cex.main = 1.4)


#=====================================================================================
#
#  Code chunk5.3.7-15
#
#=====================================================================================


EvaluationGeneScreening1 = corPredictionSuccess(
        corPrediction = CorPrediction1, 
        corTestSet=GS2,
        topNumber=seq(from=20, to=500, length=30) )
par(mfrow=c(2,2))
listcomp = EvaluationGeneScreening1$meancorTestSetOverall
matplot(x = listcomp$topNumber,
        y = listcomp[,-1], 
        main="Predicting positive and negative correlations",
        ylab="mean cor, test data", 
        xlab="top number of genes in the training data")
listcomp= EvaluationGeneScreening1$meancorTestSetPositive
matplot(x = listcomp$topNumber,
        y = listcomp[,-1], 
        main="Predicting positive correlations",
        ylab="mean cor, test data", 
        xlab="top number of genes in the training data")
listcomp= EvaluationGeneScreening1$meancorTestSetNegative
matplot(x = listcomp$topNumber,
        y = listcomp[,-1], 
        main = "Predicting negative correlations",
        ylab = "mean cor, test data", 
        xlab = "top number of genes in the training data")


#=====================================================================================
#
#  Code chunk5.3.7-16
#
#=====================================================================================


relativeCorPredictionSuccess(corPredictionNew = NS1$cor.Weighted,
                             corPredictionStandard = GS1, 
                             corTestSet=GS2,
                             topNumber=c(10,20,50,100,200,500) )


#=====================================================================================
#
#  Code chunk5.3.7-17
#
#=====================================================================================


# Create a data frame holding the results of gene screening
GeneResultsNetworkScreening=data.frame(GeneName=row.names(NS1), NS1)
# Write the data frame into a file
write.table(GeneResultsNetworkScreening, file="GeneResultsNetworkScreening.csv",
            row.names=F,sep=",")
# Output of eigengene information:
datMEy = data.frame(y, datME)
eigengeneSignificance = cor(datMEy, y)
eigengeneSignificance[1,1] = (1+max(eigengeneSignificance[-1, 1]))/2
eigengeneSignificance.pvalue = corPvalueStudent(eigengeneSignificance, nSamples = length(y))
namesME=names(datMEy)
# Form a summary data frame
out1=data.frame(t(data.frame(eigengeneSignificance,
                             eigengeneSignificance.pvalue, namesME, t(datMEy))))
# Set appropriate row names
dimnames(out1)[[1]][1]="EigengeneSignificance"
dimnames(out1)[[1]][2]="EigengeneSignificancePvalue"
dimnames(out1)[[1]][3]="ModuleEigengeneName"
dimnames(out1)[[1]][-c(1:3)]=dimnames(datExpr)[[1]]
# Write the data frame into a file
write.table(out1, file="MEResultsNetworkScreening.csv", row.names=TRUE, col.names = TRUE, sep=",")
# Display the first few rows:
head(out1)


#=====================================================================================
#
#  Code chunk5.3.7-18
#
#=====================================================================================


# Write out gene information
GeneName=dimnames(datExpr)[[2]]
GeneSummary=data.frame(GeneName, truemodule, SignalGeneIndicator,  NS1)
write.table(GeneSummary, file="GeneSummaryTutorial.csv", row.names=F,sep=",")
# here we output the module eigengenes and trait y without eigengene significances
datTraits=data.frame(ArrayName, datMEy)
dimnames(datTraits)[[2]][2:length(namesME)]=paste("Trait",  
                                                  dimnames(datTraits)[[2]][2:length(namesME)], 
                                                  sep=".")
write.table(datTraits, file="TraitsTutorial.csv", row.names=F,sep=",")
rm(datTraits)
# here we output the simulated gene expression data
MicroarrayData=data.frame(GeneName, t(datExpr))
names(MicroarrayData)[-1]=ArrayName
write.table(MicroarrayData, file="MicroarrayDataTutorial.csv", row.names=F,sep=",")
rm(MicroarrayData)


#=====================================================================================
#
#  Code chunk5.3.7-19
#
#=====================================================================================


# Perform network screening
NS1GS=networkScreeningGS(datExpr=datExpr, datME = datME, GS=GS1)
# Organize its results for easier plotting
GSprediction1=data.frame(GS1,NS1GS$GS.Weighted)
GS.Weighted=NS1GS$GS.Weighted
# Plot a comparison between standard gene significance and network-weighted gene significance

par(mfrow=c(1,1))
verboseScatterplot(GS1, GS.Weighted, 
                   main="Weighted gene significance vs. the standard GS\n",
                   col=truemodule)
abline(0,1)


#=====================================================================================
#
#  Code chunk5.3.7-20
#
#=====================================================================================


EvaluationGeneScreeningGS = corPredictionSuccess(corPrediction=GSprediction1, corTestSet=GS2,
                                                 topNumber=seq(from=20, to=500, length=30) )

par(mfrow=c(2,2))
listcomp= EvaluationGeneScreeningGS$meancorTestSetOverall
matplot(x=listcomp$topNumber,
        y=listcomp[,-1],
        main="Predicting positive and negative correlations",
        ylab="mean cor, test data",
        xlab="top number of genes in the training data")
listcomp= EvaluationGeneScreeningGS$meancorTestSetPositive
matplot(x=listcomp$topNumber, 
        y=listcomp[,-1], 
        main="Predicting positive correlations",
        ylab="mean cor, test data",
        xlab="top number of genes in the training data")
listcomp= EvaluationGeneScreeningGS$meancorTestSetNegative
matplot(x=listcomp$topNumber,
        y=listcomp[,-1],
        main="Predicting negative correlations",
        ylab="mean cor, test data",
        xlab="top number of genes in the training data")


#=====================================================================================
#
#  Code chunk5.3.7-21
#
#=====================================================================================


collectGarbage()
save.image("Simulated-Screening.RData")


#5.3.8 Visualization of gene networks=========


#=====================================================================================
#
#  Code chunk5.3.8-1
#
#=====================================================================================

library(cluster)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Load the previously saved data
load("Simulated-RelatingToExt.RData") 
load("Simulated-Screening.RData")


#=====================================================================================
#
#  Code chunk5.3.8-2
#
#=====================================================================================


cmd1=cmdscale(as.dist(dissTOM),2)

par(mfrow=c(1,1))
plot(cmd1, col=as.character(colorh1),  main="MDS plot",
     xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")


#=====================================================================================
#
#  Code chunk5.3.8-3
#
#=====================================================================================


power=6
color1=colorDynamicTOM
restGenes= (color1 != "grey")
diss1=1-TOMsimilarityFromExpr( datExpr[, restGenes], power = 6 )
hier1=hclust(as.dist(diss1), method="average" )
diag(diss1) = NA

TOMplot(diss1^4, hier1, as.character(color1[restGenes]),
        main = "TOM heatmap plot, module genes" )


#=====================================================================================
#
#  Code chunk5.3.8-4
#
#=====================================================================================


power=6
color1=colorDynamicTOM
restGenes= (color1 != "grey")
diss1=1-adjacency( datExpr[, restGenes], power = 6 )
hier1=hclust(as.dist(diss1), method="average" )
diag(diss1) = NA

TOMplot(diss1^4, hier1, as.character(color1[restGenes]),
        main = "Adjacency heatmap plot, module genes" )


#=====================================================================================
#
#  Code chunk5.3.8-5
#
#=====================================================================================



topList=rank(NS1$p.Weighted,ties.method="first")<=30
gene.names= names(datExpr)[topList]
# The following shows the correlations between the top genes
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
                   networkType="signed", useTOM=FALSE,
                   power=1, main="signed correlations")


#=====================================================================================
#
#  Code chunk5.3.8-6
#
#=====================================================================================



# The following shows the correlations between the top genes
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
                   networkType="unsigned", useTOM=FALSE,
                   power=1, main="signed correlations")


#=====================================================================================
#
#  Code chunk5.3.8-7
#
#=====================================================================================



# The following shows the TOM heatmap in a signed network
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
                   networkType="signed", useTOM=TRUE,
                   power=12, main="C. TOM in a signed network")
# The following shows the TOM heatmap in a unsigned network
plotNetworkHeatmap(datExpr, plotGenes = gene.names,
                   networkType="unsigned", useTOM=TRUE,
                   power=6, main="D. TOM in an unsigned network")


