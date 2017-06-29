#compare sets of transcriptional profiles in several ways:
#    - by overall similarity eg using Spearman correlation
#    - by similarity of signatures - usually done using GSEA (my preferred approach)

# =====================================================================================
# 
#  0 check and install all cran and bioconductor packages if necessary
# 
# =====================================================================================

list.of.cran.packages<- c("WGCNA","easypackages","stringr")
new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(easypackages)

source("https://bioconductor.org/biocLite.R")
list.of.bio.packages <- c("AnnotationDbi", "impute", "org.Mm.eg.db","GO.db",
                          "edgeR","preprocessCore")
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)

libraries(list.of.bio.packages,list.of.cran.packages)
########################################################################################
#
#  A. Data input and cleaning
#  Yang Hu Jun 29 2016
#
########################################################################################

# =====================================================================================
# 
#  1.a Expression data---stem cell stages
# 
# =====================================================================================

# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
if (Sys.info()[['sysname']]=="Darwin"){
    setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
    setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}

#GSE29397_series_matrix.txt come from Vassena R 2011 paper,
#Recorded human stem cell from Oocyte to Blastocyst stages.

#To import data into R with an unknown number of columns:
no_col <- max(count.fields("GSE29397_series_matrix.txt", sep = "\t")) 
#count.fields counts number of column per row

MicroArrayData0 = read.table("GSE29397_series_matrix.txt",sep = "\t",
                            fill=TRUE,col.names=1:no_col) #fill empty

#remove header rows. Header rows have "!".
header<-grepl("!",MicroArrayData0[,1]) #grepl returns a logical vector. 
MicroArrayData<- MicroArrayData0[!header,] #remove header rows


#rename rows and columns
#for columns:
col_names <- MicroArrayData0[MicroArrayData0[,1]=="!Sample_title",]
col_names <- droplevels(col_names)
col_names<- as.factor(unlist(col_names))
colnames(MicroArrayData) <- col_names

#for rows:
rownames(MicroArrayData) <-MicroArrayData[,1]
#remove row 1
MicroArrayData<- MicroArrayData[-1,] 

# -----------------------------------------------------------------
#  Loading annotation data (MINiML)---stem cell stages

anno0 = read.table("GPL6244-tbl-1.txt",sep = "\t",fill=TRUE) #fill empty
anno<- str_split_fixed(anno0[,10], " // ",3) #split column 10
anno<-cbind(anno0[,1],anno[,2]) #extract ID_REF and Gene name only
anno <- anno[!anno[,2]=="",] #remove empty columns with no gene names
colnames(anno) <-c("!Sample_title","gene")

anno_MicroArrayData <- merge(anno, MicroArrayData)
anno_MicroArrayData <- anno_MicroArrayData[order(anno_MicroArrayData[,"gene"]),]
anno_MicroArrayData <- anno_MicroArrayData[,-1]
collectGarbage();
# =====================================================================================
# 
#  1.b Expression data---rnaseq_rowcounts
#  consult to http://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html
# =====================================================================================
#Reading in the data

countdata <- read.csv("rnaseq_rowcounts.csv")
colnames(countdata)[1] <- "gene"
rownames(countdata) <- countdata[,"gene"]
countdata<-countdata[,-1]
# -----------------------------------------------------------------
#Filtering to remove lowly expressed genes
# Obtain CPMs
myCPM <- cpm(countdata)
# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
table(rowSums(thresh))
# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
# Subset the rows of countdata to keep the more highly expressed genes
myCPM.keep <- myCPM[keep,]


# (DEG visulization, optional)-------------------------------------------------------------
# Library sizes and distribution plots
# The las argument rotates the axis names
expr0 <- DGEList(countdata)
expr1 <- DGEList(myCPM.keep)

par(mfrow=c(2,1))

barplot(expr0$samples$lib.size,names=colnames(expr),las=2)
title("Barplot of library sizes")
barplot(expr1$samples$lib.size,names=colnames(expr),las=2)

# Get log2 counts per million
logcounts0 <- cpm(expr0,log=TRUE)
logcounts1 <- cpm(expr1,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts0, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts0),col="blue")
title("Boxplots of logCPMs (unnormalised)")
boxplot(logcounts1, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts1),col="blue")
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

#Multidimensional scaling plots
plotMDS(expr0)
plotMDS(expr1)

#put gene column back
myCPM.keep <- cbind(rownames(myCPM.keep),myCPM.keep) 
colnames(myCPM.keep)[1] <- "gene"
collectGarbage();
# =====================================================================================
# 
#  1.c overall similarity
# 
# =====================================================================================

rnaseq_MicroArrayData <- merge(myCPM.keep, anno_MicroArrayData)














x# Load the WGCNA package
library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Read in the female liver data set
rnaseqData = read.csv("rnaseq_rowcounts.csv");
# Take a quick look at what is in the data set:
dim(rnaseqData);
names(rnaseqData);


# =====================================================================================
# 
#  1.a Loading expression data
#  remove the auxiliary data and transpose
# 
# =====================================================================================
# the data files contain extra information about the surveyed probes we do not need.
# We now remove the auxiliary data and transpose the expression data for further analysis.

datExpr0 = as.data.frame(t(rnaseqData[, -1]));
names(datExpr0) = rnaseqData$X;
rownames(datExpr0) = names(rnaseqData)[-1];


# =====================================================================================
# 
#  1.b Checking data for excessive missing values and identification of outlier microarray
#   samples
# 
# =====================================================================================


gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
# If the last statement returns TRUE, all genes have passed the cuts. 
# If not, we remove the offending genes and samples from the data:

# =====================================================================================
# 
#  1.b Checking data for excessive missing values and identification of outlier microarray samples
#   remove the offending genes and samples from the data:
# 
# =====================================================================================


if (!gsg$allOK)
{
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0) 
        printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0) 
        printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


# =====================================================================================
#
#  1.b Checking data for excessive missing values and identification of outlier microarray samples
#  1.b (Continue) cluster the samples to find outlier
# 
# =====================================================================================

# in contrast to clustering genes that will come later) to see if there are any obvious outliers.
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
# pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


# =====================================================================================
# 
#  1.b Checking data for excessive missing values and identification of outlier microarray samples
#  1.b (Continue) remove outlier
# 
# =====================================================================================


# Plot a line to show the cut
abline(h = 200000, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
# clust 0 contains the samples we want to keep.
keepSamples = (clust==0)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# The variable datExpr now contains the expression data ready for network analysis.

# =====================================================================================
# 
#  1.c Loading clinical trait data
# 
# =====================================================================================


traitData = read.csv("ClinicalTraits.csv",row.names = 1);
dim(traitData)
names(traitData)

# remove columns that hold information we do not need.
allTraits = traitData[, -2];
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.

femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$Cell);
datTraits = as.data.frame(allTraits[traitRows, -1]);
rownames(datTraits) = allTraits[traitRows, 1];


collectGarbage();


# =====================================================================================
#
#  1.c Loading clinical trait data
#  1.c (Continue) visualize how the clinical traits relate to the sample dendrogram
# 
# =====================================================================================


# Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "complete")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")


# =====================================================================================
#
#  1.c Loading clinical trait data
#  1.c (Continue) save the relevant expression and trait data for use in the next steps of the tutorial.
# 
# =====================================================================================


save(datExpr, datTraits, file = "Danwei_StemCell-01-dataInput.RData")

##############################################################################
# 
# 
#  2.a Automatic network construction and module detection
# 
# 
##############################################################################



# =====================================================================================
# 
#  0. Preliminaries: setting up the R session
# 
# =====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
workingDir = ".";
setwd(workingDir); list.files()
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "Danwei_StemCell-01-dataInput.RData");
# The variable lnames contains the names of loaded variables.
lnames


# =====================================================================================
# 
#  2.a.1 Choosing the soft-thresholding power: analysis of network topology
# 
# =====================================================================================


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Error in datk[c(startG:endG), ] = foreach(t = actualThreads, .combine = rbind) %dopar%  : 
# number of items to replace is not a multiple of replacement length
# Use R instead of RStudio
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


# =====================================================================================
# 
#  2.a.2 One-step network construction and module detection
# 
# =====================================================================================
# Constructing the gene network and identifying modules is now a simple function call:


``
#Error: REAL() can only be applied to a 'numeric', not a 'integer'
# A word of caution.
# default value for blockwiseModules may be not approprate for other data.

# A second word of caution concerning block size.
# maxBlockSize =5000
# Note that if this code were to be used to analyze a data set with more than 5000 probes,
# the function blockwiseModules will split the data set into several blocks.
# This will break some of the plotting code below

# Analyze larger data sets need to do one of the following:
# 16GB should handle up to 20000 probes;

# If a computer with large-enough memory is not available, 
# the reader should follow Section 2.c, Dealing with large datasets


# To see how many modules were identified and what the module sizes are

table(net$colors)
# and indicates that there are 18 modules,labeled 1 through 18 in order of descending size,
# with sizes ranging from 609 to 34 genes. 
# The label 0 is reserved for genes outside of all modules.

# =====================================================================================
# 
#  2.a.2 One-step network construction and module detection
#  2.a.2 (Continue) The dendrogram
# 
# =====================================================================================


# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


# =====================================================================================
# 
#  2.a.2 One-step network construction and module detection
#  2.a.2 (Continue) save the module assignment and module eigengene information necessary
#  for subsequent analysis.
# 
# =====================================================================================


moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "Danwei_StemCell-02-networkConstruction-auto.RData")

##############################################################################
# 
# 
#  2.b Step-by-step network construction and module detection
# 
# 
##############################################################################


# =====================================================================================
# 
#  0 Preliminaries: setting up the R session
# 
# =====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
if (Sys.info()[['sysname']]=="Darwin"){
    setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
    setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell");getwd();list.files()}
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "Danwei_StemCell-01-dataInput.RData");
# The variable lnames contains the names of loaded variables.
lnames


# =====================================================================================
# 
#  2.b.1 Choosing the soft-thresholding power: analysis of network topology
# 
# =====================================================================================


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


# =====================================================================================
# 
#  2.b.2 Co-expression similarity and adjacency
# 
# =====================================================================================


softPower = 6;
adjacency = adjacency(datExpr, power = softPower);


# =====================================================================================
# 
#  2.b.3 Topological Overlap Matrix (TOM)
# 
# =====================================================================================

# To minimize effects of noise and spurious associations, 
# we transform the adjacency into Topological Overlap Matrix,
# and calculate the corresponding dissimilarity:
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM


# =====================================================================================
# 
#  2.b.4 Clustering using TOM
# 
# =====================================================================================


# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);


# =====================================================================================
# 
#  2.b.4 Clustering using TOM
#  2.b.4 (Continue) Dynamic Tree Cut for branch cutting
# 
# =====================================================================================


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


# =====================================================================================
# 
#  2.b.4 Clustering using TOM
#  2.b.4 (Continue) plot the module assignment under the gene dendrogram:
# 
# =====================================================================================


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


# =====================================================================================
# 
#  2.b.5 Merging of modules whose expression profiles are very similar
# 
# =====================================================================================


# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge

# =====================================================================================
# 
#  2.b.5 Merging of modules whose expression profiles are very similar
#  2.b.5 (Continue) merging of modules
# 
# =====================================================================================


MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


# =====================================================================================
# 
#  2.b.5 Merging of modules whose expression profiles are very similar
#  2.b.5 (Continue) dendrogram with the original and merged module colors underneath
# 
# =====================================================================================


sizeGrWindow(12, 9)
# pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# dev.off()


# =====================================================================================
# 
#  2.b.5 Merging of modules whose expression profiles are very similar
#  2.b.5 (Continue) We save the relevant variables for use in subsequent parts of the tutorial:
# 
# =====================================================================================


# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "Danwei_StemCell-02-networkConstruction-stepByStep.RData")

##############################################################################
# 
# 
#  3. Relating modules to external information and identifying important
#  
# 
# 
##############################################################################


# =====================================================================================
# 
#  0 Preliminaries: setting up the R session and loading results of previous parts
# 
# =====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
if (Sys.info()[['sysname']]=="Darwin"){
    setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
    setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell");getwd();list.files()}
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "Danwei_StemCell-01-dataInput.RData");
# The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Danwei_StemCell-02-networkConstruction-auto.RData");
lnames


# =====================================================================================
# 
#  3.a Quantifying module-trait associations
# 
# =====================================================================================
# In this analysis we would like to identify modules 
# that are significantly associated with the measured clinical traits.

# we simply correlate eigengenes with external
# traits and look for the most significant associations:

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


# =====================================================================================
# 
#  3.a Quantifying module-trait associations
#  (Continue) color code each association by the correlation value
# 
# =====================================================================================


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


# =====================================================================================
# 
#  3.b Gene relationship to trait and important modules: 
#  Gene Significance and Module Membership
# 
# =====================================================================================


# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"
names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");


# =====================================================================================
# 
#  3.c Intramodular analysis: identifying genes with high GS and MM
# 
# =====================================================================================
# Using the GS and MM measures, we can identify genes 
# that have a high significance for weight as well as high module
# membership in interesting modules.


module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


# =====================================================================================
# 
#  3.d Summary output of network analysis results
# 
# =====================================================================================

# Our expression data are only annotated by probe ID names:
# the command will return all probe IDs included in the analysis.
names(datExpr)


# =====================================================================================
# 
#  3.d Summary output of network analysis results
#  3.d (Continue) Return probe IDs belonging to the brown module.
# 
# =====================================================================================


names(datExpr)[moduleColors=="brown"]


# =====================================================================================
# 
#  3.d Summary output of network analysis results
#  3.d (Continue) interpretation of the results
# 
# =====================================================================================
# To facilitate interpretation of the results, we use a probe
# annotation file provided by the manufacturer of the expression arrays
# to connect probe IDs to gene names and
# universally recognized identification numbers (Entrez codes).

annot = read.csv(file = "GeneAnnotation.csv");
dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.


# =====================================================================================
# 
#  3.d Summary output of network analysis results
#  3.d (Continue) create a data frame holding informations
# 
# =====================================================================================

# We now create a data frame holding the following information for all probes


# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = probes,
                       geneSymbol = annot$gene_symbol[probes2annot],
                       LocusLinkID = annot$LocusLinkID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                           MMPvalue[, modOrder[mod]]);
    names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                         paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
geneInfo = geneInfo0[geneOrder, ]


# =====================================================================================
# 
#  3.d Summary output of network analysis results
#  3.d (Continue) This data frame can be written into a text-format spreadsheet,
# 
# =====================================================================================


write.csv(geneInfo, file = "geneInfo.csv")
fix(geneInfo)
##############################################################################
# 
# 
#  4. Interfacing network analysis with other data such as functional
#  annotation and gene ontology  
# 
# 
##############################################################################


# =====================================================================================
# 
#  0 Preliminaries: setting up the R session and loading results of previous parts
# 
# =====================================================================================


#=====================================================================================
#
#  0 Preliminaries: setting up the R session
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
if (Sys.info()[['sysname']]=="Darwin"){
    setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
    setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell");getwd();list.files()}
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "Danwei_StemCell-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Danwei_StemCell-02-networkConstruction-auto.RData");
lnames


#=====================================================================================
#
#  4.a Output gene lists for use with online software and services
#
#=====================================================================================


# Read in the probe annotation
annot = read.csv(file = "GeneAnnotation.csv");
# Match probes in the data set to the probe IDs in the annotation file 
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# Get the corresponding Locuis Link IDs
allLLIDs = annot$LocusLinkID[probes2annot];
# $ Choose interesting modules
intModules = c("brown", "red", "salmon")
for (module in intModules)
{
    # Select module probes
    modGenes = (moduleColors==module)
    # Get their entrez ID codes
    modLLIDs = allLLIDs[modGenes];
    # Write them into a file
    fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
    write.table(as.data.frame(modLLIDs), file = fileName,
                row.names = FALSE, col.names = FALSE)
}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("LocusLinkIDs-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)


#=====================================================================================
#
#  4.b Enrichment analysis directly within R
#
#=====================================================================================
#The function takes a vector of module labels,
#and the Entrez (a.k.a. Locus Link) codes for the genes
#whose labels are given.

GOenr = GOenrichmentAnalysis(moduleColors, allLLIDs, organism = "mouse", nBestP = 10);

#The function runs for awhile and returns a long list,
#the most interesting component of which is

tab = GOenr$bestPTerms[[4]]$enrichment
#This is an enrichment table containing
#the 10 best terms for each module present in moduleColors.

#Names of thecolumns within the table can be accessed by
names(tab)

#it is best to save the table into a file 
#and open it using their favorite tool
write.table(tab, file = "GOEnrichmentTable.csv", sep = ",", quote = TRUE, row.names = FALSE)


#=====================================================================================
#
#  4.b Enrichment analysis directly within R
#  display above result directly on screen
#
#=====================================================================================


keepCols = c(1, 2, 5, 6, 7, 12, 13);
screenTab = tab[, keepCols];
# Round the numeric columns to 2 decimal places:
numCols = c(3, 4);
screenTab[, numCols] = signif(apply(screenTab[, numCols], 2, as.numeric), 2)
# Truncate the the term name to at most 40 characters
screenTab[, 7] = substring(screenTab[, 7], 1, 40)
# Shorten the column names:
colnames(screenTab) = c("module", "size", "p-val", "Bonf", "nInTerm", "ont", "term name");
rownames(screenTab) = NULL;
# Set the width of R's output. The reader should play with this number to obtain satisfactory output.
options(width=95)
# Finally, display the enrichment table:
screenTab

##############################################################################
# 
# 
#  5. Network visualization using WGCNA functions
# 
# 
##############################################################################

#=====================================================================================
#
#  0 Preliminaries: setting up the R session and loading results of previous parts
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
if (Sys.info()[['sysname']]=="Darwin"){
    setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
    setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell");getwd();list.files()}

# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "Danwei_StemCell-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Danwei_StemCell-02-networkConstruction-auto.RData");
lnames
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


#=====================================================================================
#
#  5 Visualization of networks within R
#  5.a Visualizing the gene network by heatmap
#
#=====================================================================================

#This code can be executed only if the network
#was calculated using a single-block approach

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6); # it takes some time
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
#this may stuck

#=====================================================================================
#
#  5.a Visualizing the gene network by heatmap
#  with less gene=400
#
#=====================================================================================


nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")


#=====================================================================================
#
#  5.b Visualizing the network of eigengenes
#
#=====================================================================================


# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# Isolate weight from the clinical traits
weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, weight))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)


#=====================================================================================
#
#  5.b Visualizing the network of eigengenes
#  To split the dendrogram and heatmap plots
#
#=====================================================================================


# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)

##############################################################################
# 
# 
#  6. Exporting a gene network to external visualization software
# 
# 
##############################################################################

#=====================================================================================
#
#  0 Preliminaries: setting up the R session and loading results of previous parts
#
#=====================================================================================


# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
if (Sys.info()[['sysname']]=="Darwin"){
    setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
    setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell");getwd();list.files()}

# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
lnames = load(file = "Danwei_StemCell-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Danwei_StemCell-02-networkConstruction-auto.RData");
lnames


#=====================================================================================
#
#  6 Exporting network data to network visualization software
#  6.a Exporting to VisANT
#
#=====================================================================================


# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power = 6);#It takes some time.
# Read in the annotation file
annot = read.csv(file = "GeneAnnotation.csv");
# Select module
module = "brown";
# Select module probes
probes = names(datExpr)
inModule = (moduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )


#=====================================================================================
#
#  6.a Exporting to VisANT
#  restrict the genes in the output to say the 30 top hub genes in the module:
#
#=====================================================================================


nTop = 30;
IMConn = softConnectivity(datExpr[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("VisANTInput-", module, "-top30.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(annot$substanceBXH, annot$gene_symbol) )
#To provide an example of a VisANT visualization, 
#we loaded the file produced by the above code in VisANT.

#=====================================================================================
#
#  6.b Exporting to Cytoscape
#
#=====================================================================================


# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = 6); #It takes some time.
# Read in the annotation file
annot = read.csv(file = "GeneAnnotation.csv");
# Select modules
modules = c("brown", "red");
# Select module probes
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])



