#compare sets of transcriptional profiles in several ways:
#    - by overall similarity eg using Spearman correlation
#    - by similarity of signatures - usually done using GSEA (my preferred approach)
#    - by first computing modules (using WCGNA) in each dataset,
#        then comparing overlap of modules using hypergeometric test

########################################################################
#
#  0 check and install all cran and bioconductor packages if necessary
# 
# ######################################################################

list.of.cran.packages<- c("easypackages")
new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(easypackages)

source("https://bioconductor.org/biocLite.R")
list.of.bio.packages <- c("WGCNA")
#                          "AnnotationDbi", "impute")
#                          "org.Mm.eg.db","GO.db", "preprocessCore")
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)

libraries(list.of.bio.packages,list.of.cran.packages)

########################################################################################
#
# 5.1 WGCNA: Data input, cleaning and pre-processing
#  https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf
#
########################################################################################
#====5.1.1 Reading in the data (required)=============================
# Display the current working directory
getwd()
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
if (Sys.info()[['sysname']]=="Darwin"){
        setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Read in rnaseq_MicroArrayData
rnaseq_MicroArrayData = read.csv("rnaseq_MicroArrayData.csv",row.names = 1)
#colnames are changed, be careful. 
#space is replaced by "." ; "x" is add in front of number

# Take a quick look at what is in the data set:
dim(rnaseq_MicroArrayData)
names(rnaseq_MicroArrayData)
rownames(rnaseq_MicroArrayData[1:10,])

#====5.1.2  Checking data for excessive missing values (required)=============
#           and identification of outlier microarray samples ===================

# the data files contain extra information about the surveyed probes we do not need.
# We now remove the auxiliary data and transpose the expression data for further analysis.
datExpr0 = as.data.frame(t(rnaseq_MicroArrayData))
names(datExpr0[,1:10]) 
rownames(datExpr0) 

gsg = goodSamplesGenes(datExpr0, verbose = 3)

gsg$allOK
# If TRUE, all genes have passed the cuts. 


#If not
#-----remove the offending genes and samples (Alternative)------------------------------------------
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

#======cluster the samples to find outlier (required)===================

# in contrast to clustering genes that will come later) to see if there are any obvious outliers.
sampleTree = hclust(dist(datExpr0), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
# pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

#----- remove outlier (Alternative)------------------------------------
# Plot a line to show the cut
abline(h = 150, col = "red")
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
# clust 0 contains the samples we want to keep.
keepSamples = (clust==0)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# The variable datExpr now contains the expression data ready for network analysis.


#====5.1.3  Loading clinical trait data (required)===================
traitData = read.csv("ClinicalTraits.csv",row.names = 1)
dim(traitData)
head(traitData)

# remove columns that hold information we do not need.
allTraits = traitData
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.

StemCellSamples = rownames(datExpr)
traitRows = match(StemCellSamples, allTraits$Cell)
datTraits = as.data.frame(allTraits[traitRows, ])
rownames(datTraits) = allTraits[, 1]

collectGarbage()

#======visualize dendrogram(required)================================
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits[,"Protocol_number"], signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")
save(datExpr, datTraits, file = "Danwei_StemCell-01-dataInput.RData")


########################################################################################
#
# 5.2      WGCNA: Network construction and module detection
# 5.2.a    Automatic, one-step network construction and module detection
# https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-auto.pdf
#
########################################################################################
#====5.2.a.1 Reading in the data (required)=============================
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
if (Sys.info()[['sysname']]=="Darwin"){
        setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
#enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "Danwei_StemCell-01-dataInput.RData");
# The variable lnames contains the names of loaded variables.
lnames

#====5.2.a.2  Choosing the soft-thresholding power: (required)===============
# analysis of network topology

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Error in datk[c(startG:endG), ] = foreach(t = actualThreads, .combine = rbind) %dopar%  : 
# number of items to replace is not a multiple of replacement length
# Use R instead of RStudio
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.0,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



#====5.2.a.3 One-step network construction and module detection(required)===============
# Constructing the gene network and identifying modules is now a simple function call:
net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "StemCellTOM", 
                       verbose = 3)
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

#==visulization ======
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# save the module assignment
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs
geneTree = net$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree, 
     file = "Danwei_StemCell-02-networkConstruction-auto.RData")

########################################################################################
#
# 5.2.b    Step-by-step network construction and module detection
# https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-man.pdf
#
########################################################################################
#====5.2.b.1 Reading in the data (required)=============================
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
if (Sys.info()[['sysname']]=="Darwin"){
        setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);

# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
#enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "Danwei_StemCell-01-dataInput.RData");
# The variable lnames contains the names of loaded variables.
lnames

#====5.2.b.2  Choosing the soft-thresholding power: (required)===============
# analysis of network topology

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Error in datk[c(startG:endG), ] = foreach(t = actualThreads, .combine = rbind) %dopar%  : 
# number of items to replace is not a multiple of replacement length
# Use R instead of RStudio
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.0,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#====5.2.b.3 Co-expression similarity and adjacency(Required)===========
softPower = 6
adjacency = adjacency(datExpr, power = softPower)


#====5.2.b.3 Topological Overlap Matrix (TOM)(Required)===========
# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency) #takes forever
dissTOM = 1-TOM


#====5.2.b.4 Clustering using TOM(Required)===============

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# Dynamic Tree Cut for branch cutting
# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

#plot the module assignment under the gene dendrogram:

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

#====5.2.b.5 Merging of modules whose expression profiles are very similar(Required)==========
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# We choose a height cut of 0.25, corresponding to correlation of 0.75, to merge
MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

#To see what the merging did to our module colors, 
#plot the gene dendrogram again, 
#with the original and merged module colors underneath

sizeGrWindow(12, 9)
# pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# dev.off()

#We save the relevant variables for use in subsequent parts of the tutorial:
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "Danwei_StemCell-02-networkConstruction-stepByStep.RData")

#######################################################################################
#
# 5.3      WGCNA: Relating modules to external clinical traits and identifying important genes
# https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-03-relateModsToExt.pdf
#
########################################################################################

# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
if (Sys.info()[['sysname']]=="Darwin"){
        setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Load the expression and trait data saved in the first part
lnames = load(file = "Danwei_StemCell-01-dataInput.RData")
# The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Danwei_StemCell-02-networkConstruction-auto.RData");
lnames


#===5.3.a Quantifying module-trait associations==============
# In this analysis we would like to identify modules 
# that are significantly associated with the measured clinical traits.

# we simply correlate eigengenes with external
# traits and look for the most significant associations:

# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#==color code each association by the correlation value=====
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
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

#===5.3.b Gene relationship to trait and important modules: 
#  Gene Significance and Module Membership ================
# Define variable Protocol containing the Protocol column of datTrait
Protocol = as.data.frame(datTraits$Protocol_number)
names(Protocol) = "Protocol"
#names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(datExpr, Protocol, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(Protocol), sep="")
names(GSPvalue) = paste("p.GS.", names(Protocol), sep="")


#===5.3.c Intramodular analysis: identifying genes with high GS and MM
# Using the GS and MM measures, we can identify genes 
# that have a high significance for Protocol as well as high module
# membership in interesting modules.


module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body Protocol",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


#===5.3.d Summary output of network analysis results=========

# Our expression data are only annotated by probe ID names:------
# the command will return all probe IDs included in the analysis.
names(datExpr)
names(datExpr)[moduleColors=="brown"]
# To facilitate interpretation of the results, we use a probe
# annotation file provided by the manufacturer of the expression arrays
# to connect probe IDs to gene names and
# universally recognized identification numbers (Entrez codes).

#annot = read.csv(file = "GeneAnnotation.csv");
#dim(annot)
#names(annot)
#probes = names(datExpr)
#probes2annot = match(probes, annot$substanceBXH)
#   The following is the number or probes without annotation:
#sum(is.na(probes2annot))
#   Should return 0.

# We now create a data frame holding the following information for all probes
# Create the starting data frame
geneInfo0 = data.frame(
        #                       substanceBXH = probes,
        geneSymbol = colnames(datExpr),
        #                       LocusLinkID = annot$LocusLinkID[probes2annot],
        moduleColor = moduleColors,
        geneTraitSignificance,
        GSPvalue)
# Order modules by their significance for Protocol
modOrder = order(-abs(cor(MEs, Protocol, use = "p")));
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
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$Protocol));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo.csv")
fix(geneInfo)


#######################################################################################
#
# 5.4   WGCNA: Interfacing network analysis with other data such as functional
#       annotation and gene ontology 
#       https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-04-Interfacing.pdf
#
########################################################################################

# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
if (Sys.info()[['sysname']]=="Darwin"){
        setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
# Load the expression and trait data saved in the first part
lnames = load(file = "Danwei_StemCell-01-dataInput.RData")
# The variable lnames contains the names of loaded variables.
lnames
# Load network data saved in the second part.
lnames = load(file = "Danwei_StemCell-02-networkConstruction-auto.RData");
lnames

#===5.4.a Output gene lists for use with online software and services

# https://www.biostars.org/p/10457/
#  http://www.ensembl.org/biomart/martview
#  chose filters/gene/input external references ID list/ HGNC symbols

#  Read in the probe annotation
library(biomaRt)
ensembl=useMart("ENSEMBL_MART_ENSEMBL",dataset = "hsapiens_gene_ensembl")  # using ensembl database data
filter <- listFilters(ensembl)  # check which filters are available
grep("locus",filter[,1])
filter[link,1]

Attributes<-listAttributes(ensembl) # check attributes are available to select.More information on ensembl data base
Attributes[grep("link",Attributes[,1]),]
mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), mart = ensembl)
head(mapping)

#annot = read.csv(file = "GeneAnnotation.csv");
#  Match probes in the data set to the probe IDs in the annotation file 
#probes = names(datExpr)
#probes2annot = match(probes, annot$substanceBXH)
#  Get the corresponding Locuis Link IDs
#allLLIDs = annot$LocusLinkID[probes2annot];
#   $ Choose interesting modules
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