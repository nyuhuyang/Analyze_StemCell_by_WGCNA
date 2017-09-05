#Contant:
#  5.1 Network analysis of expression data from XX:
#        5.1.1 Data input and cleaning 
#        5.1.2 Network construction and module detection
#                5.1.2a. Automatic, one-step network construction and module detection
#                5.1.2b Step-by-step network construction and module detection
#                5.1.2c Dealing with large datasets: block-wise network construction and module detection
#        5.1.3 Relating modules to external clinical traits and identifying important genes
#        5.1.4 Interfacing network analysis with other data such as functional annotation and gene ontology  (required)=============================
#        5.1.5 Network visualization using WGCNA functions
#        5.1.6 Export of networks to external software
#  5.2 Consensus analysis of naive state Stem Cell and Danwei's expression data

########################################################################
#
#   0 preparation and parameter adjustion
#  
#  ######################################################################

list.of.cran.packages<- c("easypackages","WGCNA")
new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(easypackages)

source("https://bioconductor.org/biocLite.R")
list.of.bio.packages <- c("biomaRt","impute")
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)

libraries(list.of.bio.packages,list.of.cran.packages)

#  Display the current working directory=======================
getwd()

if (Sys.info()[['sysname']]=="Darwin"){
        setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
if (Sys.info()[['sysname']]=="Linux"){
        setwd("/pbtech_mounts/homes030/yah2014/R/Danwei_StemCell/dataset");getwd();list.files()}
#  The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
#  Allow multi-threading within WGCNA. This helps speed up certain calculations.
#  At present this call is necessary for the code to work.
#  Any error here may be ignored but you may want to update WGCNA if you see one.
#  Caution: skip this line if you run RStudio or other third-party R environments. 
#  See note above.

if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) disableWGCNAThreads()
if (is.na(Sys.getenv("RSTUDIO", unset = NA)))  enableWGCNAThreads()

#  extract the top 5000 most variant genes for WGCNA studies.
#  set up the extract
#extract <- TRUE
extract <- FALSE
########################################################################################
#
#  5.1 Network analysis of expression data from XX:
#   finding modules related to cell_names
#
#   https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
#   labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html
########################################################################################
#====5.1.1 Data input and cleaning (required)=============================

#=====================================================================================
#
#  Code chunk5.1.1-1
#
#=====================================================================================
list_files <- read.csv("list_files.csv",header = T)
list_files
matrix <- read.csv(as.character(list_files$dataset[1]))
#prime_hESC <- read.csv("./Huang_Cell_2014/prime_hESC.csv")
#matrix <- merge(matrix,prime_hESC,by ="X")
#write.csv(matrix,"Xie_dataset.csv")
for(i in c(2:9)){ matrix1 <- read.csv(as.character(list_files$dataset[i]))
matrix <- merge(matrix,matrix1,by ="X") }
dim(matrix)
head(matrix[,1:3])
rownames(matrix) <-matrix[,1]
matrix <- matrix[,-1]
#  ==normalize==
matrix.scale <- scale(matrix)
#  check that we get mean of 0 and sd of 1
allmean <- mean(colMeans(matrix))  #  faster version of apply(scaled.dat, 2, mean)
matrix.scale <- matrix.scale +allmean 

#  colnames
sample_names <- read.csv("sample_names.csv",row.names = 1,header = T)
colnames(matrix.scale) <- sample_names$file_names
#  Read in rnaseq_MicroArrayData
Xie = matrix.scale[,c(28:45)] # 
#  Take a quick look at what is in the data set:
dim(Xie)
colnames(Xie)


#=====================================================================================
#
#  Code chunk5.1.1-2
#
#=====================================================================================


datExpr0 = as.data.frame(t(Xie))
head(names(datExpr0))
rownames(datExpr0)


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
        #  Optionally, print the gene and sample names that were removed:
        if (sum(!gsg$goodGenes)>0) 
                printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
        if (sum(!gsg$goodSamples)>0) 
                printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
        #  Remove the offending genes and samples from the data:
        datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


#=====================================================================================
#
#  Code chunk5.1.1-5
#
#=====================================================================================
#  replace all datExpr0 with Xie.Data

sampleTree = hclust(dist(datExpr0), method = "average")
#  Plot the sample tree: Open a graphic output window of size 12 by 9 inches
#  The user should change the dimensions if the window is too large or too small.
#sizeGrWindow(12,9)
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


#  Plot a line to show the cut
abline(h = 100, col = "red")
#  Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 100, minSize = 10)
table(clust)
#  clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#  extract the top 5000 most variant genes for WGCNA studies.
#  transpose matrix to correlate genes in the following
if(extract==TRUE){datExpr = datExpr[,order(apply(datExpr,2,mad), decreasing = T)[1:10000]]}


#-------------------------------------------------------------------------
#
#  Code chunk5.1.1-7
#
#-------------------------------------------------------------------------


traitData = read.csv("./Huang_Cell_2014/ClinicalTraits.csv",row.names = 1)
dim(traitData)
names(traitData)

#  remove columns that hold information we do not need.
#  convert character to numeric in r
allTraits = apply(traitData,2,function(x) as.numeric(as.factor(x)))
allTraits <- as.data.frame(allTraits)
names(allTraits)
rownames(allTraits) <- traitData$file_names
head(allTraits,10)
#  Form a data frame analogous to expression data that will hold the clinical traits.

#replace all XieSamples with Xie.Samples


Xie.Samples = rownames(datExpr)
traitRows = match(Xie.Samples, traitData$file_names)
datTraits = allTraits[traitRows, ]
#rownames(datTraits) = allTraits[traitRows, 1]

collectGarbage()

#-------------------------------------------------------------------------
#
#  Code chunk5.1.1-8
#
#-------------------------------------------------------------------------


#  Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
#  Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)
#  Error in numbers2colors(datTraits[, "cell_names"], signed = FALSE) :
#  'x' must be numeric. For a factor, please use as.numeric(x) in the call.

#  Plot the sample dendrogram and the colors underneath.
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


save(datExpr, datTraits, file = "Xie-01-dataInput.RData")

#====5.1.2 Network construction and module detection (required)=============================

#====5.1.2a. Automatic, one-step network construction and module detection
#  Skip

#====5.1.2b Step-by-step network construction and module detection======
#=====================================================================================
#
#  Code chunk5.1.2b-1
#
#=====================================================================================

#  Load the data saved in the first part
lnames = load(file = "Xie-01-dataInput.RData")
#  The variable lnames contains the names of loaded variables.
lnames
#  extract the top 5000 most variant genes for WGCNA studies.
#  transpose matrix to correlate genes in the following
if(extract==TRUE){datExpr = datExpr[,order(apply(datExpr,2,mad), decreasing = T)[1:10000]]}

#=====================================================================================
#
#  Code chunk5.1.2b-2
#  labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html
#=====================================================================================


#  Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
#  Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#  Plot the results:

par(mfrow = c(1,2))
par(mar=c(2,2,2,2))
par(oma=c(2,2,2,2))
cex1 = 0.9
#  Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
#  this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
#  Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#  chose power for Unsigned and signed hybrid networks networks 
#  https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html
#cuts <- c(0, 20, 30, 40,60, Inf)
#chose.power <- c(10, 9,8,7,6)
#power <- chose.power[findInterval(nrow(datExpr), cuts)]
#power

#  We first built adjacency matrices for each dataset using a soft power
#  threshold of 60 @ Kevin Huang Cell Stem Cell 2014
power=60
#=====================================================================================
#
#  Code chunk5.1.2b-3
#
#=====================================================================================
#  generated similarity matrices S based on Pearson correlations between all gene pairs"
S <- 0.5 + 0.5*cor(datExpr, method="spearman") #Sij = 0.5 + 0.5 × cor(i,j)


#  (1) Co-expression similarity and adjacency
#  built adjacency matrices for each dataset using a soft power threshold of 60
adjacency = adjacency.fromSimilarity(S, power = power,type = "signed") #type = "signed"! #power =softPower
SubGeneNames <-rownames(adjacency)
collectGarbage()
#=====================================================================================
#
#  Code chunk5.1.2b-4f
#  https://www.researchgate.net/post/What_do_adjacency_matrix_and_Topology_Overlap_Matrix_from_WGCNA_package_tell_about_the_data
#=====================================================================================
#  To minimize effects of noise and spurious associations
#  Turn adjacency into topological overlap
#  calculate the corresponding dissimilarity
TOM = TOMsimilarity(adjacency)
rownames(TOM) <- SubGeneNames
colnames(TOM) <- SubGeneNames #same like above
dissTOM = 1-TOM
collectGarbage()
#=====================================================================================
#
#  Code chunk5.1.2b-5
#
#=====================================================================================


#  Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")
#  Plot the resulting clustering tree (dendrogram)
par(mfrow = c(1,1))
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)


#=====================================================================================
#
#  Code chunk5.1.2b-6
#
#=====================================================================================

#  dynamicTreeCut
#  We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
#  Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)


#=====================================================================================
#
#  Code chunk5.1.2b-7
#
#=====================================================================================

#  Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
#  Plot the dendrogram and colors underneath

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


#=====================================================================================
#
#  Code chunk5.1.2b-8
#
#=====================================================================================

#  Merging of modules whose expression profiles are very similar
#  Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
#  moduleEigengenes represents the module expressions of the q-th module by the module eigengene E
#  The eigengene E can be thought of as a weighted average expression profile.

MEs = MEList$eigengenes
#  Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)
#  Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
#  Plot the result

plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")


#=====================================================================================
#
#  Code chunk5.1.2b-9
#
#=====================================================================================

#  choose a height cut of 0.25, corresponding to correlation of 0.75, to merge
MEDissThres = 0.25
#  Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
#  Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
#  The merged module colors
mergedColors = merge$colors
#  Eigengenes of the new merged modules:
mergedMEs = merge$newMEs


#=====================================================================================
#
#  Code chunk5.1.2b-10
#
#=====================================================================================


#  To see what the merging did to our module colors, 
#  we plot the gene dendrogram again, 
#  with the original and merged module colors underneath
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#dev.off()
#plotDendroAndColors(geneTree,mergedColors,"Merged dynamic",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
#=====================================================================================
#
#  Code chunk5.1.2b-11
#
#=====================================================================================


#  Rename to moduleColors
moduleColors = mergedColors
#  Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
#  Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "Xie-02-networkConstruction-stepByStep.RData")

#  Extract modules
#module_colors <- setdiff(unique(mergedColors), "grey")
#for (color in module_colors){
#        module=SubGeneNames[which(mergedColors==color)]
#        write.table(module, paste("Xie.module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
#}
#  Look at expression patterns of these genes, as they are clustered

module.order <- unlist(tapply(1:ncol(datExpr),as.factor(mergedColors),I)) #dynamicColors
m<-t(t(datExpr[,module.order])/apply(datExpr[,module.order],2,max))
par(mfrow = c(1,1))
par(oma = c(10,2,1,1))
cex1 = 0.9
heatmap(t(m),zlim=c(0,1),col=gray.colors(100),Rowv=NA,Colv=NA,labRow=NA,scale="none",RowSideColors=mergedColors[module.order])
#===> jump to chunk5.1.6-2

#===5.1.2c Dealing with large datasets: block-wise network construction and module detection====
#  Skip

#######################################################################################
#====5.1.3 Relating modules to external clinical traits and identifying important genes (required)=============================
#=====================================================================================
#
#  Code chunk5.1.3-1
#
#=====================================================================================

#  The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
#  Load the expression and trait data saved in the first part
lnames = load(file = "Xie-01-dataInput.RData")
#  The variable lnames contains the names of loaded variables.
lnames
#  Load network data saved in the second part.
lnames = load(file = "Xie-02-networkConstruction-stepByStep.RData")
lnames


#=====================================================================================
#
#  Code chunk5.1.3-2
#
#=====================================================================================


#  Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#  Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p") #method = "spearman"
#  Remove columns from dataframe where some of values are NA
moduleTraitCor <- moduleTraitCor[ , apply(moduleTraitCor, 2, function(x) !any(is.na(x)))]

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


#=====================================================================================
#
#  Code chunk5.1.3-3
#
#=====================================================================================



#  Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mfrow =c(1,1))
par(mar = c(6, 8.5, 3, 3))
#  Display the correlation values within a heatmap plot
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
#  WGCNA::greenWhiteRed: this palette is not suitable for people
#  with green-red color blindness (the most common kind of color blindness).
#  Consider using the function blueWhiteRed instead.

#=====================================================================================
#
#  Code chunk5.1.3-4
#
#=====================================================================================


#  Define variable cell_names containing the cell_names column of datTrait
cell_names = as.data.frame(datTraits$cell_names)
names(cell_names) = "cell_names"
#  names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

geneTraitSignificance = as.data.frame(cor(datExpr, cell_names, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(cell_names), sep="")
names(GSPvalue) = paste("p.GS.", names(cell_names), sep="")

#=====================================================================================
#
#  Code chunk5.1.3-5
#
#=====================================================================================


module = "bisque4" #  based on Module-trait relationships
column = match(module, modNames)
moduleGenes = moduleColors==module


par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for cell stage",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module) #yellow

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
names(datExpr)[moduleColors==module]

#---------------------------------------------------------------#
#
#  Code chunk5.1.3-8
#
#---------------------------------------------------------------#
#  Skip this 

#annot = read.csv(file = "GeneAnnotation.csv")
#dim(annot)
#names(annot)
#probes = names(datExpr)
#probes2annot = match(probes, annot$substanceBXH)
#  The following is the number or probes without annotation:
#sum(is.na(probes2annot))
#  Should return 0.


#=====================================================================================
#
#  Code chunk5.1.3-9
#
#=====================================================================================
#  We now create a data frame holding the following information for all probes
#  Create the starting data frame
geneInfo0 = data.frame(
        #      substanceBXH = probes,
        geneSymbol = colnames(datExpr),
        #      LocusLinkID = annot$LocusLinkID[probes2annot],
        moduleColor = moduleColors,
        geneTraitSignificance,
        GSPvalue)
#  Order modules by their significance for cell_names
modOrder = order(-abs(cor(MEs, cell_names, use = "p")))
#  Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
        oldNames = names(geneInfo0)
        geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                               MMPvalue[, modOrder[mod]])
        names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                             paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
head(geneInfo0[1:3,])
#  Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.cell_names))
geneInfo = geneInfo0[geneOrder, ]


#=====================================================================================
#
#  Code chunk5.1.3-10
#
#=====================================================================================
write.csv(geneInfo, file = "Xie.geneInfo.csv")

#====5.1.4 Interfacing network analysis with other data such as functional annotation and gene ontology  (required)=============================
#  Skip

#===5.1.5 Network visualization using WGCNA functions===========


#=====================================================================================
#
#  Code chunk5.1.5-1
#
#=====================================================================================
#  The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
#  Load the expression and trait data saved in the first part
lnames = load(file = "Xie-01-dataInput.RData")
#  The variable lnames contains the names of loaded variables.
lnames
#  Load network data saved in the second part.
#lnames = load(file = "Xie-02-networkConstruction-auto.RData")
#  or----------------
lnames = load(file = "Xie-02-networkConstruction-stepByStep.RData")

lnames
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


#=====================================================================================
#
#  Code chunk5.1.5-2
#
#=====================================================================================

#  Calculate topological overlap anew: this could be done more efficiently by saving the TOM
#  calculated during module detection, but let us do it again here.
#dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 7)
#  Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
#plotTOM = dissTOM^7
#  Set diagonal to NA for a nicer plot
#diag(plotTOM) = NA

#  Call the plot function

#TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes") #takes long time
#  WARNING: On some computers, this code can take a while to run (20 minutes??). 
#  I suggest you skip it.

#=====================================================================================
#
#  Code chunk5.1.5-3
#
#=====================================================================================

#  select 400 genes
#nSelect = 400
#  For reproducibility, we set the random seed
#set.seed(10)
#select = sample(nGenes, size = nSelect)
#selectTOM = dissTOM[select, select]
#  There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
#selectTree = hclust(as.dist(selectTOM), method = "average")
#selectColors = moduleColors[select]


#  Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
#  the color palette setting the diagonal to NA also improves the clarity of the plot
#plotDiss = selectTOM^7
#diag(plotDiss) = NA
#TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")


#=====================================================================================
#
#  Code chunk5.1.5-4
#
#=====================================================================================


#  Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
#  Isolate weight from the clinical traits
cell_names = as.data.frame(datTraits$cell_names)
names(cell_names) = "cell_names"
#  Add the cell_names to existing module eigengenes
MEs = orderMEs(cbind(MEs, cell_names))
MET = orderMEs(MEs)
#  Plot the relationships among the eigengenes and the trait

par(cex = 0.9)
par(oma = c(2,2,2,2))
par(oma = c(2,2,2,2))
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", 
                      marDendro = c(0,4,1,2), 
                      marHeatmap = c(3,4,1,2), 
                      cex.lab = .8, 
                      xLabelsAngle = 90)


#=====================================================================================
#
#  Code chunk5.1.5-5
#
#=====================================================================================


#  To split the dendrogram and heatmap plots, 
#  we can use the following code

par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram",
                      marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
#  Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap",
                      marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE,
                      xLabelsAngle = 90)

#===5.1.6 Export of networks to external software=============
# Skip






















