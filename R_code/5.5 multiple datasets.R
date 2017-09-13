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
#  5.3 Analysis of simulated data
#  5.4 module preservation
#  5.5 Meta-analyses of data from two (or more) microarray data sets
########################################################################
#
#   0 preparation and parameter adjustion
#  
#######################################################################

list.of.cran.packages<- c("easypackages","flashClust","Hmisc",
                          "ggfortify")
new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(easypackages)

source("https://bioconductor.org/biocLite.R")
list.of.bio.packages <- c("biomaRt","impute","WGCNA")
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)

libraries(list.of.bio.packages,list.of.cran.packages)

#  Display the current working directory=======================
getwd()

if (Sys.info()[['sysname']]=="Darwin"){
        setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset/WGCNA/metaAnalysisFiles");getwd();list.files()}
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
#   5.5 Meta-analyses of data from two (or more) microarray data sets
#   
#       https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/JMiller/
########################################################################################


#=======5.5.1 load data ================
list_files <- read.csv("list_files.csv",header = T)
list_files
matrix <- read.csv(as.character(list_files$dataset[1]))

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
#  rownames
gene.names <- rownames(matrix)
#  colnames
sample.names <- read.csv("sample_names.csv",row.names = 1,header = T)
colnames(matrix.scale) <- sample.names$file_names

#  authors group
authors <- unique(sample.names$authors)
authors <- authors[!(authors %in% c("Gu Y et al",
                                    "Theunissen et al",
                                    "Tachibana M et al",
                                    "Chan et al",
                                    "Hanna et al",
                                    ""))]

#=======5.5.2 clean data ================
#  We work with 6 sets:
nSets = length(authors)
#  Form multi-set expression data: columns starting from 9 contain actual expression data.
multiExpr = vector(mode = "list", length = nSets)

for (set in 1:nSets){
        exp = matrix.scale[,sample.names$stages == "naive" &
                                   sample.names$authors ==authors[set]]
        multiExpr[set] = list(data = as.data.frame((exp)))
}


#=======5.5.3 Correlating general network properties================
softPower = 60 # (Read WGCNA tutorial to learn how to pick your power)
multiExpr.rank = vector(mode = "list", length = nSets)
for (set in 1:nSets){
        multiExpr.rank[[set]] = rank(rowMeans(multiExpr[[set]]))
}
random5000= sample(gene.names,5000)
multiConn.rank = vector(mode = "list", length = nSets)
for (set in 1:nSets){
        multiConn.rank[[set]] = rank(softConnectivity(t(multiExpr[[set]][random5000,]),
                                                      type="signed",
                                                      power=softPower))
}

png(file = "Plots/generalNetworkProperties.png", 
    width = 1112, height = 807)
par(mfrow=c(3,2))
verboseScatterplot(multiExpr.rank[[1]],
                   multiExpr.rank[[2]],
                   xlab=paste0("Ranked Expression", authors[1]),
                   ylab=paste0("Ranked Expression", authors[2]))
verboseScatterplot(multiConn.rank[[1]],
                   multiConn.rank[[2]],
                   xlab=paste0("Ranked Connectivity", authors[1]),
                   ylab=paste0("Ranked Connectivity", authors[2]))                  
verboseScatterplot(multiExpr.rank[[3]],
                   multiExpr.rank[[4]],
                   xlab=paste0("Ranked Expression", authors[3]),
                   ylab=paste0("Ranked Expression", authors[4]))
verboseScatterplot(multiConn.rank[[3]],
                   multiConn.rank[[4]],
                   xlab=paste0("Ranked Connectivity", authors[3]),
                   ylab=paste0("Ranked Connectivity", authors[4]))                  
verboseScatterplot(multiExpr.rank[[5]],
                   multiExpr.rank[[6]],
                   xlab=paste0("Ranked Expression", authors[5]),
                   ylab=paste0("Ranked Expression", authors[6]))
verboseScatterplot(multiConn.rank[[5]],
                   multiConn.rank[[6]],
                   xlab=paste0("Ranked Connectivity", authors[5]),
                   ylab=paste0("Ranked Connectivity", authors[6]))                  
dev.off()

#  =======5.5.4 Run WGCNA on the data sets================
#  So computational reasons and for simplicity 
#  we first will choose the top 5000 most expressed probes in data set
multiExpr.data = vector(mode = "list", length = nSets)
for (set in 1:nSets){
        keepGenesExpr= rank(-rowMeans(multiExpr[[set]]))<=5000
        multiExpr.data[set] <- list(data = multiExpr[[set]][keepGenesExpr,])
}


#  we will calculate all of the necessary values to run WGCNA.
#  Initialize an appropriate array to hold the adjacencies
adjacencies = array(0, dim = c(nSets, 5000, 5000))
dissTOMA = array(0, dim = c(nSets, 5000, 5000))
geneTree = vector("list",nSets) 
#  Calculate adjacencies, TOM and tree in each individual data set
system.time( {
for (set in 1:nSets){
        adjacencies[set, , ] <- adjacency(t(multiExpr.data[[set]]),
                                        power=softPower,
                                        type="signed")
        diag(adjacencies[set, , ])=0
        
        dissTOMA[set, , ] <- 1-TOMsimilarity(adjacencies[set, , ],
                                             TOMType="signed")
        geneTree[[set]] <- flashClust(as.dist(dissTOMA[set, , ]),
                                         method="average")
}
})

#save.image("adj_TOM_tree.RData")  #  (Section will take ~5-15 minutes to run)
#  NOTE: This file will be ~700GB and is not required, but if your computer crashes, you can type:
#' load("tutorial.RData") ' to restart at this point. 

png(file = "Plots/Gene clustering on TOM-based dissimilarity1.png", 
    width = 1112, height = 807)
par(mfrow=c(3,1))
for (set in 1:3){
        plot(geneTree[[set]],xlab="",sub="",
             main=paste0("Gene clustering on TOM-based dissimilarity ",authors[set]),
             labels=FALSE,
             hang=0.04)
        }
dev.off()

png(file = "Plots/Gene clustering on TOM-based dissimilarity2.png", 
    width = 1112, height = 807)
par(mfrow=c(3,1))
for (set in 4:6){
        plot(geneTree[[set]],xlab="",sub="",
             main=paste0("Gene clustering on TOM-based dissimilarity ",authors[set]),
             labels=FALSE,
             hang=0.04)
}
dev.off()

# set A1 
m.Colorh=NULL
ds=1
tree =NULL
for (set in 1:6){
        tree = cutreeHybrid(dendro = geneTree[[set]], pamStage=FALSE,
                            minClusterSize = (30), cutHeight = 0.99, 
                            deepSplit = ds, distM = dissTOMA[set, , ])
        m.Colorh=cbind(m.Colorh,labels2colors(tree$labels));
}

#png(file = "Plots/Module_choices.png",width = 1112, height = 807)

plotDendroAndColors(geneTree[[1]], 
                    m.Colorh, 
                    c("module colors",authors[2:6]), 
                    main = "Module choices",
                    dendroLabels=FALSE)
dev.off()
modulesA1 =  m.Colorh[,1] # (Chosen based on plot below)


#  Next we calculate the principle components for visualizations

PCs1A    = moduleEigengenes(t(multiExpr.data[[1]]),
                            colors=modulesA1) 
ME_1A    = PCs1A$eigengenes
distPC1A = 1-abs(cor(ME_1A,use="p"))
distPC1A = ifelse(is.na(distPC1A), 0, distPC1A)
pcTree1A = hclust(as.dist(distPC1A),method="a") 
MDS_1A   = cmdscale(as.dist(distPC1A),2)
colorsA1 = names(table(modulesA1))

#save.image("multiplemodules.RData")
#load("multiplemodules.RData")
png(file = "Plots/ModuleEigengeneVisualizations.png", 
    width = 1112, height = 807)
par(mfrow=c(2,1), mar=c(2, 4, 1, 1) + 0.1, cex=1)
plot(pcTree1A, xlab="",ylab="",main="",sub="")

plot(MDS_1A, col= colorsA1,  main="MDS plot", cex=2, pch=19)

#ordergenes = geneTree[[set]]$order
#autoplot(as.matrix(multiExpr.data[[set]]),
#         rlabels= modulesA1[ordergenes], 
#         clabels= colnames(multiExpr.data[[set]]), 
#         rcols=modulesA1[ordergenes])

#for (which.module in names(table(modulesA1))){
#        ME = ME_1A[, paste("ME",which.module, sep="")] 
#        barplot(ME, col=which.module, main="", cex.main=2, 
#                ylab="eigengene expression",xlab="array sample") 
#} 

dev.off()

##=======5.5.5 Qualitatively and quantitatively measure network preservation at the module level.
png(file = "Plots/Final_modules.png", 
    width = 1112, height = 807)

plotDendroAndColors(geneTree[[1]],
                    modulesA1,
                    "Modules", 
                    dendroLabels=FALSE,
                    hang=0.03, 
                    addGuide=TRUE, 
                    guideHang=0.05, 
                    main=paste0("Gene dendrogram and module colors " ,authors[1]))
plotDendroAndColors(geneTree[[2]], 
                    modulesA1, 
                    "Modules",
                    dendroLabels=FALSE,
                    hang=0.03, 
                    addGuide=TRUE, 
                    guideHang=0.05, 
                    main=paste0("Gene dendrogram and module colors " ,authors[2]))
dev.off()

#modulePreservation
multiExpr.module =  list(A1=list(data=t(multiExpr.data[[1]])),
                        A2=list(data=t(multiExpr.data[[2]])),
                        A3=list(data=t(multiExpr.data[[3]])),                   
                        A4=list(data=t(multiExpr.data[[4]])),                   
                        A5=list(data=t(multiExpr.data[[5]])),                   
                        A6=list(data6=t(multiExpr.data[[6]]))) 
multiColor = list(A1 = modulesA1)

system.time( {
        mp = modulePreservation(multiExpr.module, multiColor,
                                networkType = "signed", 
                                referenceNetworks = 1,
                                nPermutations = 30,#200, 30sec/permutations
                                maxGoldModuleSize=100,
                                randomSeed = 1,
                                quickCor = 0,
                                verbose = 3)
} ) # this takes long time

stats = mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2
Zsummary <- stats[order(-stats[,2]),c(1:2)]
write.csv(Zsummary,"Zsummary.csv")

##=======5.5.6 Module membership (kME) and its use in comparing networks
##=======5.5.6.1 get the kME values,
#along with their associated p-values for A1
#and will then output the resulting table to a file ("kMEtable1.csv")
geneModuleMembership1 = signedKME(t(multiExpr.data[[1]]), ME_1A)
colnames(geneModuleMembership1)=paste("PC",colorsA1,".cor",sep=""); 

MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(multiExpr.data[[1]])[[2]]); 
colnames(MMPvalue1)=paste("PC",colorsA1,".pval",sep="");

Gene       = rownames(multiExpr.data[[1]])
kMEtable1  = cbind(Gene,Gene,modulesA1)
for (i in 1:length(colorsA1))
        kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i])
colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1), colnames(MMPvalue1))))

write.csv(kMEtable1,"kMEtable1.csv",row.names=FALSE)


#Now repeat for A2, using the module assignments from A1 to determine kME values.

# First calculate MEs for A2, since we haven't done that yet
PCs2A = moduleEigengenes(t(multiExpr.data[[2]]),
                         colors=modulesA1) 
ME_2A = PCs2A$eigengenes

geneModuleMembership2 = signedKME(t(multiExpr.data[[2]]), ME_2A)
colnames(geneModuleMembership1)=paste("PC",colorsA1,".cor",sep=""); 

MMPvalue2=corPvalueStudent(as.matrix(geneModuleMembership2),dim(multiExpr.data[[2]])[[2]]); 
colnames(MMPvalue2)=paste("PC",colorsA1,".pval",sep="");

kMEtable2  = cbind(Gene,Gene,modulesA1)
for (i in 1:length(colorsA1))
        kMEtable2 = cbind(kMEtable2, geneModuleMembership2[,i], MMPvalue2[,i])
colnames(kMEtable2)=colnames(kMEtable1)

write.csv(kMEtable2,"kMEtable2.csv",row.names=FALSE)

##=======5.5.6.2 compare the resulting networks. 
##=======5.5.6.2a plot the kME values
pdf("Plots/all_kMEtable2_vs_kMEtable1.pdf",height=8,width=8)
for (c in 1:length(colorsA1)){
        verboseScatterplot(geneModuleMembership2[,c],geneModuleMembership1[,c],main=colorsA1[c],
                           xlab="kME in A2",ylab="kME in A1")
}; dev.off()

pdf("Plots/inModule_kMEtable2_vs_kMEtable1.pdf",height=8,width=8)
for (c in 1:length(colorsA1)){
        inMod = modulesA1== colorsA1[c]
        verboseScatterplot(geneModuleMembership2[inMod,c],geneModuleMembership1[inMod,c],main=colorsA1[c],
                           xlab="kME in A2",ylab="kME in A1")
}; dev.off()
#save.image("tutorial.RData") #(optional line of code)

#=========5.5.6.2b determine which genes are hubs in both networks 
topGenesKME = NULL
for (c in 1:length(colorsA1)){
        kMErank1    = rank(-geneModuleMembership1[,c])
        kMErank2    = rank(-geneModuleMembership2[,c])
        maxKMErank  = rank(apply(cbind(kMErank1,kMErank2+.00001),1,max))
        topGenesKME = cbind(topGenesKME,Gene[maxKMErank<=10])
}; colnames(topGenesKME) = colorsA1
topGenesKME

##=======5.5.7 Comparing networks and annotating modules using programs outside of R
source("tutorialFunctions.R")  
for (co in colorsA1[colorsA1!="grey"])
        visantPrepOverall(modulesA1, co, t(multiExpr.data[[1]]), 
                          rownames(multiExpr.data[[1]]), 500, softPower, TRUE)


#=======5.5.7.1 output data from our network for import into VisANT
for (co in colorsA1[colorsA1!="grey"])
        visantPrepOverall(modulesA1, co, t(multiExpr.data[[1]]),
                          rownames(multiExpr.data[[1]]), 500, softPower, TRUE)

#=======5.5.7.2 Hub genes specific to one network
multiExpr.data1_2 = t(cbind(multiExpr.data[[1]],multiExpr.data[[2]]))
i1 = 1:dim(multiExpr.data[[1]])[[2]];
i2 = (1:dim(multiExpr.data[[2]])[[2]])+length(i1)
for (co in colorsA1[colorsA1!="grey"])
        visantPrep(modulesA1, co, i1, i2, multiExpr.data1_2, rownames(multiExpr.data[[1]]), 500, softPower, TRUE)

#=======5.5.7.3 EASE to annotate modules based on Gene Ontology
folder = "geneLists_A12/"
for (c in colorsA1){ 
        fn = paste(folder, c, ".txt",sep=""); 
        write.geneList(Gene[modulesA1==c], fn) 
};  
write(Gene,paste(folder,"all.txt",sep=""))

#=======5.5.7.4 annotate modules based enrichment for user-defined lists
enrichments = userListEnrichment(Gene, modulesA1, c("exampleListInput.csv","exampleMMInput.csv"), 
                                 c("cellType","humanModules"), "enrichment.csv")


##======5.5.8 Using phenotypic information to determine differentially expressed genes and modules
#

enrichments$ovGenes$'green -- cellType_Astrocyte_Cahoy_all'
# The output is...
# [1] "ACSBG1"        "ADD3"          "ADHFE1"        "AGL"           "AGT"         
# [2] "AGTRL1"        "AGXT2L1"       "AHCYL1"        "AK3"           "ALDH1A1"
# (etc.)

====5.5.8.1 phenotypic information
region = rep("CA1",32);
region[c(1,4,6,11,12,15,16,17,22,24,25,26,28,29,31,32)] = "CA3"
age = c(81,72,86,90,88,90,90,74,83,73,73,70,85,85,75,90,72,70,90,84,75,85,80,
        86,85,84,81,88,80,90,83,74)

#====5.5.8.2 find the genes that show the highest differential expression
# Find the top 3 region genes
var    = list(region=="CA1", region=="CA3")
datReg = t(apply(multiExpr.data[[1]],1,t.test.l))
colnames(datReg)=c("MeanCA1","MeanCA3","SD_CA1","SD_CA3","PvalRegion")
datReg[order(datReg[,5])[1:3],]
#           MeanCA1    MeanCA3    SD_CA1    SD_CA3   PvalRegion
# NRIP3  10509.8842 22481.7087 2951.5061 6239.7370 6.739892e-07
# SAT1    4656.6850  2964.6736  885.3605  658.5044 1.335158e-06
# KCNH3   3814.7772  1840.0098 1049.7946  749.2281 1.493591e-06

# Find the top 3 age genes
var    = age
datAge = t(apply(multiExpr.data[[1]],1,cor.test.l))
colnames(datAge)=c("CorrAge","PvalAge")
datAge[order(datAge[,2])[1:3],]
#             CorrAge      PvalAge
# DDX42     0.7116889 4.947339e-06
# SCD       0.6961764 9.660017e-06
# KIAA0907 -0.6818583 1.728760e-05


#====5.5.8.3  find all of the modules that are either region or age related
# Find the region-related modules
var     = list(region=="CA1", region=="CA3")
datRegM = t(apply(t(ME_1A),1,t.test.l))
colnames(datRegM)=c("MeanCA1","MeanCA3","SD_CA1","SD_CA3","PvalRegion")
datRegM[datRegM[,5]<0.02,]
#             MeanCA1     MeanCA3    SD_CA1    SD_CA3  PvalRegion
# MEgreen  0.07875933 -0.07875933 0.1869713 0.1359238 0.011049496
# MEpink  -0.08035560  0.08035560 0.1241285 0.1936073 0.009706674

# Find the age-related modules
var     = age
datAgeM = t(apply(t(ME_1A),1,cor.test.l))
colnames(datAgeM)=c("CorrAge","PvalAge")
datAgeM[datAgeM[,2]<0.02,]
#              CorrAge     PvalAge
# MEblack    0.4485706 0.010023898
# MEmagenta -0.4502347 0.009717422
# MEred     -0.4303339 0.013952662


#====5.5.8.4 We can now visualize these results
pdf("RegionAgePlots.pdf",width=16,height=4)
par(mfrow=c(1,4))
verboseBoxplot(as.numeric(multiExpr.data[[1]]["NRIP3",]), region,
               main="NRIP3 expression -", las=2, xlab="Region", ylab="")
verboseScatterplot(age, as.numeric(multiExpr.data[[1]]["DDX42",]), 
                   main="DDX42 expression -", las=2, abline=TRUE, xlab="Age", ylab="")
verboseBoxplot(as.numeric(ME_1A[,"MEgreen"]), region,
               main="Green ME expr. -", las=2, xlab="Region", ylab="")
verboseScatterplot(age, as.numeric(ME_1A[,"MEmagenta"]), 
                   main="Magenta ME expr. -", las=2, abline=TRUE, xlab="Age", ylab="")
dev.off()


#===5.5.8.5  visually compare how a gene or module relates to phenotype across data sets
region2 = rep("CA1",31);
region2[c(1,4,7,8,9,11,15,18,20,21,22,26,28,29,30)] = "CA3"
var     = list(region2=="CA1", region2=="CA3")
datReg2 = t(apply(multiExpr.data[[2]],1,t.test.l))
colnames(datReg2)=c("MeanCA1","MeanCA3","SD_CA1","SD_CA3","PvalRegion")
datRegM2 = t(apply(t(ME_2A),1,t.test.l))
colnames(datRegM2)=c("MeanCA1","MeanCA3","SD_CA1","SD_CA3","PvalRegion")

pdf("RegionAgePlots12.pdf",width=16,height=4)
par(mfrow=c(1,4))
verboseBoxplot(as.numeric(multiExpr.data[[1]]["NRIP3",]), region,
               main="NRIP3 expression (A1) -", las=2, xlab="Region (A1)", ylab="")
verboseBoxplot(as.numeric(multiExpr.data[[2]]["NRIP3",]), region2,
               main="NRIP3 expression (A2) -", las=2, xlab="Region (A2)", ylab="")
verboseBoxplot(as.numeric(ME_1A[,"MEgreen"]), region,
               main="Green ME expr. (A1) -", las=2, xlab="Region (A1)", ylab="")
verboseBoxplot(as.numeric(ME_2A[,"MEgreen"]), region2,
               main="Green ME expr. (A2) -", las=2, xlab="Region (A2)", ylab="")
dev.off()


##======5.5.9 Comparing networks with different module definitions
#======5.5.9.1  make a network using data set B1
# (This section will take 5-15 minutes to run)

multiExpr.data[[3]] = (collapseRows(multiExpr.data[[3]],genesA,probesA))[[1]]
GeneAB = sort(intersect(rownames(multiExpr.data[[3]]g),Gene)) # There are 3396 genes in this network

multiExpr.data[[3]] = multiExpr.data[[3]]g[GeneAB,]
adjacencyB2 = adjacency(t(multiExpr.data[[3]]g),power=softPower,type="signed");
diag(adjacencyB2)=0
dissTOMB2   = 1-TOMsimilarity(adjacencyB2, TOMType="signed")
geneTreeB2  = flashClust(as.dist(dissTOMB2), method="average")

mColorh=NULL
for (ds in 0:3){
        tree = cutreeHybrid(dendro = geneTreeB2, pamStage=FALSE, minClusterSize = (30-3*ds), 
                            cutHeight = 0.99, deepSplit = ds, distM = dissTOMB2)
        mColorh=cbind(mColorh,labels2colors(tree$labels));
}
pdf("Module_choices_B2.pdf", height=10,width=25); 
plotDendroAndColors(geneTreeB2, mColorh, paste("dpSplt =",0:3), main = "",dendroLabels=FALSE);
dev.off()
modulesB2 =  mColorh[,1] # (Deepsplit=0 Chosen based on plot below)
colorsB2 = names(table(modulesB2))

#=====5.5.9.2 We can now check module overlap
modulesB2_new = matchModules(Gene, modulesA1, GeneAB, modulesB2)
# Old - New labels: black - purple   # Your screen output will look like this
# Old - New labels: blue - brown
# Old - New labels: brown - greenyellow  # etc.
enrichmentsB2A1 = userListEnrichment(GeneAB,modulesB2,"kMEtable1.csv","A1","enrichmentB2_A1.csv")
enrichmentsB2A1$sigOverlaps  # To show the significant overlaps on the screen
#      InputCategories UserDefinedCategories CorrectedPvalues      
# pOut "blue"          "A1_brown"            "2.41636641189253e-59"
# pOut "turquoise"     "A1_yellow"           "1.51940862185272e-43"
# pOut "green"         "A1_brown"            "3.67833377321179e-41"  # etc.

pdf("A1_B2_dendrogram_plots.pdf",height=5,width=15)
plotDendroAndColors(geneTreeA1, modulesA1, "A1", dendroLabels=FALSE, hang=0.03, 
                    addGuide=TRUE, guideHang=0.05, main="A1 dendrogram")
plotDendroAndColors(geneTreeB2, modulesB2_new, "B2_new", dendroLabels=FALSE, hang=0.03, 
                    addGuide=TRUE, guideHang=0.05, main="B2 dendrogram")  
dev.off()

#=====5.5.9.3 Comparing module annotations


