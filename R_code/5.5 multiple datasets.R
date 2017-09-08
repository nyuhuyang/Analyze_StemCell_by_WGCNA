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
#   5.5 Meta-analyses of data from two (or more) microarray data sets
#   
#   https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/JMiller/#   https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/ModulePreservation/Tutorials/
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

#Next we will determine modules based on data set A1 
m.Colorh=NULL
for (ds in 0:3){
        tree = cutreeHybrid(dendro = geneTree[[1]], pamStage=FALSE,
                            minClusterSize = (30-3*ds), cutHeight = 0.99, 
                            deepSplit = ds, distM = dissTOMA[1, , ])
        m.Colorh=cbind(m.Colorh,labels2colors(tree$labels));
}

png(file = "Plots/Module_choices.png",width = 1112, height = 807)

plotDendroAndColors(geneTree[[1]], m.Colorh, 
                    paste("dpSplt =",0:3), 
                    main = "Module choices",
                    dendroLabels=FALSE);
dev.off()
A1.modules =  m.Colorh[,1] # (Chosen based on plot below)


#  Next we calculate the principle components for visualizations

PCs1A    = moduleEigengenes(t(multiExpr.data[[set]]),
                            colors=A1.modules) 
ME_1A    = PCs1A$eigengenes
distPC1A = 1-abs(cor(ME_1A,use="p"))
distPC1A = ifelse(is.na(distPC1A), 0, distPC1A)
pcTree1A = hclust(as.dist(distPC1A),method="a") 
MDS_1A   = cmdscale(as.dist(distPC1A),2)
colorsA1 = names(table(A1.modules))

save.image("multiplemodules.RData")
#load("multiplemodules.RData")
png(file = "Plots/ModuleEigengeneVisualizations.png", 
    width = 1112, height = 807)
par(mfrow=c(2,1), mar=c(2, 4, 1, 1) + 0.1, cex=1)
plot(pcTree1A, xlab="",ylab="",main="",sub="")

plot(MDS_1A, col= colorsA1,  main="MDS plot", cex=2, pch=19)

#ordergenes = geneTree[[set]]$order
#autoplot(as.matrix(multiExpr.data[[set]]),
#         rlabels= A1.modules[ordergenes], 
#         clabels= colnames(multiExpr.data[[set]]), 
#         rcols=A1.modules[ordergenes])

#for (which.module in names(table(A1.modules))){
#        ME = ME_1A[, paste("ME",which.module, sep="")] 
#        barplot(ME, col=which.module, main="", cex.main=2, 
#                ylab="eigengene expression",xlab="array sample") 
#} 

dev.off()

##=======5.5.5 Qualitatively and quantitatively measure network preservation at the module level.
png(file = "Plots/Final_modules.png", 
    width = 1112, height = 807)

plotDendroAndColors(geneTree[[1]],
                    A1.modules,
                    "Modules", 
                    dendroLabels=FALSE,
                    hang=0.03, 
                    addGuide=TRUE, 
                    guideHang=0.05, 
                    main=paste0("Gene dendrogram and module colors " ,authors[1]))
plotDendroAndColors(geneTree[[2]], 
                    A1.modules, 
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
multiColor = list(A1 = A1.modules)

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
