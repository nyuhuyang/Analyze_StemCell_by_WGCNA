
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

list.of.cran.packages<- c("easypackages","stringr",
                          "dplyr","gridExtra")
new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(easypackages)

source("https://bioconductor.org/biocLite.R")
list.of.bio.packages <- c("oligo","AnnotationDbi", "impute","GEOquery")
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)


libraries(list.of.bio.packages,list.of.cran.packages)


# #####################################################################
# 
#  1.1 Load microarray data---stem cell stages
#  https://kasperdanielhansen.github.io/genbioconductor/html/GEOquery.html (coursera)
#  https://f1000research.com/articles/5-1384/v1
# ####################################################################
#=========Reading in the data (required)===================
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
if (Sys.info()[['sysname']]=="Darwin"){
        setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
# ===========1.1.1 load processed data (Recommend)=======================
eList <- getGEO("GSE18290") #Vassena R 2011 GSE29397, Xie D 2010 GSE18290
eList
names(eList)
eData <- eList[[3]] #GSE29397 chose [[1]], GSE18290 chose [[3]]
eData
names(pData(eData))#there is usually a lot of unnecessary stuff here
pData(eData)$title #enough cell information from title
sampleNames(eData) <-pData(eData)$title
#saveRDS(eData,"Vassena_eData")
par(mfrow=c(1,1))
par(mar=c(16,4,3,1))
boxplot(eData, ylab = "log2 counts",las = 3,cex=1.3,
        main ="Microarray for early human pluripotent stem cells (GSE29397)")
#normData <- rma(rawData) #  already normalized, doesn't work 
#normlize


# ----------1.1.2 Loading raw data (Alternative)--------------------------------
#too slow, enough data from GSE29397
eList2 <- getGEOSuppFiles("GSE29397")
eList2
list.files("GSE29397")
untar("GSE29397/GSE29397_RAW.tar", exdir = "GSE29397/CEL")
list.files("GSE29397/CEL") 
celfiles <- list.files("GSE29397/CEL", full = TRUE)
rawData <- read.celfiles(celfiles)
rawData
head(exprs(rawData))
head(pData(rawData))
# clean up the phenotype information
filename <- sampleNames(rawData)
pData(rawData)$filename <- filename
sampleNames <- sub(".*_", "", filename) #remove all contain before "_"
sampleNames <- sub(".CEL.gz$", "", sampleNames)
sampleNames(rawData) <- sampleNames
pData(rawData)$group <- ifelse(grepl("^OSA", sampleNames(rawData)),
                               "OSA", "Control")
head(pData(rawData))
saveRDS(rawData,"rawData")
rawData <-readRDS("rawData")
#-(Alternative)--------------------------
#boxplot(rawData) # stuck and didn't work
normData <- rma(rawData) # it works
normData
head(exprs(normData))

par(mfrow=c(1,2))
boxplot(eData)
boxplot(normData)

#The MAplot also allows summarization, so groups can be compared more easily:
MAplot(rawData[, 1:4], pairs=TRUE)

# =======1.1.3 Loading raw data (Alternative)=====================
#must be done from raw data Xie D 2010 GSE18290
eList2 <- getGEOSuppFiles("GSE18290")
eList2
list.files("GSE18290")
untar("GSE18290/GSE18290_RAW.tar", exdir = "GSE18290/CEL")
list.files("GSE18290/CEL")
celfiles <- list.files("GSE18290/CEL", full = TRUE)

#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18290
human_embryo <-paste0("GSE18290/CEL/GSM",456643:456660,".CEL.gz")
celfiles %in% human_embryo
rawData <- read.celfiles(human_embryo)
rawData
head(exprs(rawData))
head(pData(rawData))
# clean up the phenotype information
pData(rawData)$filename <- sampleNames(rawData)
sampleNames(rawData) <- pData(eData)$title #Xie D 2010 GSE18290
head(pData(rawData))

boxplot(rawData) # before normalization
normData <- rma(rawData) # it works
normData
head(exprs(normData))

par(mfrow=c(1,2))
boxplot(rawData)
boxplot(normData)

#The MAplot also allows summarization, so groups can be compared more easily:
#MAplot(rawData[, 1:4], pairs=TRUE)

saveRDS(normData,"Xie_normData")

# #####################################################################
# 
#  1.2 prepare microarray data---stem cell stages
#  https://github.com/benilton/oligoOld/wiki/Getting-the-grips-with-the-oligo-Package
#
# ####################################################################

#======1.2.1 Loading annotation data (Recommend)=================
annotation(rawData)
platf <- getGEO(annotation(eData), AnnotGPL=TRUE)
anot <- data.frame(attr(dataTable(platf), "table"))
anot1<- str_split_fixed(anot[,"Gene.symbol"], "///",2) #split column Gene.symbol
anno<-cbind.data.frame(anot[,"ID"],anot1[,1])
anno <- anno[!anno[,2]=="",] #remove empty columns with no gene names
colnames(anno) <-c("Probe_ID","Gene_name")
#===========================================================================#

#---------1.2.2 Loading annotation data from biomaRT (Alternative)--------------
# https://www.biostars.org/p/10457/
#  http://www.ensembl.org/biomart/martview
#  chose filters/gene/input external references ID list/ HGNC symbols

allprobe = featureNames(eData)
library(biomaRt)
# downlaod annotationfrom biomAT
mart <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart)
annot <-getBM(attributes = c("ensembl_gene_id", #Gene stable ID
                     "external_gene_name", #Gene name
                     "go_id",   #GO term accession
                     "affy_hugene_1_0_st_v1"), #AFFY HuGene 1 0 st v1
      filters = "affy_hugene_1_0_st_v1", values = allprobe,  #save time with filter
      mart = ensembl)
table(annot$affy_hugene_1_0_st_v1 %in% allprobe)
table(table(annot$affy_hugene_1_0_st_v1))

gene2annot = match(allprobe, annot$affy_hugene_1_0_st_v1)
gene2annot <- gene2annot[!is.na(gene2annot)] #Remove all NA values from a vector
# Get the corresponding GO IDs
anno = annot[gene2annot,c("affy_hugene_1_0_st_v1","external_gene_name")] 
anyNA(anno$affy_hugene_1_0_st_v1) # test empty elements in the vector
#allLLIDs <- allLLIDs[allLLIDs !=""] #Remove the empty element in the vector
colnames(anno) <-c("Probe_ID","Gene_name")

#===============1.2.6 replace affymetrix ID to gene ID (Required)==================##
#Probe_ID <- as.numeric(rownames(exprs(eData))) 
Probe_ID <- rownames(exprs(eData)) 
MicroArrayData<- cbind.data.frame(Probe_ID,exprs(eData))

anno_MicroArrayData <- merge(anno, MicroArrayData, all = TRUE)
anno_MicroArrayData <- anno_MicroArrayData[!is.na(anno_MicroArrayData[,3]),] #remove NA row in exprs
anno_MicroArrayData <- anno_MicroArrayData[!is.na(anno_MicroArrayData[,"Gene_name"]),] #remove NA row in genes

dim(anno_MicroArrayData)
anno_MicroArrayData <- anno_MicroArrayData[order(anno_MicroArrayData[,"Gene_name"]),]
anno_MicroArrayData <- distinct(anno_MicroArrayData,Gene_name,.keep_all = TRUE) #remove duplicated rows
dim(anno_MicroArrayData)

rownames(anno_MicroArrayData) <- anno_MicroArrayData[,"Gene_name"]
anno_MicroArrayData <- anno_MicroArrayData[,-1] #remove "PROBES"
#write.csv(anno_MicroArrayData[,-1],"Vassena_dataset.csv") #remove "gene name"
write.csv(anno_MicroArrayData[,-1],"Xie_dataset.csv") #remove "gene name"

# ===========1.2.7 get gene list that have higher std (top 25%) (essential)===========
# we can find the standard deviation for each chip
# http://jura.wi.mit.edu/bio/education/bioinfo2007/arrays/array_exercises_1R.html#ratios

sd = apply(anno_MicroArrayData[,-1], 1, sd) #run sd w/o gene_ID
sd = sd[order(sd, decreasing=TRUE)]
Top_Std_gene <-names(sd)[1:(length(sd)*0.25)]


#save files
write.csv(anno,"anno")
write.csv(anno_MicroArrayData,"anno_MicroArrayData") 
write.csv(Top_Std_gene,"Top_Std_gene")  #===> next page
# ===========================================================================#


# ------------1.2.8 visulization sd ~counts (Alternative)-------------------------------------------------------------
median  = apply(anno_MicroArrayData[,-1], 1, median)
par(mfrow=c(1,1))
plot(median,sd,xlab="log2(count size)", ylab="standard deviation")
lines(lowess(median,sd), col="blue")# lowess line (x,y)

# ------1.1.3 Loading raw CEL data (locally, alternative)-----------------------------------------------------------
# Display the current working directory
getwd();
MacPath="/Users/yah2014/Dropbox/Public/Olivier/R/"
PCPath="C:/Users/User/Dropbox/Public/Olivier/R"
# If necessary, change the path below to the directory where the data files are stored. 
if (Sys.info()[['sysname']]=="Darwin"){
        setwd(paste0(MacPath,"/Danwei_StemCell/dataset"))}
if (Sys.info()[['sysname']]=="Windows"){
        setwd(paste0(PCPath,"/Danwei_StemCell/dataset"))}
getwd();list.files()

#following instruction
#http://jura.wi.mit.edu/bio/education/bioinfo2007/arrays/array_exercises_1R.html
#Read the CEL files (first command below) and then summarize and normalize with MAS5.
#This could take a few minutes.
affy.data <- ReadAffy()
eset.mas5 = mas5(affy.data)

#getting the expression matrix (probesets/genes in rows, chips in columns)
exprSet.nologs = exprs(eset.mas5)
# List the column (chip) names
colnames(exprSet.nologs)
# Rename the column names if we want
colnames(exprSet.nologs) = c("1-cell.1", "1-cell.2", "1-cell.3",
                             "2-cell.1", "2-cell.2", "2-cell.3",
                             "4-cell.1", "4-cell.2", "4-cell.3",
                             "8-cell.1", "8-cell.2", "8-cell.3",
                             "morula.1", "morula.2", "morula.3",
                             "blastocyst.1", "blastocyst.2","blastocyst.3")

#we'll do all logs with base 2
exprSet = log(exprSet.nologs, 2)

#To print out our expression matrix (as with most data), we can use a command like
write.table(exprSet, file="Su_mas5_matrix.txt", quote=F, sep="\t")

#-----------------------------------------------------------------


#  ----------1.2.3 Loading annotation data (locally,alternative)---stem cell stages----------

# Display the current working directory
getwd()
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
if (Sys.info()[['sysname']]=="Darwin"){
        setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}


anno0 = read.table("GPL6244-tbl-1.txt",sep = "\t",fill=TRUE) #fill empty
anno<- str_split_fixed(anno0[,10], " // ",3) #split column 10
anno<-cbind.data.frame(anno0[,1],anno[,2]) #extract ID_REF and Gene name only
anno <- anno[!anno[,2]=="",] #remove empty columns with no gene names
colnames(anno) <-c("PROBES","gene_ID")
# -----------------------------------------------------------------




