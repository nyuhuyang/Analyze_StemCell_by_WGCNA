# =====================================================================================
# 
#  0 check and install all cran and bioconductor packages if necessary
# 
# =====================================================================================

list.of.cran.packages<- c("WGCNA","easypackages","stringr","dplyr")
new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(easypackages)

source("https://bioconductor.org/biocLite.R")
list.of.bio.packages <- c("AnnotationDbi", "impute","GEOquery","oligo",
                          "edgeR","preprocessCore")
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)


libraries(list.of.bio.packages,list.of.cran.packages)
memory.limit(size=25000)
# =====================================================================================
# 
#  1.a Load Expression data---stem cell stages
# 
# =====================================================================================
# ------load processed data----------------------------------------------------------
eList <- getGEO("GSE29397") #it takes some time
eList
names(eList)
eData <- eList[[1]] #extract expressionset
eData
names(pData(eData))#there is usually a lot of unnecessary stuff here
pData(eData)$title #enough cell information from title
sampleNames(eData) <-pData(eData)$title
boxplot(eData)
#normData <- rma(rawData) #  already normalized, doesn't work 
#normData
# -----------------------------------------------------------------



# ------Loading raw data (optional)-----------------------------------------------------------
#too slow, enough data from above processed data
eList2 <- getGEOSuppFiles("GSE29397")
eList2
list.files("GSE29397")
untar("GSE29397/GSE29397_RAW.tar", exdir = "GSE29397/CEL")
list.files("GSE29397/CEL")
celfiles <- list.files("GSE29397/CEL", full = TRUE)
rawData <- read.celfiles(celfiles)
rawData
exprs(rawData)
max(exprs(rawData))
# clean up the phenotype information
filename <- sampleNames(rawData)
pData(rawData)$filename <- filename
sampleNames <- sub(".*_", "", filename) #remove all contain before "_"
sampleNames <- sub(".CEL.gz$", "", sampleNames)
sampleNames(rawData) <- sampleNames
pData(rawData)$group <- ifelse(grepl("^OSA", sampleNames(rawData)),
                               "OSA", "Control")
pData(rawData)
#Normalization
boxplot(rawData) # stuck and didn't work
normData <- rma(rawData) # it works
normData
boxplot(normData)
par(mfrow=c(1,2))
boxplot(eData)
boxplot(normData)
# -----------------------------------------------------------------


#  ----------Loading annotation data (MINiML)---stem cell stages----------

# Display the current working directory
getwd();
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


#  ----------Loading annotation data (GEOquery, optional)---stem cell stages----------
annotation(eData)
platf <- getGEO(annotation(eData), AnnotGPL=TRUE)
anot <- data.frame(attr(dataTable(platf), "table"))
anot1<- str_split_fixed(anot[,"Gene.symbol"], "///",2) #split column 10
anno<-cbind.data.frame(anot[,"ID"],anot1[,1])
anno <- anno[!anno[,2]=="",] #remove empty columns with no gene names
colnames(anno) <-c("PROBES","gene_ID")
# -----------------------------------------------------------------



# ------replace affymetrix ID to gene ID-----------------------------------------------------------
MicroArrayData <- exprs(eData)
PROBES <- as.numeric(rownames(MicroArrayData))
MicroArrayData<- cbind.data.frame(PROBES,MicroArrayData)

anno_MicroArrayData <- merge(anno, MicroArrayData, all = TRUE)
anno_MicroArrayData <- anno_MicroArrayData[!is.na(anno_MicroArrayData[,3]),] #remove NA row in exprs
anno_MicroArrayData <- anno_MicroArrayData[!is.na(anno_MicroArrayData[,"gene_ID"]),] #remove NA row in genes

dim(anno_MicroArrayData)
anno_MicroArrayData <- anno_MicroArrayData[order(anno_MicroArrayData[,"gene_ID"]),]
anno_MicroArrayData <- distinct(anno_MicroArrayData,gene_ID,.keep_all = TRUE) #remove duplicated rows
dim(anno_MicroArrayData)

rownames(anno_MicroArrayData) <- anno_MicroArrayData[,"gene_ID"]
anno_MicroArrayData <- anno_MicroArrayData[,-1] #remove "PROBES"


# =====================================================================================
# 
#  1.b Expression data---rnaseq_rowcounts
#  consult to http://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html
# =====================================================================================
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
if (Sys.info()[['sysname']]=="Darwin"){
        setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}

#Reading in the data

countdata <- read.csv("rnaseq_rowcounts.csv")
colnames(countdata)[1] <- "gene_ID"
rownames(countdata) <- countdata[,"gene_ID"]
countdata<-countdata[,-1]
# -----------------------------------------------------------------
#Filtering to remove lowly expressed genes
# Obtain CPMs
myCPM <- cpm(countdata)
# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 1
table(rowSums(thresh))
# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 9
# Subset the rows of countdata to keep the more highly expressed genes
myCPM.keep <- myCPM[keep,]


# (DEG visulization, optional)-------------------------------------------------------------
# Library sizes and distribution plots
# The las argument rotates the axis names
expr0 <- DGEList(countdata)
expr1 <- DGEList(myCPM.keep)

par(mfrow=c(1,2))

# Get log2 counts per million
#counts0 <- cpm(expr0,log=FALSE)
logcounts1 <- cpm(expr0,log=TRUE)
logcounts2 <- cpm(expr1,log=TRUE)
# Check distributions of samples using boxplots
#boxplot(counts0, xlab="", ylab="counts per million",las=2) # test with all((counts0)>=0)
#abline(h=median(counts0),col="blue")
#title("CPMs (unnormalised)")

boxplot(logcounts1+5.2, xlab="", ylab="Log2 counts per million",las=2) # test with all((logcounts1+5.2)>0)
abline(h=median(logcounts1+5.2),col="blue")
title("logCPMs (unnormalised)")

boxplot(logcounts2+2, xlab="", ylab="Log2 counts per million",las=2) # test with all((logcounts2+2)>0)
abline(h=median(logcounts2+2),col="blue")
title("logCPMs (normalised)")


#Multidimensional scaling plots
plotMDS(expr0)
plotMDS(expr1)

#put gene column back
logcounts <- cbind.data.frame(rownames(logcounts2),logcounts2+4) 
colnames(logcounts)[1] <- "gene_ID"

#merge 1.a and 1.b
rnaseq_MicroArrayData <- merge(logcounts, anno_MicroArrayData)
dim(rnaseq_MicroArrayData)
rownames(rnaseq_MicroArrayData) <- rnaseq_MicroArrayData[,"gene_ID"]
rnaseq_MicroArrayData <- rnaseq_MicroArrayData[,-1] #remove "gene_ID"
par(mfrow=c(1,1))
par(oma=c(12,2,0,2))
boxplot(rnaseq_MicroArrayData, xlab="", ylab="counts per million",
        main = "boxplot for log Counts Per Million of all samples",las=2) # test with all((counts0)>=0)
abline(h=median(as.matrix(rnaseq_MicroArrayData)),col="blue")
write.csv(rnaseq_MicroArrayData,"rnaseq_MicroArrayData.csv")
# =====================================================================================
# 
#  2. overall similarity
# 
# =====================================================================================
par(mfrow=c(1,1))
par(oma=c(2,2,2,2))
c <- cor(rnaseq_MicroArrayData, method="spearman")
d <- as.dist(1-c)
hc <- hclust(d, method="complete")
p<-plot(hc,col = "black",cex = 1,xlab="", sub="",
        main="Cluster Dendrogram, method='spearman'")

library(pheatmap)
pheatmap(c,cex=1.05,
         main ="Spearman correlation between all samples")