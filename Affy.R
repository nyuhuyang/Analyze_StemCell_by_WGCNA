#compare sets of transcriptional profiles in several ways:
#    - by overall similarity eg using Spearman correlation
#    - by similarity of signatures - usually done using GSEA (my preferred approach)

# =====================================================================================
# 
#  0 check and install all cran and bioconductor packages if necessary
# 
# =====================================================================================

list.of.cran.packages <- c("WGCNA","easypackages")
new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(easypackages)

source("https://bioconductor.org/biocLite.R")
list.of.bio.packages <- c("affy", "oligo", "limma")
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)
libraries(list.of.bio.packages,list.of.cran.packages)

# =====================================================================================
# 
#  1.a Loading expression data (Affymetrix)
# 
# =====================================================================================


# Display the current working directory
getwd();
MacPath="/Users/yah2014/Dropbox/Public/Olivier/R/"
PCPath="C:/Users/User/Dropbox/Public/Olivier/R"
# If necessary, change the path below to the directory where the data files are stored. 
if (Sys.info()[['sysname']]=="Darwin"){
    setwd(paste0(MacPath,"Danwei_StemCell/dataset"))}
if (Sys.info()[['sysname']]=="Windows"){
    setwd(paste0(PCPath,"Danwei_StemCell/dataset"))}
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

# Load the WGCNA package
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
#  1.a Loading expression data (Danwei)
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
#   samples (Danwei)
# 
# =====================================================================================


gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
# If the last statement returns TRUE, all genes have passed the cuts. 
# If not, we remove the offending genes and samples from the data:

# =====================================================================================
# 
#  1.b Checking data for excessive missing values and identification of outlier microarray
#   samples (Danwei)
#   remove the offending genes and samples from the data:
# 
# =====================================================================================


if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0) 
        printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0) 
        printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
dim(datExpr0)

# =====================================================================================
# 
#  1.a Loading expression data (Danwei)
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
#  1.a Loading expression data (Danwei)
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
#   samples (Danwei)
# 
# =====================================================================================


gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
# If the last statement returns TRUE, all genes have passed the cuts. 
# If not, we remove the offending genes and samples from the data:

# =====================================================================================
# 
#  1.b Checking data for excessive missing values and identification of outlier microarray
#   samples (Danwei)
#   remove the offending genes and samples from the data:
# 
# =====================================================================================


if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0) 
        printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0) 
        printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
dim(datExpr0)

