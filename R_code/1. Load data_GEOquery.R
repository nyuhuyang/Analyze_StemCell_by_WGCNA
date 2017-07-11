
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

list.of.cran.packages<- c("WGCNA","easypackages","stringr",
                          "dplyr","pheatmap","gridExtra")
new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(easypackages)

source("https://bioconductor.org/biocLite.R")
list.of.bio.packages <- c("AnnotationDbi", "impute","GEOquery","oligo",
                          "edgeR","DESeq2","preprocessCore")
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)


libraries(list.of.bio.packages,list.of.cran.packages)
memory.limit(size=25000)

# #####################################################################
# 
#  1.a Load microarray data---stem cell stages
#  https://kasperdanielhansen.github.io/genbioconductor/html/GEOquery.html (coursera)
#
# ####################################################################

# ===========load processed data (essential)============================================
eList <- getGEO("GSE29397") #it takes some time
eList
names(eList)
eData <- eList[[1]] #extract expressionset
eData
names(pData(eData))#there is usually a lot of unnecessary stuff here
pData(eData)$title #enough cell information from title
sampleNames(eData) <-pData(eData)$title
par(mfrow=c(1,1))
boxplot(eData)
#normData <- rma(rawData) #  already normalized, doesn't work 
#normData
#===========================================================================#


# ------Loading raw data (alternative)-----------------------------------------------------------
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

# ------Loading raw CEL data (locally, alternative)-----------------------------------------------------------
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
# -----------------------------------------------------------------


# #####################################################################
# 
#  1.b prepare microarray data---stem cell stages
#  (home made)
#
# ####################################################################

#==================Loading annotation data (essential)---stem cell stages============
annotation(eData)
platf <- getGEO(annotation(eData), AnnotGPL=TRUE)
anot <- data.frame(attr(dataTable(platf), "table"))
anot1<- str_split_fixed(anot[,"Gene.symbol"], "///",2) #split column 10
anno<-cbind.data.frame(anot[,"ID"],anot1[,1])
anno <- anno[!anno[,2]=="",] #remove empty columns with no gene names
colnames(anno) <-c("PROBES","gene_ID")
#===========================================================================#


#  ----------Loading annotation data (locally,alternative)---stem cell stages----------

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
# -----------------------------------------------------------------


#===============#replace affymetrix ID to gene ID (essential)==================##
PROBES <- as.numeric(rownames(exprs(eData)))
MicroArrayData<- cbind.data.frame(PROBES,exprs(eData))

anno_MicroArrayData <- merge(anno, MicroArrayData, all = TRUE)
anno_MicroArrayData <- anno_MicroArrayData[!is.na(anno_MicroArrayData[,3]),] #remove NA row in exprs
anno_MicroArrayData <- anno_MicroArrayData[!is.na(anno_MicroArrayData[,"gene_ID"]),] #remove NA row in genes

dim(anno_MicroArrayData)
anno_MicroArrayData <- anno_MicroArrayData[order(anno_MicroArrayData[,"gene_ID"]),]
anno_MicroArrayData <- distinct(anno_MicroArrayData,gene_ID,.keep_all = TRUE) #remove duplicated rows
dim(anno_MicroArrayData)

rownames(anno_MicroArrayData) <- anno_MicroArrayData[,"gene_ID"]
anno_MicroArrayData <- anno_MicroArrayData[,-1] #remove "PROBES"

# ===========get gene list that have higher std (top 25%) (essential)===========
# we can find the standard deviation for each chip
# http://jura.wi.mit.edu/bio/education/bioinfo2007/arrays/array_exercises_1R.html#ratios

sd = apply(anno_MicroArrayData[,-1], 1, sd) #run sd w/o gene_ID
sd = sd[order(sd, decreasing=TRUE)]
Top_Std_gene <-names(sd)[1:(length(sd)*0.25)]
# ===========================================================================#


# ------------visulization (optional)-------------------------------------------------------------
median  = apply(anno_MicroArrayData[,-1], 1, median)
par(mfrow=c(1,1))
plot(median,sd,xlab="log2(count size)", ylab="standard deviation")
lines(lowess(median,sd), col="blue")# lowess line (x,y)







# ######################################################################
# 
#  1.c loading Expression data---rnaseq_rowcounts
#  consult to http://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html
#
# ######################################################################

#======================Reading in the data (essential)======================
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
if (Sys.info()[['sysname']]=="Darwin"){
        setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}

countdata <- read.csv("rnaseq_rowcounts.csv")
colnames(countdata)[1] <- "gene_ID"
rownames(countdata) <- countdata[,"gene_ID"]
countdata<-countdata[,-1]
#Filtering to remove lowly expressed genes
# Obtain CPMs
myCPM <- cpm(countdata)
# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 1
table(rowSums(thresh) >= 3)
# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 3
# Subset the rows of countdata to keep the more highly expressed genes
myCPM.keep <- myCPM[keep,]

#creat DGE------------------------
expr0 <- DGEList(countdata)
group <- as.factor(c("E8","E8","E8","PGL","PGL","PGL","SCIL","SCIL","SCIL"))
expr0$samples$group <- group

expr1 <- expr0
expr1$counts <-myCPM.keep

# Get log2 counts per million
logcounts1 <- cpm(expr0,log=TRUE)
logcounts2 <- cpm(expr1,log=TRUE)
#normalising gene expression distributions
expr1 <- calcNormFactors(expr1, method ="TMM")
expr1$samples
logcounts3 <- cpm(expr1,log=TRUE)
#===========================================================================#



# ######################################################################
# 
#  1.d prepare Expression data---rnaseq_rowcounts
#  consult to Law c 's F1000Research 2016 http://dx.doi.org/10.12688/f1000research.9005.1
#
# ######################################################################

#===========plot for the low reads range===========
par(mfrow=c(1,2))
par(oma=c(2,2,2,2))

#
#plot(counts0,cex=.5)
plot(logcounts1[,1:2], cex=.5)
title("raw data")
plot(logcounts3[,1:2], cex=.5) 
title("filterd data")


#===========density polt===========
library(RColorBrewer)
nsamples <- ncol(logcounts2)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2),cex=1.3)
plot(density(logcounts1[,1]), col=col[1], lwd=2, ylim=c(0,0.4), las=2,
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
        den <- density(logcounts1[,i])
        lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(logcounts1), text.col=col, bty="n")


plot(density(logcounts3[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
        den <- density(logcounts2[,i])
        lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(logcounts2), text.col=col, bty="n")

#===========boxplot (raw vs normalised)===========
par(mfrow=c(1,3))
par(oma=c(5,2,2,2))
#distributions of samples using boxplots
# The las argument rotates the axis names

boxplot(logcounts1, xlab="", ylab="Log2 counts per million",las=2) # test with all((logcounts1+5.2)>0)
abline(h=median(logcounts1),col="blue")
title("logCPMs (raw)")

boxplot(logcounts2, xlab="", ylab="Log2 counts per million",las=2) # test with all((logcounts2+2)>0)
abline(h=median(logcounts2),col="blue")
title("logCPMs (trimed,un-normalisd)")


boxplot(logcounts3, xlab="", ylab="Log2 counts per million",las=2) # test with all((logcounts2+2)>0)
abline(h=median(logcounts2),col="blue")
title("logCPMs (normalised)")

#===========Unsupervised clustering of samples===========
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(logcounts2, labels=group, col=col.group)
title(main="A. Sample groups")
plotMDS(logcounts3, labels=group, col=col.group)
title(main="B. Sample groups")

#-----Visulaze Mean-variance relationships (optional)-----------------
sd      = apply(logcounts2, 1, sd)
median  = apply(logcounts2, 1, median)
par(mfrow=c(1,1))

plot(median,sd,xlab="log2(count size)", ylab="standard deviation")
lines(lowess(median,sd), col="blue")# lowess line (x,y)

######################################################################
# 
#  1.e merge rnaseq with MicroArrayData
#  (home made)
#######################################################################

#=============    merge 1.a and 1.b (essential)======================
#put gene column back-
logcounts <- cbind.data.frame(rownames(logcounts3),logcounts3) 
colnames(logcounts)[1] <- "gene_ID"


rnaseq_MicroArrayData <- merge(logcounts, anno_MicroArrayData)
dim(rnaseq_MicroArrayData)
rownames(rnaseq_MicroArrayData) <- rnaseq_MicroArrayData[,"gene_ID"]
rnaseq_MicroArrayData <- rnaseq_MicroArrayData[,-1] #remove "gene_ID"
write.csv(rnaseq_MicroArrayData,"rnaseq_MicroArrayData.csv")




#===========get gene list that have higher std (top 25%) (essential)===========
all_gene <- rownames(rnaseq_MicroArrayData)
rnaseq_MicroArrayData_std <-rnaseq_MicroArrayData[all_gene %in% Top_Std_gene,]

par(mfrow=c(1,1))
par(oma=c(16,2,0,2), cex=1.2)
boxplot(rnaseq_MicroArrayData_std, xlab="", ylab="log2 counts per million",
        main = "boxplot for log Counts Per Million of all samples",las=2) # test with all((counts0)>=0)

#rnaseq_MicroArrayData[Top_Std_gene,] will creat NA

dim(rnaseq_MicroArrayData_std)
write.csv(rnaseq_MicroArrayData_std,"rnaseq_MicroArrayData_std.csv")
#==========================================================================#

######################################################################
# 
#  1.e merge rnaseq with MicroArrayData
#  (home made)

#======Visulaze Mean-variance relationships (optional)===========
par(mfrow=c(1,2))
sd1      = apply(rnaseq_MicroArrayData[,10:ncol(rnaseq_MicroArrayData_std)], 1, sd)
median1  = apply(rnaseq_MicroArrayData[,10:ncol(rnaseq_MicroArrayData_std)], 1, median)

sd2      = apply(rnaseq_MicroArrayData[,1:9], 1, sd)
median2  = apply(rnaseq_MicroArrayData[,1:9], 1, median)


plot(median1,sd1,xlab="log2(count size)", ylab="standard deviation",
     main = "Mean-variance relationships for Mircroarray",las=2)
lines(lowess(median1,sd1), col="blue")# lowess line (x,y)
abline(h=median(sd1[all_gene %in% Top_Std_gene]),col="red")



plot(median2,sd2,xlab="log2(count size)", ylab="",
     main = "Mean-variance relationships for RNA-seq",las=2)
lines(lowess(median2,sd2), col="blue")# lowess line (x,y)
abline(h=median(sd2[all_gene %in% Top_Std_gene]),col="red")


#----find out highly expression genes (optional)-----------
#we can find the median for each RNA-Seq
median = apply(rnaseq_MicroArrayData[,1:9], 1, median)
median = median[order(median, decreasing=TRUE)]
Top_expressed_genes <-names(median)
par(oma=c(2,2,2,2), cex=1.3)
boxplot(t(rnaseq_MicroArrayData[Top_expressed_genes[1:20],1:9]),
        ylab="log2 counts per million",main = "Top 15 highly expressed genes in RNA-seq data",las=2)

# ------------visulization filteration (optional)-------------------------------------------------------------
#check the low reads range
par(mfrow=c(1,2))
par(oma=c(2,2,2,2))

plot(rnaseq_MicroArrayData[,1:2],cex=.5)
plot(rnaseq_MicroArrayData_std[,1:2],cex=.5)
