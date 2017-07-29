
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
list.of.bio.packages <- c("DESeq2")
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)


libraries(list.of.bio.packages,list.of.cran.packages)


# ######################################################################
# 
#  2.1 Load RNA-seq data and pre-process --DESeq (Recommend)
#  http://genomicsclass.github.io/book/pages/rnaseq_gene_level.html
# ######################################################################
#=========Reading in the data (required)===================
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
if (Sys.info()[['sysname']]=="Darwin"){
        setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}

countdata <- read.csv("rnaseq_rowcounts.csv")
rownames(countdata) <- countdata[,1]
countdata<-countdata[,-1]
#create dds
group <- data.frame(group = c("E8","E8","E8","PGL","PGL","PGL","SCIL","SCIL","SCIL"))
dds <- DESeqDataSetFromMatrix(countdata,
                              colData = group,
                              design  = ~group) # must have ~
dds
dds <- estimateSizeFactors(dds)
rld <- rlog(dds, blind=FALSE)
rld
#filter <0 counts(Alternative)----
rs <- rowSums(assay(rld)) 
rld <-rld[rs > 0,] 
rld
#=============
head(assay(rld))
saveRDS(rld,"rld")  #===> next page

# ######################################################################
# 
#  2.2 RNA-seq DESeq visulization (Recommend)
# #http://genomicsclass.github.io/book/pages/rnaseq_gene_level.html
# ######################################################################

#---2.2.1 Quality Control (Alternative)---------------------
#--test the geometric mean of all samples(Alternative)----------------
par(mfrow=c(1,1))
loggeomeans <- rowMeans(log(counts(dds)))
hist(log(counts(dds)[,1]) - loggeomeans, 
     col="grey", main="geometric mean of all samples", xlab="", breaks=40)

#--Normalization for sequencing depth(Alternative)----------------
sizeFactors(dds)
colSums(counts(dds))
plot(sizeFactors(dds), colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))

#--Examine not normalized and normalized samples(Alternative)---------
rs <- rowSums(counts(dds)) #filter
log.norm <- normTransform(dds) #class "DESeq2"
log.norm.counts <- log2(counts(dds, normalized=TRUE) + 1) #class "matrix"

par(mfrow=c(1,3))
boxplot(log2(counts(dds)[rs > 0,] + 1)) # not normalized
boxplot(log2(counts(dds, normalized=TRUE) + 1)[rs > 0,]) # normalized
boxplot(log2(assay(log.norm)[rs > 0,] + 1)) # normalized 
#The same matrix as above is stored in assay(log.norm).


#====Stabilizing count variance(required)=======================
par(mfrow=c(1,3))
plot(log.norm.counts[,1:2], cex=.1)
#rlog perform better when the size factors vary widely
plot(assay(rld)[,1:2], cex=.1)
#varianceStabilizingTransformation is much faster
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
plot(assay(vsd)[,1:2], cex=.1)


#===standard deviation of rows over the mean (required)==========
library(vsn)
# genes with high variance for the log come from the genes with lowest mean
meanSdPlot(log.norm.counts, ranks=FALSE) 
#For the rlog:
meanSdPlot(assay(rld), ranks=FALSE)
#For the VST:
meanSdPlot(assay(vsd), ranks=FALSE)

#===examining relationships between samples (Recommend)==========
plotPCA(log.norm, intgroup="group")
#Using the rlog:
plotPCA(rld, intgroup="group")
#Using the VST:
plotPCA(vsd, intgroup="group")

# for nicer plot
(data <- plotPCA(rld, intgroup=c("group"), returnData=TRUE)) 
(percentVar <- 100*round(attr(data, "percentVar"),2))
makeLab <- function(x,pc) paste0("PC",pc,": ",x,"% variance")
ggplot(data, aes(PC1,PC2,col=dex,shape=cell)) + geom_point() +
        xlab(makeLab(percentVar[1],1)) + ylab(makeLab(percentVar[2],2))


# ######################################################################
# 
#  2.3 loading Expression data and pre-processing using dge (Alternative)
#  consult to http://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html
#             https://www.bioconductor.org/help/workflows/RNAseq123/
#             https://bioconductor.org/packages/devel/bioc/vignettes/DEFormats/inst/doc/DEFormats.html
# 
# ######################################################################

#----------------------Reading in the data (required)----------------------
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
if (Sys.info()[['sysname']]--"Darwin"){
        setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
if (Sys.info()[['sysname']]--"Windows"){
        setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}

countdata <- read.csv("rnaseq_rowcounts.csv")
colnames(countdata)[1] <- "gene_ID"
rownames(countdata) <- countdata[,"gene_ID"]
countdata<-countdata[,-1]
#Filtering to remove lowly expressed genes
library(edgeR)
# Obtain CPMs

#----filtering----=
myCPM <- cpm(countdata)
# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 1
table(rowSums(thresh) >= 3)
# we would like to keep genes that have at least 3 TRUES in each row of thresh
keep <- rowSums(thresh) >= 3
# Subset the rows of countdata to keep the more highly expressed genes
myCPM.keep <- myCPM[keep,]

#--creat DGE------------
dge0 <- DGEList(countdata)
group <- as.factor(c("E8","E8","E8","PGL","PGL","PGL","SCIL","SCIL","SCIL"))
dge0$samples$group <- group

dge1 <- dge0
dge1$counts <-myCPM.keep #filter

#--normalization--------
dge2 <- calcNormFactors(dge1, method ="TMM") #normalized
dge2$samples

# Get log2 counts per million
logcounts0 <- cpm(dge0,log=TRUE) # raw for boxplot only
logcounts1 <- cpm(dge1,log=TRUE) # trimed for boxplot only
logcounts2 <- cpm(dge2,log=TRUE) # trimed and normalized

#--Convert DGEList to DESeqDataSet(Required)--------------
library(DEFormats)
dds = as.DESeqDataSet(dge0)


# ######################################################################
# 
#  2.4 RNA-seq deg visulization
# consult to Law c 's F1000Research 2016 http://dx.doi.org/10.12688/f1000research.9005.1
# ######################################################################

#--plot for the low reads range(recommended)----=
par(mfrow=c(1,2))
par(oma=c(2,2,2,2))

#
#plot(counts0,cex=.5)
plot(logcounts0[,1:2], cex=.5)
title("raw data")
plot(logcounts2[,1:2], cex=.5) 
title("filterd data")


#--density polt(recommended)----=
library(RColorBrewer)
nsamples <- ncol(logcounts0)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2),cex=1.3)
plot(density(logcounts0[,1]), col=col[1], lwd=2, ylim=c(0,0.4), las=2,
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
        den <- density(logcounts0[,i])
        lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(logcounts0), text.col=col, bty="n")


plot(density(logcounts2[,1]), col=col[1], lwd=2, ylim=c(0,0.21), las=2,
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
        den <- density(logcounts2[,i])
        lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", colnames(logcounts2), text.col=col, bty="n")

#--boxplot (raw vs normalised)(recommended)----------=
par(mfrow=c(1,3))
par(oma=c(5,2,2,2))
#distributions of samples using boxplots
# The las argument rotates the axis names

boxplot(logcounts0, xlab="", ylab="Log2 counts per million",
        las=2,col=col,main="logCPMs (raw)") # test with all((logcounts1+5.2)>0)
abline(h=median(logcounts0),col="blue")


boxplot(logcounts1, xlab="", ylab="Log2 counts per million",
        las=2, col=col, main="logCPMs (trimed,un-normalisd)") # test with all((logcounts2+2)>0)
abline(h=median(logcounts1),col="blue")



boxplot(logcounts2, xlab="", ylab="Log2 counts per million",
        las=2,col=col, main="logCPMs (normalised)") # test with all((logcounts2+2)>0)
abline(h=median(logcounts2),col="blue")


#--Unsupervised clustering of samples(recommended)------------
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <- brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(logcounts1, labels=group, col=col.group)
title(main="A. Sample groups")
plotMDS(logcounts2, labels=group, col=col.group)
title(main="B. Sample groups")

#-----Visulaze Mean-variance relationships (Required)-----------------
par(mfrow=c(1,2))

sd1      = apply(logcounts2, 1, sd) 
median1  = apply(logcounts2, 1, median)

sd2      = apply(rld.matrix, 1, sd)
median2  = apply(rld.matrix, 1, median)


plot(median1,sd1,xlab="log2(count size)", ylab="standard deviation")
lines(lowess(median1,sd1), col="blue")# lowess line (x,y)


plot(median2,sd2,xlab="log2(count size)", ylab="standard deviation")
lines(lowess(median2,sd2), col="blue")# lowess line (x,y)

# ######################################################################
# 
#  2.5 Differential expression analysis --voom (Alternative)
# # https://www.bioconductor.org/help/workflows/arrays/
# Page 70 https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf
# https://www.bioconductor.org/help/workflows/RNAseq123/
# ######################################################################

#Creating a design matrix (required)------------=
design <- model.matrix(~0+group) # describe model to be fit
colnames(design) <- gsub("group", "", colnames(design))
design

# limma-trend (Alternative)-------------
#simplest and most robust approach to differential exis to use limma-trend
fit <- lmFit(logcounts2, design)  # fit each probeset to model
efit <- eBayes(fit)        # empirical Bayes adjustment
topTable(efit, coef=2)      # table of differentially expressed probesets

# voom (recommended)----------------=
#When the library sizes are quite variable between samples,
#then the voom approach is theoretically more powerful than limma-trend.

#Creating a contrasts matrix
contr.matrix <- makeContrasts(
        E8vsPGL = E8-PGL, 
        E8vsSCIL = E8 - SCIL, 
        PGLvsSCIL = PGL - SCIL, 
        levels = colnames(design))
contr.matrix

#Removing heteroscedascity from count data
par(mfrow=c(1,2))
v <- voom(dge1, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean???variance trend")
# I can't extract count matrix from voom object. Give up
# have to use DESeq

