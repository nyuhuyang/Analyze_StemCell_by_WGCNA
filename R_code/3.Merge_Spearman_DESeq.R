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
list.of.bio.packages <- c()
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)


libraries(list.of.bio.packages,list.of.cran.packages)

######################################################################
# 
#  3.1 merge rnaseq with MicroArrayData
#  (home made)
#######################################################################
#=========Reading in the data (required)===================
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
if (Sys.info()[['sysname']]=="Darwin"){
        setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}

#=============    merge microarray and rna-seq (required)======================
#put gene column back-
rld <- readRDS("rld")
rnaseq <- assay(rld)
anno <- read.csv("anno")
anno_MicroArrayData <- read.csv("anno_MicroArrayData",row.names = 1)

rnaseq <- cbind.data.frame(rownames(rnaseq),rnaseq) 
colnames(rnaseq)[1] <- "gene_ID"

dim(rnaseq)
dim(anno_MicroArrayData)

rnaseq_MicroArrayData <- merge(rnaseq, anno_MicroArrayData)
dim(rnaseq_MicroArrayData)
rownames(rnaseq_MicroArrayData) <- rnaseq_MicroArrayData[,"gene_ID"]
rnaseq_MicroArrayData <- rnaseq_MicroArrayData[,-1] #remove "gene_ID"
write.csv(rnaseq_MicroArrayData,"rnaseq_MicroArrayData.csv")




#===========get gene list that have higher std (top 25%) (required)===========
all_gene <- rownames(rnaseq_MicroArrayData)

sd = apply(anno_MicroArrayData[,-1], 1, sd) #run sd w/o gene_ID
sd = sd[order(sd, decreasing=TRUE)]
Top_Std_gene <-names(sd)[1:(length(sd)*0.25)]

rnaseq_MicroArrayData_std <-rnaseq_MicroArrayData[all_gene %in% Top_Std_gene,]

par(mfrow=c(1,1))
par(oma=c(16,2,0,2), cex=1.2)
boxplot(rnaseq_MicroArrayData_std, xlab="", ylab="log2 counts per million",
        main = "boxplot for log Counts Per Million of all samples",las=2) # test with all((counts0)>=0)

#rnaseq_MicroArrayData[Top_Std_gene,] will creat NA

dim(rnaseq_MicroArrayData_std)
write.csv(rnaseq_MicroArrayData_std,"rnaseq_MicroArrayData_std.csv")

#===> next page

######################################################################
# 
#  3.2 Explore data (Alternative)
#  (home made)

#======Visulaze Mean-variance relationships (Recommend)===========
sd1      = apply(rnaseq_MicroArrayData[,10:ncol(rnaseq_MicroArrayData_std)], 1, sd)
median1  = apply(rnaseq_MicroArrayData[,10:ncol(rnaseq_MicroArrayData_std)], 1, median)

sd2      = apply(rnaseq_MicroArrayData[,1:9], 1, sd)
median2  = apply(rnaseq_MicroArrayData[,1:9], 1, median)

par(mfrow=c(1,2))
par(oma=c(2,2,2,2), cex=1.2)

plot(median1,sd1,xlab="log2(count size)", ylab="standard deviation",
     main = "Mean vs. variance for Mircroarray",las=2)
lines(lowess(median1,sd1), col="blue")# lowess line (x,y)
abline(h=median(sd1[all_gene %in% Top_Std_gene]),col="red")



plot(median2,sd2,xlab="log2(count size)", ylab="",
     main = "Mean vs. variancefor RNA-seq",las=2)
lines(lowess(median2,sd2), col="blue")# lowess line (x,y)
abline(h=median(sd2[all_gene %in% Top_Std_gene]),col="red")


#----find out highly expression genes (Alternative)-----------
#we can find the median for each RNA-Seq
par(mfrow=c(1,1))
median = apply(rnaseq_MicroArrayData[,1:9], 1, median)
median = median[order(median, decreasing=TRUE)]
Top_expressed_genes <-names(median)
par(oma=c(2,2,2,2), cex=1.3)
boxplot(t(rnaseq_MicroArrayData[Top_expressed_genes[1:20],1:9]),
        ylab="log2 counts per million",main = "Top 15 highly expressed genes in RNA-seq data",las=2)

# ------------visulization filteration (Alternative)-------------------------------------------------------------
#check the low reads range
par(mfrow=c(1,2))
par(oma=c(2,2,2,2))

plot(rnaseq_MicroArrayData[,1:2],cex=.5)
plot(rnaseq_MicroArrayData_std[,1:2],cex=.5)
