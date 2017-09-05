#compare sets of transcriptional profiles in several ways:
#    - by overall similarity eg using Spearman correlation
#    - by similarity of signatures - usually done using GSEA (my preferred approach)
#    - by first computing modules (using WCGNA) in each dataset,
#        then comparing overlap of modules using hypergeometric test
########################################################################
#
#  0 preparation and parameter adjustion
# 
# ######################################################################

list.of.cran.packages<- c("easypackages","WGCNA")
new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(easypackages)

source("https://bioconductor.org/biocLite.R")
list.of.bio.packages <- c()
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)


libraries(list.of.bio.packages,list.of.cran.packages)


########################################################################################
# 
#  6. GSEA (required)
#  http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
# 
########################################################################################

# =====6.1 load csv and prepare clusters (required)==================================================
# =====6.1.1 load csv (required)==================================================
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
if (Sys.info()[['sysname']]=="Darwin"){
        setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}

list_files <- read.csv("list_files.csv",header = T)
list_files
matrix <- read.csv(as.character(list_files$dataset[1]))
#prime_hESC <- read.csv("./Huang_Cell_2014/prime_hESC.csv")
#matrix <- merge(matrix,prime_hESC,by ="X")
#write.csv(matrix,"Xie_dataset.csv")
matrix1 <- read.csv(as.character(list_files$dataset[9]))
matrix <- merge(matrix,matrix1,by ="X") 
dim(matrix)
head(matrix[,1:3])
rownames(matrix) <-matrix[,1]
matrix <- matrix[,-1]
#==normalize=======
matrix.scale <- scale(matrix)
# check that we get mean of 0 and sd of 1
allmean <- mean(colMeans(matrix))  # faster version of apply(scaled.dat, 2, mean)
matrix.scale <- matrix.scale +allmean 

#boxplot
par(mfrow=c(1,2))
boxplot(matrix)
title("Boxplot Before Normalization")
boxplot(matrix.scale)
title("boxplot after Normalization")



#====6.1.2 Select test group (Required)=================================
#===extract primed hESCs and Danwei's data ====
sample_names <- read.csv("sample_names.csv",row.names = 1,header = T)
#sample_names <-sample_names[c(25:27,147:155),c("cell_names","file_names")]
sample_names <-sample_names[c(147:155),c("cell_names","file_names")]
matrix.selected <-matrix.scale[,as.character(sample_names$file_names)]
colnames(matrix.selected)
#====6.1.3 Prepare GSEA files (Required)=========

# ===Create txt Expression Data for GSEA====
GSEA <- as.data.frame(matrix.selected[,1:2])# didn't works for [,1]
GSEA[,1]<-NA
GSEA_exp<- cbind.data.frame(GSEA[,1],matrix.selected)
colnames(GSEA_exp)[1] <-"DESCRIPTION" #required

write.table(GSEA_exp, "GSEA_exp_less.txt", sep="\t")
# Insert "NAME" at [1,1] using excel

# ===Create phentype labels for GSEA====

ncol<-ncol(matrix.selected)
#GSEA.cls = data.frame(c(ncol,4,1,rep("",ncol-3)),
#                 c("#","hESCs","E8","PGL","SCIL",rep("",ncol-5)),
#                 c(rep("hESCs",3),
#                   rep("E8",3),
#                   rep("PGL",3),
#                   rep("SCIL",3)),
#                 fix.empty.names = FALSE)
GSEA.cls = data.frame(c(ncol,3,1,rep("",ncol-3)),
                 c("#","E8","PGL","SCIL",rep("",ncol-4)),
                 c(rep("E8",3),
                   rep("PGL",3),
                   rep("SCIL",3)),
                 fix.empty.names = FALSE)
write.table(t(GSEA.cls), "sample_less.table.cls", sep=" ")
# open with texteditor
# delete the first row;
# delete all "
# delete the space before and after each row.
# http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats

#====> next section










#====6.1.2 Select test group (Required)=================================
#====6.1.2-a Calculate Spearman correlation (obsolete)=================================
#========extract 5000 most variant genes=============
#extract the top 5000 most variant genes for Correlation studies.
#transpose matrix to correlate genes in the following
#matrix_5000 = matrix.scale[order(apply(matrix.scale,1,mad), decreasing = T)[1:5000],]

#===extract primed hESCs and Danwei's data ====
sample_names <- read.csv("sample_names.csv",row.names = 1,header = T)
#sample_names <-sample_names[c(1:27,147:155),c("cell_names","file_names")]
sample_names <-sample_names[c(25:27,147:155),c("cell_names","file_names")]
matrix.selected <-matrix[,as.character(sample_names$file_names)]
#====Spearman correlation (required)=================================
c <- cor(matrix.selected, method="spearman")
diag(c) <-NA
colnames(c) <- sample_names$cell_names
rownames(c) <- sample_names$cell_names

#====6.1.3 Prepare GSEA files (Alternative)=================================

#====hcluster====
par(mfrow=c(1,1))
par(oma=c(2,2,2,2))
d <- as.dist(1-c)
hc <- hclust(d, method="complete")
p<-plot(hc,col = "black",cex = 1,xlab="", sub="",
        #        main="Cluster Dendrogram for Vassane and Danwei samples")
        main="Cluster Dendrogram for hESCs and Danwei samples")
# Plot a line to show the cut
abline(h = 0.010, col = "red");
clust = cutreeStatic(hc, cutHeight = 0.011, minSize = 3)
table(clust)

#----6.1.3-b Select group by cluster (Alternative)----
clust = cutreeStatic(hc, cutHeight = 0.011, minSize = 3)
table(clust)
cluster <- data.frame(labels=hc$labels,
                      clust=clust,
                      order=hc$order,
                      class=NA)
cluster[cluster$clust==1,"class"]<-"primed hESCs"
cluster[cluster$clust==2,"class"]<-"E8"
cluster[cluster$clust==3,"class"]<-"PBL"
cluster[cluster$clust==4,"class"]<-"SCIL"


#define cls 
cls <-function(class){
        return(as.character(cluster[cluster[,"class"]==class,"labels"]))
}
cls("primed hESCs") #test

# =====5.1 Preparing Data Files for GSEA(required)==================================================
# Preparing Data Files for GSEA
# Create txt Expression Data for GSEA

rnaseq_MicroArrayData_GSEA <- matrix.scale[,1:2]# didn't works for [,1]
rnaseq_MicroArrayData_GSEA[,"DESCRIPTION"]<-NA

test_group <-"primed hESCs"      #clust==1
#test_group <-"E8"                 #clust==2
#test_group <-"PBL"            #clust==3
#test_group <-"SCIL"               #clust==4

rnaseq_MicroArrayData_GSEA<- cbind.data.frame(
        rnaseq_MicroArrayData_GSEA[,"DESCRIPTION"],
        rnaseq_MicroArrayData[,cls("primed hESCs")],
        rnaseq_MicroArrayData[,cls(test_group)] #change test_group
)
colnames(rnaseq_MicroArrayData_GSEA)[1] <-"DESCRIPTION"


write.table(rnaseq_MicroArrayData_GSEA, "rnaseq_MicroArrayData_GSEA.txt", sep="\t")
# Insert "NAME" at [1,1] using excel
#===========================================
# Create phentype labels for GSEA

ncol<-ncol(rnaseq_MicroArrayData_GSEA)
df3 = data.frame(c(ncol-1,2,1,rep("",ncol-4)),
                 c("#","Danwei",test_group,rep("",ncol-4)),#change test_group
                 c(rep("Danwei",sum(cluster$class=="Danwei")),
                   rep(test_group,sum(cluster$class==test_group))),#change test_group
                 fix.empty.names = FALSE)
write.table(t(df3), "sample.table.cls", sep=" ")
#open with texteditor
#delete the first row;
#delete all "
#delete the space before each row.
#http://recipes.genomespace.org/view/15#collapse_5

#===> next page


