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

list.of.cran.packages<- c("pheatmap")
new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(easypackages)

source("https://bioconductor.org/biocLite.R")
list.of.bio.packages <- c("WGCNA")
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)


libraries(list.of.bio.packages,list.of.cran.packages)

######################################################################
# 
#  4.1 overall similarity using Spearman correlation
#  (home made)
#######################################################################
#====4.1.1 Reading in the data (required)===================
# Display the current working directory
getwd();
# If necessary, change the path below to the directory where the data files are stored. 
# "." means current directory. On Windows use a forward slash / instead of the usual \.
if (Sys.info()[['sysname']]=="Darwin"){
        setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}

list_files <- read.csv("list_files.csv",row.names = 1)
# Read in rnaseq_MicroArrayData or rnaseq_MicroArrayData_std
rnaseq_MicroArrayData_std = read.csv("rnaseq_MicroArrayData_std.csv",row.names = 1);
#colnames are changed, be careful. space is replaced by .
# Take a quick look at what is in the data set:
dim(rnaseq_MicroArrayData_std)
head(rnaseq_MicroArrayData_std[,1:3])

par(mfrow=c(1,1))
par(oma=c(2,2,2,2))

#====4.1.2 Spearman correlation (required)=================================
c <- cor(rnaseq_MicroArrayData_std, method="spearman") # or rnaseq_MicroArrayData_std
diag(c) <-NA
pheatmap(c,cex=1.05,
#         cluster_rows=T,
         cluster_cols = T,
         main ="Spearman correlation between all samples")

########################################################################################
# 
#  5. GSEA (required)
#  http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
# 
########################################################################################

# =====5.1 load csv and prepare clusters (required)==================================================
#GSEA requires expression data of all genes
rnaseq_MicroArrayData = read.csv("rnaseq_MicroArrayData.csv",row.names = 1)
#colnames are changed, be careful; space is replaced by .
# Take a quick look at what is in the data set:
dim(rnaseq_MicroArrayData)
head(rnaseq_MicroArrayData[,11:13])

par(mfrow=c(1,1))
par(oma=c(2,2,2,2))
d <- as.dist(1-c)

hc <- hclust(d, method="complete")

p<-plot(hc,col = "black",cex = 1,xlab="", sub="",
        main="Cluster Dendrogram, method='complete'")
# Plot a line to show the cut
abline(h = 0.3, col = "red");
# Determine cluster under the line
clust = cutreeStatic(hc, cutHeight = 0.3, minSize = 3)
table(clust)
cluster <- data.frame(labels=hc$labels,
                      clust=clust,
                      order=hc$order,
                      class=NA)
cluster[cluster$clust==1,"class"]<-"Oocyte_to_6_Cells"
cluster[cluster$clust==2,"class"]<-"Danwei"
cluster[cluster$clust==3,"class"]<-"Morula_and_Blastocyst"
cluster[cluster$clust==4,"class"]<-"ES_Cell_line"
cluster[cluster$clust==5,"class"]<-"6~8_Cells"

#define cls 
cls <-function(class){
        return(cluster[cluster[,"class"]==class,"labels"])
}
cls("Oocyte_to_6_Cells") #test

# =====5.1 Preparing Data Files for GSEA(required)==================================================
# Preparing Data Files for GSEA
# Create txt Expression Data for GSEA
rnaseq_MicroArrayData_GSEA <- rnaseq_MicroArrayData[,1:2]# didn't works for [,1]
rnaseq_MicroArrayData_GSEA[,"DESCRIPTION"]<-NA

test_group <-"Oocyte_to_6_Cells"        #clust==1
#test_group <-"Danwei"                   #clust==2
#test_group <-"Morula_and_Blastocyst"    #clust==3
#test_group <-"ES_Cell_line"             #clust==4
#test_group <-"6~8_Cells"                #clust==5

rnaseq_MicroArrayData_GSEA<- cbind.data.frame(
        rnaseq_MicroArrayData_GSEA[,"DESCRIPTION"],
        rnaseq_MicroArrayData[,cls("Danwei")],
        rnaseq_MicroArrayData[,cls(test_group)] #change test_group
)
colnames(rnaseq_MicroArrayData_GSEA)[1] <-"DESCRIPTION"


write.table(rnaseq_MicroArrayData_GSEA, "rnaseq_MicroArrayData_GSEA.txt", sep="\t")
# Insert "NAME" at [1,1] using excel

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


