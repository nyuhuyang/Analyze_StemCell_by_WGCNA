#compare sets of transcriptional profiles in several ways:
#    - by overall similarity eg using Spearman correlation
#    - by similarity of signatures - usually done using GSEA (my preferred approach)
#    - by first computing modules (using WCGNA) in each dataset,
#        then comparing overlap of modules using hypergeometric test

########################################################################
#  0 preparation and parameter adjustion
# 
# ######################################################################

list.of.cran.packages<- c("easypackages","ggfortify","grid","pheatmap")
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
#  4. overall similarity using Spearman correlation
#  
#######################################################################
#====4.1 Reading in the data (required)===================
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
for(i in c(2:9)){ matrix1 <- read.csv(as.character(list_files$dataset[i]))
matrix <- merge(matrix,matrix1,by ="X") }
dim(matrix)
head(matrix[,1:3])
rownames(matrix) <-matrix[,1]
matrix <- matrix[,-1]
#==normalize
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

#========extract 5000 most variant genes=============
#extract the top 5000 most variant genes for Correlation studies.
#transpose matrix to correlate genes in the following
matrix_5000 = matrix.scale[order(apply(matrix.scale,1,mad), decreasing = T)[1:5000],]

#===subset and group====
sample_names <- read.csv("sample_names.csv",row.names = 1,header = T)
naive_cells <- sample_names[sample_names$stage=="naive",c("cell_names","authors","stages","file_names")]
naive_matrix <- matrix_5000[,as.numeric(rownames(naive_cells))] # keep the order
colnames(naive_matrix) <- naive_cells$cell_names
#====4.2 Calculate Spearman correlation (required)=================================
c <- cor(naive_matrix, method="spearman")
diag(c) <-NA
colnames(c) <- naive_cells$cell_names
rownames(c) <- naive_cells$cell_names

pheatmap(c,cex=.8,
         cluster_rows=F,
         cluster_cols = F,
         fontsize_row = 10,
         fontsize_col = 10,
         fontsize =30,
         main ="Pre-implantation correlation by Spearman")
#====4.3 Comparision between dataset (required)=================================
#====4.3-a Vassena vs.Xie==========
pheatmap(c[1:27,28:48], 
         cluster_rows=F,
         cluster_cols = F,
         fontsize_row = 15,
         fontsize_col = 15,
         fontsize =15,
         main ="Vassena et al 2011 vs Xie et al 2010")
#====4.3-b Vassena vs. All other labs ==========
#====4.3-b1 make annotation_col and gap_col============
annotation_col <- data.frame(database =naive_cells$authors[49:ncol(c)])
rownames(annotation_col) <- naive_cells$cell_names[49:ncol(c)]
annotation_col <-droplevels(annotation_col)

gaps_col <- as.data.frame(table(annotation_col))
gaps_col
gaps_col <- gaps_col[c(5,6,3,7,1,4,2),] #switch sequence
gaps_col$gap <-NA #insert empty column
gaps_col$gap[1] <-gaps_col$Freq[1]
for(i in 2:length(gaps_col$gap)){
        gaps_col$gap[i] <- gaps_col$gap[i-1] + gaps_col$Freq[i]}
gaps_col
#====4.3-b2 pheatmap===========
#Vassena 2011 against all
pheatmap(c[1:27,49:ncol(c)], 
         cluster_rows=F,
         cluster_cols = F,
         fontsize_row = 15,
         fontsize_col = 11,
         fontsize =15,
         gaps_col=gaps_col$gap,
         annotation_col = annotation_col,
         main ="Spearman correlation: Vassena 2011 against all samples")

#Vassena 2011 against Danwei
upside_down <- unlist(t(c[27:1,93:ncol(c)]))
pheatmap(c[1:27,93:ncol(c)], 
         cluster_rows=F,
         cluster_cols = F,
         fontsize_row = 15,
         fontsize_col = 15,
         fontsize =15,
         #         gaps_col=gaps_col$gap,
         #         annotation_col = annotation_col,
         main ="Spearman correlation: Vassena 2011 against Danwei")
grid.text(round(upside_down,digits =3),
          x=seq(0.06, 0.76, length.out=9), 
          y=rep(seq(0.07, 0.93, 0.032)+0.032, each=9)) #rep =length(X)

#Xie 2010 against All
pheatmap(c[28:48,49:ncol(c)],
         cluster_rows=F,
         cluster_cols = F,
         fontsize_row = 15,
         fontsize_col = 11,
         fontsize =15,
         gaps_col=gaps_col$gap,
         annotation_col = annotation_col,
         main ="Spearman correlation: xie 2010 against all samples")

upside_down <- unlist(t(c[48:28,93:ncol(c)]))
pheatmap(c[28:48,93:ncol(c)],
         cluster_rows=F,
         cluster_cols = F,
         fontsize_row = 15,
         fontsize_col = 15,
         fontsize =15,
         #         gaps_col=gaps_col$gap,
         #         annotation_col = annotation_col,
         main ="Spearman correlation: Xie 2010 against Danwei")
grid.text(round(upside_down,digits =3),
          x=seq(0.06, 0.76, length.out=9), 
          y=rep(seq(0.07, 0.93, 0.042)+0.032, each=9)) #rep =length(X)


#==4.4 average pheatmap===========
c <- cor(naive_matrix, method="spearman") #repeat to avoid NA
colnames(c)
EScell_name <-droplevels(unique(naive_cells$cell_names[1:48])) #average_c rowname
EScell_papers_name <-c(as.character(EScell_name),
                       as.character(list_files$name[-c(1:2)])) #average_c colname
average_c <- data.frame(matrix(0, #create 0 matrix
                               nrow = length(EScell_name), #number of EScell name
                               ncol = length(EScell_papers_name))) #number of papers
rownames(average_c) <- EScell_name
colnames(average_c) <- EScell_papers_name
average_c

for(row in 1:nrow(average_c)){
        for(col in 1:length(EScell_name)){
                average_c[row,col] <-mean(c[rownames(c)==rownames(average_c)[row],
                                            colnames(c)==colnames(average_c)[col]])}
}
for(row in 1:nrow(average_c)){
        for(col in (length(EScell_name)+1):ncol(average_c)){
                average_c[row,col] <-mean(c[rownames(c)==EScell_name[row],
                                            naive_cells$authors==EScell_papers_name[col]])}
}
average_c

upside_down <- unlist(t(average_c[9:1,(length(EScell_name)+1):ncol(average_c)]))
pheatmap(average_c[1:9,(length(EScell_name)+1):ncol(average_c)],
         cluster_rows=F,
         cluster_cols = F,
         fontsize_row = 20,
         fontsize_col = 20,
         fontsize =15,
         main ="Spearman correlation: Vassena 2011 against all samples(Average)")
grid.text(round(upside_down,digits =3),
          x=seq(0.06, 0.73, length.out=7), 
          y=rep(seq(0.28, 0.95, 0.076)+0.03, each=7))


upside_down <- unlist(t(average_c[nrow(average_c):10,(length(EScell_name)+1):ncol(average_c)]))
pheatmap(average_c[10:nrow(average_c),(length(EScell_name)+1):ncol(average_c)],
         cluster_rows=F,
         cluster_cols = F,
         fontsize =20,
         main ="Spearman: Xie 2010 against all samples(Average)")
grid.text(round(upside_down,digits =3),
          x=seq(0.06, 0.71, length.out=7), 
          y=rep(seq(0.28, 0.89, 0.098)+0.03, each=7))



