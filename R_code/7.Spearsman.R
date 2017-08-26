######################################################################
# 
#  7.1 overall similarity using Spearman correlation
#  (home made)
#######################################################################
#====7.1.1 Reading in the data (required)===================
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
#===subset and group====
sample_names <- read.csv("sample_names.csv",row.names = 1,header = T)
naive_cells <- sample_names[sample_names$stage=="naive",c("cell_names","authors","stages","file_names")]
naive_matrix <- matrix[,as.numeric(rownames(naive_cells))] # keep the order
colnames(naive_matrix) <- naive_cells$cell_names
#====7.1.2 Spearman correlation (required)=================================
c <- cor(naive_matrix, method="spearman")
diag(c) <-NA
colnames(c) <- naive_cells$cell_names
rownames(c) <- naive_cells$cell_names
library(pheatmap)
pheatmap(c,cex=.8,
         cluster_rows=F,
         cluster_cols = F,
         fontsize_row = 10,
         fontsize_col = 10,
         fontsize =25,
         main ="Pre-implantation correlation by Spearman")
#====7.1.3 Spearman correlation smaller size (required)=================================
#=========make annotation_col and gap_col============
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
#pheatmap===========
pheatmap(c[1:27,49:ncol(c)],
         cluster_rows=F,
         cluster_cols = F,
         fontsize_row = 15,
         fontsize_col = 11,
         fontsize =15,
         gaps_col=gaps_col$gap,
         annotation_col = annotation_col,
         main ="Pre-implantation correlation against Vassena et al 2011")

pheatmap(c[28:48,49:ncol(c)],
         cluster_rows=F,
         cluster_cols = F,
         fontsize_row = 15,
         fontsize_col = 11,
         fontsize =15,
         gaps_col=gaps_col$gap,
         annotation_col = annotation_col,
         main ="Pre-implantation correlation against Xie et al 2010")
#average pheatmap===========
c <- cor(naive_matrix, method="spearman") #repeat to avoid NA
colnames(c)
EScell_name <-droplevels(unique(naive_cells$cell_names[1:48])) #average_c rowname
EScell_papers_name <-c(as.character(EScell_name),
                       as.character(list_files$name[-c(1:2)])) #average_c colname
average_c <- data.frame(matrix(NA, 
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
library(grid)
upside_down <- unlist(t(average_c[9:1,(length(EScell_name)+1):ncol(average_c)]))
pheatmap(average_c[1:9,(length(EScell_name)+1):ncol(average_c)],
         cluster_rows=F,
         cluster_cols = F,
         fontsize_row = 20,
         fontsize_col = 20,
         fontsize =15,
         main ="Pre-implantation correlation against Vassena at al 2011")
grid.text(round(upside_down,digits =3),
          x=seq(0.06, 0.73, length.out=7), 
          y=rep(seq(0.28, 0.95, 0.076)+0.03, each=7))


upside_down <- unlist(t(average_c[nrow(average_c):10,(length(EScell_name)+1):ncol(average_c)]))
pheatmap(average_c[10:nrow(average_c),(length(EScell_name)+1):ncol(average_c)],
         cluster_rows=F,
         cluster_cols = F,
         fontsize =20,
         main ="Pre-implantation correlation against Xie et al 2010")
grid.text(round(upside_down,digits =3),
          x=seq(0.06, 0.71, length.out=7), 
          y=rep(seq(0.28, 0.89, 0.098)+0.03, each=7))
