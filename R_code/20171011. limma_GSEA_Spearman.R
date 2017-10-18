########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
# 
# ######################################################################

source("./Dropbox/Public/Olivier/R/Functions.R")

initialize.project("Danwei_StemCell")
cran.libraries("scales","RColorBrewer","xlsx","stringr","pheatmap","ggplot2")
bio.libraries("edgeR","limma","Glimma","gplots")

# ######################################################################
# 
#  1 Differential expression analysis
# https://www.bioconductor.org/help/workflows/arrays/
# Page 70 https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf
# https://www.bioconductor.org/help/workflows/RNAseq123/
# ######################################################################
# 1.1 Reading in count-data
Taka_RNAseq <- read.csv("Takashima2014_3_RNAseq.csv", header = TRUE)
samplenames <-c("H9_R1","H9_R2","H9_R3",
                "H9_reset_R1","H9_reset_R2","H9_reset_R3")
Taka_RNAseq2 <- Taka_RNAseq[,c("Human.Gene.Symbol",samplenames)]
Taka_RNAseq2 <- Taka_RNAseq2[!duplicated(Taka_RNAseq2$Human.Gene.Symbol),]
rownames(Taka_RNAseq2)<- Taka_RNAseq2$Human.Gene.Symbol
Taka_RNAseq2<- Taka_RNAseq2[,-1]
head(Taka_RNAseq2)

# 1.2 Data pre-processing
Taka_RNAseq2 <- cpm(Taka_RNAseq2,log=TRUE)
keep.exprs <- rowSums(Taka_RNAseq2)>=4
Taka_RNAseq3 <- Taka_RNAseq2[keep.exprs,]
Taka_RNAseq3 <- Taka_RNAseq2
table(rowSums(Taka_RNAseq3)>=4)
head(Taka_RNAseq3)
allmean <- mean(colMeans(Taka_RNAseq3))  # faster version of apply(scaled.dat, 2, mean)
Taka_RNAseq3 <- scale(Taka_RNAseq3)
Taka_RNAseq3 <- Taka_RNAseq3 +allmean

# QC box plot
par(mfrow=c(1,2))
boxplot(Taka_RNAseq2)
boxplot(Taka_RNAseq3)

# QC density plot
nsamples <- ncol(Taka_RNAseq2)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(Taka_RNAseq2[,1]), col=col[1], lwd=2, las=2, 
     main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
        den <- density(Taka_RNAseq2[,i])
        lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
plot(density(Taka_RNAseq3[,1]), col=col[1], lwd=2, las=2, 
     main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=0, lty=3)
for (i in 2:nsamples){
        den <- density(Taka_RNAseq3[,i])
        lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")

# QC clustering of samples
par(mfrow=c(1,2))
col.group <- group <- as.factor(c("conventional", "conventional", "conventional",
                                  "reset", "reset", "reset"))
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(Taka_RNAseq2, labels=group, col=col.group)
title(main="A. Sample groups")
plotMDS(Taka_RNAseq3, labels=group, col=col.group)
title(main="B. Sequencing lanes")


# 1.3 Differential expression analysis
design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design
contr.matrix <- makeContrasts(
        ResetvsConv = reset-conventional, 
        levels = colnames(design))
contr.matrix
# Removing heteroscedascity from data
par(mfrow=c(1,2))
v <- voom(Taka_RNAseq3, design, plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean variance trend")

# Examining the number of DE genes
summary(decideTests(efit))
tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

# Examining individual DE genes from top to bottom
Reset.vs.Conv <- topTreat(tfit, coef=1, n=Inf)
head(Reset.vs.Conv)
write.csv(Reset.vs.Conv,"Takashima2014_3_RNAseq_DE.csv")
#Useful graphical representations of differential expression results
par(mfrow=c(1,1))
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], 
       xlim=c(-8,13))
#glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1],
#        side.main="ENTREZID", groups=group, launch=FALSE)
Reset.vs.Conv.topgenes <- rownames(Reset.vs.Conv)[1:100]
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(v$E[Reset.vs.Conv.topgenes,], scale="row",
          labCol=group, 
          col=mycol, trace="none", density.info="none", 
          margin=c(8,6), lhei=c(2,10), dendrogram="column")

# ######################################################################
# 
#  2 Prepare GSEA files
#  http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
#
# ######################################################################

# Prepare gmx file

libraries("xlsx","stringr")
Theunissen2014 <- read.xlsx2("mmc2.xlsx", 1)
gene<- str_split_fixed(Theunissen2014$Gene, "///",2)
Theunissen2014$Gene <- gene[,1]
write.csv(Theunissen2014, "Theunissen2014.csv")
#Theunissen_Naive.vs.Prime <- positive logFC (>3.0) & P.Value <0.05
#Theunissen_Prime.vs.Naive <- negative logFC (<-2.5) & P.Value <0.05
# Prepare gmx file
#Theunissen_Naive.vs.Prime <- as.numeric(as.character(Theunissen2014$Log2_FC_Naive.vs..hESM.))>3
#Theunissen_Naive.vs.Prime <- unique(Theunissen2014$Gene[Theunissen_Naive.vs.Prime])
#write.xlsx(Theunissen2014, "Theunissen2014.xlsx") # too slow

# Takashima2014_3_RNAseq_DE.csv
# Takashima_Naive.vs.Prime <- positive logFC (>1.3) & P.Value <0.05 & adj.P.value < 0.1
# Takashima_Prime.vs.Naive <- Negative logFC & P.Value <0.05

# Combine two gene sets

# ######################################################################
# 
#  3 Prepare heatmap filter with geneset
#  YAP Induces Human Naive Pluripotency Figuer Figure 3
#
# ######################################################################

list_files <- read.csv("list_files.csv",header = T)
list_files
# 3.1 read all csv files into a list, filter, and then merge======
list_matrix <- list()
matrix <- read.csv(as.character(list_files$dataset[1]))[,1:2] #read only frist two columns
for(i in c(1:9)){ 
        list_matrix[[i]] <- read.csv(as.character(list_files$dataset[i]))
        matrix <- merge(matrix,list_matrix[[i]],by ="X") }
dim(matrix)
head(matrix[,1:4])
rownames(matrix) <-matrix[,1]
matrix <- matrix[,-c(1:2)] #remove frist two columns

# 3.2 filter with std=============

matrix.filter <- matrix[order(apply(matrix,1,mad),decreasing = T)[1:1000],]
dim(matrix.filter)
# 3.2 filter with geneset----------
NaivevsPrime_geneset <- read.csv2("NaivevsPrime_geneset.csv",header = FALSE)
NaivevsPrime <- rownames(matrix) %in% as.character(NaivevsPrime_geneset$V1)
matrix <- matrix[NaivevsPrime,]
dim(matrix)

# 3.3 subset and group for both prim and naive======
sample_names <- read.csv("sample_names.csv",row.names = 1,header = T)
# naive
samples <- sample_names[sample_names$Label !="",
                            c("cell_names",
                              "authors",
                              "stages",
                              "file_names",
                              "Label")]
matrix.label <- matrix[,as.numeric(rownames(samples))]
colnames(matrix.label) <- samples$Label

# 3.4 average and normalize========
# normalize
matrix.scale <- scale(matrix.label)
# check that we get mean of 0 and sd of 1
allmean <- mean(colMeans(matrix.label))  # faster version of apply(scaled.dat, 2, mean)
matrix.scale <- matrix.scale +allmean 
par(mfrow= c(1,2))
boxplot(matrix.label)
boxplot(matrix.scale)

#====3.5 Calculate Spearman correlation (required)=================================
c <- cor(matrix.scale, method="spearman") # or naive_matrix
diag(c) <-NA

pheatmap(c,cex=.9,
         cluster_rows= T,
         cluster_cols = T,
         fontsize_row = 8,
         fontsize_col = 10,
         fontsize =25,
         main ="Hierarchical clustering")

#==3.6 average pheatmap===========
# average
Label <- droplevels(unique(samples$Label))
matrix.average <- data.frame(matrix(0, #create 0 matrix
                                    nrow = nrow(matrix.scale), #number of gene
                                    ncol = length(Label))) #number of cell types
for(row in 1:nrow(matrix.scale)){        
        for(col in 1:length(Label)){
                matrix.average[row,col] <-mean(matrix.scale[row,
                                                            colnames(matrix.scale)==Label[col]])}
}
rownames(matrix.average) <-rownames(matrix.scale)
colnames(matrix.average) <-Label

# adjust
colnames(matrix.average)
        colnames(matrix.average[,-c(6:7,9:11)])
matrix.adjust <- matrix.average[,-c(6:7,9:11)]
c <- cor(matrix.adjust, method="pearson") # or naive_matrix
diag(c) <-NA
## Row clustering (adjust here distance/linkage methods to what you need!)
hr <- hclust(as.dist(1-cor(matrix.adjust, method="spearman")), method="complete")

par(oma=c(2,2,2,2))
heatmap.2(c,
#          dendrogram="both",
#          Rowv=as.dendrogram(hr),
#          Colv=as.dendrogram(hr),
#          hclustfun = hclust,
          xlab = NULL,
          ylab = NULL,
          cexCol = 1.4,
          cexRow = 1.4,
          key=TRUE,
          na.color="#FFFFB2",
          keysize=0.8,
          key.xlab = "Distance metric",
          trace="none",
          density.info=c("none"),
          margins = c(14,14),
          col=colorRampPalette(brewer.pal(n = 9,
                                              name ="YlOrBr"))(100)
)
# KNN===================
colnames(matrix.scale) <- make.unique(colnames(matrix.scale), sep = ".")

kclust = kmeans(t(matrix.scale),centers=6)
#kclust$cluster <- names(matrix.label)
Cell.Type <- names(matrix.label)
d = dist(t(matrix.scale), method = "euclidean") 
fit = cmdscale(d,eig=TRUE, k=2) # k is the number of dim

# with legend
p = ggplot(data.frame(t(matrix.scale)),
           aes(fit$points[,1], fit$points[,2],
               color =  Cell.Type))+
        geom_point(size=6)
p= p + ggtitle("KNN cluster filter with gene sets")#ggplot title
p <- p + theme(axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.text=element_text(size=16),
               plot.title = element_text(hjust = 0.5,size=35), #ggplot title size
               legend.text = element_text(size = 14, colour = "black"),
             legend.title = element_text(size = 18, colour = "black"),
             legend.key.size = unit(1.5, "cm")) #title in middle
p

# without legend
p = ggplot(data.frame(t(matrix.scale)),
           aes(fit$points[,1], fit$points[,2],
               color =  Cell.Type,
               label = Cell.Type))+
        geom_point(size=8)
p =  p+ geom_text()
p= p + ggtitle("KNN cluster with 1000 genes")#ggplot title
p <- p + theme(axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.text=element_text(size=16),
               plot.title = element_text(hjust = 0.5,size=35), #ggplot title size
               legend.position="none") #title in middle
p
# ######################################################################
# 
#  4 WCGNA
#
# ######################################################################
# run Rscript in server with 5. forYangHu_WGCNAcode.R
# change input csv, repalce author name every time.


