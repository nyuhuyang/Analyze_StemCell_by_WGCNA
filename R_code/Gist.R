#//
#//  Debug.R
#//  
#//
#//  Created by Yang Hu on 2017/6/20.
#//
#//

#include <stdio.h>

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
#Error in datk[c(startG:endG), ] = foreach(t = actualThreads, .combine = rbind) %dopar%  :
#number of items to replace is not a multiple of replacement length
#"The code below uses parallel computation where multiple cores are available. This works well when R is run from a terminal or from the Graphical User Interface (GUI) shipped with R itself, but at present it does not work with RStudio and possibly other third-party R environments." -WGCNA package tutorial




#Import data into R with an unknown number of columns?
#You need to use the fill argument in the read.table function
fill=TRUE
#There is nice function count.fields (see help) which counts number of column per row:
no_col <- max(count.fields("test", sep = "\t"))
data <- read.table("test",sep="\t",fill=TRUE,col.names=1:no_col)
data

#remove header which have "!"
#grepl returns a logical vector
grepl

#
grep("!Sample_title",MicroArrayData0[,1])
MicroArrayData0[,1]=="!Sample_title"



#the factor variable retains all of its original levels -- even when they do not exist in the new data frame.
droplevels

#Error in sort.list(y) : 'x' must be atomic for 'sort.list'
#Have you called 'sort' on a list?
unlist

#split column into multiple ones
str_split_fixed

#remove rows with empty value in column 1
X <- X[!X[,1]=="",]

#Increase momery
memory.limit(size=65000)

#---2017-07-01
#What is the difference between gsub and sub methods in r
#The g stands for global, as in replace globally (all):

sub('l', '*', "hello")  #=> "he*lo"
gsub('l', '*',"hello")  #=> "he**o"

#sub(".*") can remove undefined content
filename <- c("GSM726928_OOC1.CEL.gz","GSM726949_STC2_1.CEL.gz")
sampleNames <- sub(".*_", "", filename) #=> "OOC1.CEL.gz" "1.CEL.gz"
sampleNames <- sub(".CEL.gz$", "", sampleNames) #=> "OOC1" "1"


#Avoid rbind()/cbind() conversion from numeric to factor
cbind.data.frame

#merge unequal dataframes and replace missing rows with 0
merge(df1, df2, all = TRUE)

#adjust all log data to positive
logcounts0 <- cpm(expr0,log=TRUE)
min(logcounts0)
all((logcounts0+5.2)>0)

#What does %>% function mean in R?
#The dplyr package introduced the %.% operator to pass the left hand side
#as an argument of the function on the right hand side, similar to a *NIX pipe. |
df <- data.frame(
        x = sample(10, 10, rep = TRUE),
        y = sample(10, 10, rep = TRUE),
        z = sample(10, 10, rep = TRUE)
)
nrow(df)
nrow(distinct(df))
nrow(distinct(df, x, y))
nrow(distinct(df, x))

distinct(df, x)
distinct(df, y)

# Can choose to keep all other variables as well
distinct(df, x, .keep_all = TRUE)
distinct(df, y, .keep_all = TRUE)

# You can also use distinct on computed variables
distinct(df, diff = abs(x - y))

# The same behaviour applies for grouped data frames
# except that the grouping variables are always included
df <- tibble(
        g = c(1, 1, 2, 2),
        x = c(1, 1, 2, 1)
) %>% group_by(g)
df %>% distinct()
df %>% distinct(x)


#To remove the subtitle use the following:
plot(hc, xlab="", sub="")

#----------Task on 2017/07/05----------------
#try to remove genes that don't' change much in the RNAseq experiment
thresh <- myCPM > 0.5 #thresh number decide lower limits of std
keep <- rowSums(thresh) >= 9 # keep rowsum number decide individule circles

#make a heatmap out of the correlation matrix c:
#which microarray sample is each of Danwei's experiments closest to
pheatmap(cor(rnaseq_MicroArrayData, method="spearman"))
#pickup gene set with p value <0.05
#Do deseq, 

#find out linkID
#complete WGCNA


#Plot a data frame as a table in r
library(gridExtra)
title <-as.data.frame(eData$title)
dev.off()
grid.table(title)

# =====================================================================================
# 
#  4. DEseq  ----some values in assay are not integers
# 
# =====================================================================================
# read count matrix
rnaseq_MicroArrayData <-read.csv("rnaseq_MicroArrayData.csv", row.names = 1)

# create sample table for DESeq
sample.table <- data.frame( title = colnames(rnaseq_MicroArrayData),
                            group = c(rep("Danwei",9),
                                      as.character(eData$labelDescription))
)
dim(sample.table) 
## Creating a DESeqDataSet object
#create ddsHTSeq
ddsHTSeq <- DESeqDataSetFromMatrix( countData = rnaseq_MicroArrayData,
                                    colData = sample.table,
                                    design= ~ group)
ddsHTSeq


#Combining data frames with unequal number of columns [duplicate]
df1 = data.frame(A = c(ncol(rnaseq_MicroArrayData),2,1))
df2 = data.frame(B = c("#","Danwei","StemCell","StemCellLine"))
df3 = data.frame(C=c(rep("Danwei",9),as.character(eData$labelDescription)))
l <- merge(t(df1),t(df2),all=TRUE)
l <- merge(l,t(df3),all=TRUE)


#----------Task on 2017/07/07
#try to keep the genes that have higher std (top 25%) in affrymetric data
sd = apply(rnaseq_MicroArrayData[,10:ncol(rnaseq_MicroArrayData)], 1, sd)
sd = sd[order(sd, decreasing=TRUE)]
Top_Std_genes <-names(sd)[1:(length(sd)*0.25)]

# The las argument rotates the axis names
par(las=2)

#similarity of signatures
GSEA 4
#find out linkID
#complete WGCNA