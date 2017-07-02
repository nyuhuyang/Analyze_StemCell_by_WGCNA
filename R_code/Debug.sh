#//
#//  Debug.sh
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

#i
grep("!Sample_title",MicroArrayData0[,1])    MicroArrayData0[,1])=="!Sample_title"



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

