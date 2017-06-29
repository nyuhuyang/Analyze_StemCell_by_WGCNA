//
//  Debug.sh
//  
//
//  Created by Yang Hu on 2017/6/20.
//
//

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

#remove emypt column