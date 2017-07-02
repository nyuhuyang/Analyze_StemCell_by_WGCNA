# =====================================================================================
# 
#  0 check and install all cran and bioconductor packages if necessary
# 
# =====================================================================================

list.of.cran.packages<- c("WGCNA","easypackages")
new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(easypackages)

source("https://bioconductor.org/biocLite.R")
list.of.bio.packages <- c("AnnotationDbi", "impute","GEOquery","oligo","org.Hs.eg.db","annotation",
                          "edgeR","preprocessCore")
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)


libraries(list.of.bio.packages,list.of.cran.packages)
memory.limit(size=65000)
# =====================================================================================
# 
#  1.a Load Expression data---stem cell stages
# 
# =====================================================================================
# ------load processed data----------------------------------------------------------
eList <- getGEO("GSE29397") #it takes some time
eList
names(eList)
eData <- eList[[1]] #extract expressionset
eData
names(pData(eData))#there is usually a lot of unnecessary stuff here
pData(eData)$title #enough cell information from title
sampleNames(eData) <-pData(eData)$title
boxplot(eData)
#normData <- rma(rawData) #  already normalized, doesn't work 
#normData
# -----------------------------------------------------------------



# ------Loading raw data (optional)-----------------------------------------------------------
#too slow, enough data from above processed data
eList2 <- getGEOSuppFiles("GSE29397")
eList2
list.files("GSE29397")
untar("GSE29397/GSE29397_RAW.tar", exdir = "GSE29397/CEL")
list.files("GSE29397/CEL")
celfiles <- list.files("GSE29397/CEL", full = TRUE)
rawData <- read.celfiles(celfiles)
rawData
exprs(rawData)
max(exprs(rawData))
# clean up the phenotype information
filename <- sampleNames(rawData)
pData(rawData)$filename <- filename
sampleNames <- sub(".*_", "", filename) #remove all contain before "_"
sampleNames <- sub(".CEL.gz$", "", sampleNames)
sampleNames(rawData) <- sampleNames
pData(rawData)$group <- ifelse(grepl("^OSA", sampleNames(rawData)),
                               "OSA", "Control")
pData(rawData)
#Normalization
boxplot(rawData) # stuck and didn't work
normData <- rma(rawData) # it works
normData
boxplot(normData)
par(mfrow=c(1,2))
boxplot(eData)
boxplot(normData)
# -----------------------------------------------------------------


#  ----------Loading annotation data (MINiML)---stem cell stages----------
annotation(eData)
platf <- getGEO(annotation(eData), AnnotGPL=TRUE)
anot <- data.frame(attr(dataTable(platf), "table"))

gbs1 <- pData(featureData(eData))$GB_LIST
gbs1 <- unlist(strsplit(as.character(gbs1), ","))
mapped <- select(org.Hs.eg.db, gbs1, c("ENTREZID","SYMBOL"), "ACCNUM")
mapped <- unique(mapped[,2:3])
#apply(mapped[,2:3], 2, function(x) sum(!is.na(x)))
#GPL6244       Homo sapiens   hugene10sttranscriptcluster
#biocLite("hugene10sttranscriptcluster.db")
colnames(mapped) <-c("PROBES","gene_ID")

# ------replace affymetrix ID to gene ID-----------------------------------------------------------
MicroArrayData <- exprs(eData)
PROBES <- as.numeric(rownames(MicroArrayData))
MicroArrayData<- cbind.data.frame(PROBES,MicroArrayData)

anno_MicroArrayData <- merge(mapped, MicroArrayData)
anno_MicroArrayData <- anno_MicroArrayData[order(anno_MicroArrayData[,"gene"]),]
anno_MicroArrayData <- anno_MicroArrayData[,-1]

