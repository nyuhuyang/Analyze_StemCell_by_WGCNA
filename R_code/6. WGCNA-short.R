# Weighted gene co-expression network analysis with RNAseq and microarry data
# Andreas Mock
# Cancer Research UK Cambridge Institute, University of Cambridge

# by 2017-04-24

# Contents

# 1 Assembly and preprocessing of TCGA RNAseq data
# 2 Construction of co-expression network
# 3 Identification of co-expression modules
# 4 Relation of co-expression modules to sample traits
# 5 Exploration of individual genes within co-expression module
# 6 Session information
# 7 References

########################################################################
#
#  0 check and install all cran and bioconductor packages if necessary
# 
# ######################################################################

list.of.cran.packages<- c("easypackages",
                          "ape")
new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(easypackages)

source("https://bioconductor.org/biocLite.R")
list.of.bio.packages <- c("WGCNA","limma")
new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite(new.packages)

libraries(list.of.bio.packages,list.of.cran.packages)

# Display the current working directory=======================
getwd()

if (Sys.info()[['sysname']]=="Darwin"){
        setwd("/Users/yah2014/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
if (Sys.info()[['sysname']]=="Windows"){
        setwd("C:/Users/User/Dropbox/Public/Olivier/R/Danwei_StemCell/dataset");getwd();list.files()}
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

########################################################################################
#
# 6.1 Assembly and preprocessing of RNAseq data
# https://bioconductor.org/packages/devel/bioc/vignettes/CVE/inst/doc/WGCNA_from_TCGA_RNAseq.html
#
########################################################################################

#6.1.1 Data input and cleaning (required)====

# Read in rnaseq_MicroArrayData
Danwei_Vassena = read.csv("rnaseq_MicroArrayData.csv",row.names = 1)
# Take a quick look at what is in the data set:
dim(Danwei_Vassena)
names(Danwei_Vassena)

# extract MicroArrayData 

Vassena = as.data.frame(Danwei_Vassena[, -c(1:9)])
dim(Vassena)
names(Vassena)
head(rownames(Vassena))

#6.1.2 expression estimates with counts in less than 20% of cases==========
Vassena = Vassena[apply(Vassena,1,function(x) sum(x==0))<ncol(Vassena)*0.8,]
dim(Vassena)

#6.1.3 normalized RNAseq data with the voom methodology
Vassena_voom = voom(Vassena)$E

#6.1.4 extract the top 5000 most variant genes for WGCNA studies.
#transpose matrix to correlate genes in the following
WGCNA_matrix = t(Vassena_voom[order(apply(Vassena_voom,1,mad), decreasing = T)[1:5000],])


########################################################################################
#
# 6.2  Construction of co-expression network
#
########################################################################################
#6.2.1 Employ the biweight midcorrelation between the profiles of node ii and jj
#similarity measure between gene profiles: biweight midcorrelation
s = abs(bicor(WGCNA_matrix))

#a weighted adjacency is defined by raising 
#6.2.2 the co-expression similarity to a power ('soft' threshold):
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(WGCNA_matrix, powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
     type='n', main = paste('Scale independence'))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col='red')
abline(h=0.70,col='red')

#6.2.3 calculation of adjacency matrix
beta = 3
a = s^beta
#6.2.4 dissimilarity measure
w = 1-a
########################################################################################
#
# 6.3  Identification of co-expression modules
#
########################################################################################
# 6.3.1 To identify co-expression modules,
#genes are next clustered based on the dissimilarity measure

#create gene tree by average linkage hierarchical clustering 
geneTree = hclust(as.dist(w), method = 'average')

#module identification using dynamic tree cut algorithm
modules = cutreeDynamic(dendro = geneTree, distM = w, deepSplit = 4, pamRespectsDendro = FALSE,
                        minClusterSize = 30)
#assign module colours
module.colours = labels2colors(modules)

#plot the dendrogram and corresponding colour bars underneath
plotDendroAndColors(geneTree, module.colours, 'Module colours', dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main='')


#6.3.2 the first principal component of its expression matrix (eigengene)
#calculate eigengenes
MEs = moduleEigengenes(WGCNA_matrix, colors = module.colours, excludeGrey = FALSE)$eigengenes

#calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

#cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = 'average');

#plot the result with phytools package
par(mar=c(2,2,2,2))
plot.phylo(as.phylo(METree),type = 'fan',show.tip.label = FALSE, main='')
tiplabels(frame = 'circle',col='black', text=rep('',length(unique(modules))), bg = levels(as.factor(module.colours)))

########################################################################################
#
# 6.4 Relation of co-expression modules to sample traits
#
########################################################################################

#load clinical metadata. Make sure that patient barcodes are in the same format 
#create second expression matrix for which the detailed clinical data is available 
WGCNA_matrix2 = WGCNA_matrix[match(clinical$Name, rownames(WGCNA_matrix)),]

#CAVE: 1 sample of detailed clinical metadata is not in downloaded data (TCGA-GN-A269-01')
not.available = which(is.na(rownames(WGCNA_matrix2))==TRUE)
WGCNA_matrix2 = WGCNA_matrix2[-not.available,]
str(WGCNA_matrix2)

#hence it needs to be removed from clinical table for further analysis
clinical = clinical[-not.available,]

#grouping in high and low lymphocyte score (lscore)
lscore = as.numeric(clinical$LYMPHOCYTE.SCORE)
lscore[lscore<3] = 0
lscore[lscore>0] = 1

#calculate gene significance measure for lymphocyte score (lscore) - Welch's t-Test
GS_lscore = t(sapply(1:ncol(WGCNA_matrix2),function(x)c(t.test(WGCNA_matrix2[,x]~lscore,var.equal=F)$p.value,
                                                        t.test(WGCNA_matrix2[,x]~lscore,var.equal=F)$estimate[1],
                                                        t.test(WGCNA_matrix2[,x]~lscore,var.equal=F)$estimate[2])))
GS_lscore = cbind(GS.lscore, abs(GS_lscore[,2] - GS_lscore[,3]))
colnames(GS_lscore) = c('p_value','mean_high_lscore','mean_low_lscore',
                        'effect_size(high-low score)'); rownames(GS_lscore) = colnames(WGCNA_matrix2)
#reference genes = all 5000 top mad genes 
ref_genes = colnames(WGCNA_matrix2)

#create data frame for GO analysis
library(org.Hs.eg.db)
GO = toTable(org.Hs.egGO); SYMBOL = toTable(org.Hs.egSYMBOL)
GO_data_frame = data.frame(GO$go_id, GO$Evidence,SYMBOL$symbol[match(GO$gene_id,SYMBOL$gene_id)])

#create GOAllFrame object
library(AnnotationDbi)
GO_ALLFrame = GOAllFrame(GOFrame(GO_data_frame, organism = 'Homo sapiens'))

#create gene set
library(GSEABase)
gsc <- GeneSetCollection(GO_ALLFrame, setType = GOCollection())

#perform GO enrichment analysis and save results to list - this make take several minutes
library(GEOstats)
GSEAGO = vector('list',length(unique(modules)))
for(i in 0:(length(unique(modules))-1)){
        GSEAGO[[i+1]] = summary(hyperGTest(GSEAGOHyperGParams(name = 'Homo sapiens GO', 
                                                              geneSetCollection = gsc, geneIds = colnames(RNAseq)[modules==i], 
                                                              universeGeneIds = ref.genes, ontology = 'BP', pvalueCutoff = 0.05, 
                                                              conditional = FALSE, testDirection = 'over')))
        print(i)
}

cutoff_size = 100

GO_module_name = rep(NA,length(unique(modules)))
for (i in 1:length(unique(modules))){
        GO.module.name[i] = 
                GSEAGO[[i]][GSEAGO[[i]]$Size<cutoff_size,
                            ][which(GSEAGO[[i]][GSEAGO[[i]]$Size<cutoff_size,]$Count==max(GSEAGO[[i]][GSEAGO[[i]]$
                                                                                                              Size<cutoff.size,]$Count)),7]
}

GO.module.name[1] = 'module 0'
#calculate module significance
MS.lscore = as.data.frame(cbind(GS.lscore,modules))
MS.lscore$log_p_value = -log10(as.numeric(MS.lscore$p_value))
MS.lscore = ddply(MS.lscore, .(modules), summarize, mean(log_p_value), sd(log_p_value))
colnames(MS.lscore) = c('modules','pval','sd')
MS.lscore.bar = as.numeric(MS.lscore[,2])
MS.lscore.bar[MS.lscore.bar<(-log10(0.05))] = 0
names(MS.lscore.bar) = GO.module.name

METree.GO = METree
label.order = match(METree$labels,paste0('ME',labels2colors(0:(length(unique(modules))-1))))
METree.GO$labels = GO.module.name[label.order]
plotTree.wBars(as.phylo(METree.GO), MS.lscore.bar, tip.labels = TRUE, scale = 0.2)

########################################################################################
#
# 6.5 Exploration of individual genes within co-expression module
#
########################################################################################
#Calculate module membership
MM = abs(bicor(RNAseq, MEs))

#plot individual module of interest (MOI)
MOI = 3 #T cell differentiation co-expression module
plot(-log10(GS.lscore[modules==MOI,1]), MM[modules==MOI,MOI], pch=20,
     cex=(GS.lscore[modules==MOI,4]/max(GS.lscore[,4],na.rm=TRUE))*4,
     xlab='p-value (-log10) lymphocyte score', ylab='membership to module 3')
abline(v=-log10(0.05), lty=2, lwd=2)