########################################################################
#
#  0 check and install all cran and bioconductor packages if necessary
# 
# ######################################################################
list.of.cran.packages<- c("easypackages","ggfortify")
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
#  7. PCA
# 
########################################################################################

# remove primed hESCs from Vassena and Xie dataset
sample_names <- read.csv("sample_names.csv",row.names = 1,header = T)
sample_names[c(25:27,46:48),"stages"] <-""
#rename primed hESCs

sample_names$authors <- ifelse(sample_names$authors=="Tachibana M et al",
                               "primed hESCs",as.character(sample_names$authors))
#===subset naive and prime====
#deselection
naive_primed_cells <- sample_names[sample_names$stage!="",c("cell_names","authors","stages","file_names")]
#deselected_authors <- c("Takashima et al","Hanna et al")
#naive_primed_cells <- naive_primed_cells[!(naive_primed_cells$authors %in% deselected_authors),]
#naive_primed
naive_primed_matrix <- matrix.scale[,as.numeric(rownames(naive_primed_cells))] # keep the order
colnames(naive_primed_matrix) <- naive_primed_cells$cell_names
# add label
rownames(naive_primed_cells)<-paste0(naive_primed_cells$cell_names,
                                     ".",
                                     rownames(naive_primed_cells))
#PCA plot without label

autoplot(prcomp(t(naive_primed_matrix)),
         data=naive_primed_cells,
         label = FALSE,
#         frame = TRUE, frame.type = 'norm',
         colour = 'authors',
         shape= 'stages',
#        label.size = 5,
         size = 8)+
        ggtitle("Principal Component Analysis for All Samples Before normalization")+#ggplot title
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5))+ #title in middle
        guides(colour = guide_legend(override.aes = list(size=10)), #larger legend diagram 
               shape = guide_legend(override.aes = list(size=10))) #larger legend diagram 

#PCA plot with label
autoplot(prcomp(t(naive_primed_matrix)),
         data=naive_primed_cells,
         label = TRUE,
         #         frame = TRUE, frame.type = 'norm',
         colour = 'authors',
         shape= 'stages',
         label.size = 5,
         size = 0)+
        ggtitle("Principal Component Analysis for All Samples")+#ggplot title
        theme(text = element_text(size=20),     #larger text including legend title							
              plot.title = element_text(hjust = 0.5))+ #title in middle
        guides(colour = guide_legend(override.aes = list(size=10)), #larger legend diagram 
               shape = guide_legend(override.aes = list(size=10))) #larger legend diagram 
