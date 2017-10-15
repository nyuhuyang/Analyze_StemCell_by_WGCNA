
########################################################################
#
# Check, install, and load multiple packages from cran and bioconductor
# 
# ######################################################################



cran.libraries <- function(...){
        # For cran libraries
        list.of.cran.packages<- c("easypackages",...)
        new.packages <- list.of.cran.packages[!(list.of.cran.packages %in% installed.packages()[,"Package"])]
        if(length(new.packages)) install.packages(new.packages)
        library(easypackages)
        libraries(list.of.cran.packages)
        }



bio.libraries <- function(...){
        # For bioconductor libraries
        if(!("easypackages" %in% installed.packages()[,"Package"]))
                install.packages("easypackages")
        source("https://bioconductor.org/biocLite.R")
        list.of.bio.packages <- c(...)
        new.packages <- list.of.bio.packages[!(list.of.bio.packages %in% installed.packages()[,"Package"])]
        if(length(new.packages)) biocLite(new.packages)
        library(easypackages)
        libraries(list.of.bio.packages)
        }

# #####################################################################
# 
#  detect OS and set enviroment
#
# ####################################################################

initialize.project <- function(x){
        if (Sys.info()[['sysname']]=="Darwin"){
                WD <- paste0("/Users/yah2014/Dropbox/Public/Olivier/R/",x,"/datasets")
                setwd(WD);print(getwd());print(list.files())}
        if (Sys.info()[['sysname']]=="Windows"){
                WD <- paste0("C:/Users/User/Dropbox/Public/Olivier/R/",x,"/datasets")
                setwd(WD);print(getwd());print(list.files())}
        if (Sys.info()[['sysname']]=="Linux"){
                WD <- paste0("/pbtech_mounts/homes030/yah2014/R/",x,"/datasets")
                setwd(WD);print(getwd());print(list.files())}
        }
