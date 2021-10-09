#-------------------------------------------
# checkpoint
# replicate results using specific R package versions.
#-------------------------------------------
# options(download.file.method = "auto")
# options(download.file.method = "auto", download.file.extra = c("--no-check-certificate"))

options(download.file.method="libcurl")

# dir.create(file.path('~/R', ".checkpoint"), recursive = TRUE, showWarnings = FALSE)
# options(install.packages.compile.from.source = "no")
oldLibPaths <- .libPaths()

# library(checkpoint)
# options(checkpoint.mranUrl = 'http://mran.microsoft.com/snapshot')
# checkpoint(
#     snapshotDate = '2019-09-01',
#     R.version = '3.6.1',
#     checkpointLocation = 'D:/R',
#     forceInstall = F,
#     scanForPackages = T,
#     project = getwd()
# )

options('repos')
.libPaths()
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
options()$BioC_mirror
# utils::setRepositories(ind=1)

## curl
## https://curl.haxx.se/download.html


cran_packages <- function(){
    
    metr_pkgs <- c('curl', 'httr', 'devtools', 'testthat', 'kableExtra',
                   'bookdown', 'knitr', 'pander', 'tableone', 'broom',
                   'ggplot2', 'ggdendro', 'ggalluvial', 'tidyverse', 
                   'XLConnect', 'data.table', 'corrplot',
                   'pheatmap', 'factoextra', 'igraph', 'plotly', 'psych',
                   'DT', 'multcompView', 'bootstrap', 'Rtsne', 'VennDiagram',
                   'extrafont', 'animation', 'magick', 'pdftools',
                   'caret', 'mlbench', 'MLmetrics', 'skimr', 'RANN', 'randomForest', 
                   'fastAdaboost', 'gbm', 'xgboost', 'caretEnsemble', 'C50', 'earth',
                   'kernlab', 'rpartScore', 'ordinalNet', 'pryr', 'egg', 'survminer')
    # WGCNA, DGCA
    list_installed <- installed.packages()
    
    new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
    
    if(length(new_pkgs)!=0){
        
        lapply(new_pkgs, function(p) {
            install.packages(p, dependencies = T)
        })
        # install.packages(new_pkgs, dependencies = TRUE)
        print(c(new_pkgs, " packages added..."))
    }
    
    if((length(new_pkgs)<1)){
        print("No new packages added...")
    }
}

cran_packages()

## bioconductor
install.packages("BiocManager", dependencies = T)
BiocManager::version()
BiocManager::valid()
BiocManager::install(version = '3.11')

# source("http://bioconductor.org/biocLite.R")

metanr_packages <- function(){
    
    metr_pkgs <- c("Rserve", "ellipse", "scatterplot3d", "Cairo", 
                   "randomForest", "caTools", "e1071", "som", "impute", 
                   "pcaMethods", "RJSONIO", "ROCR", "globaltest", 
                   "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", 
                   "pheatmap", "SSPA", "sva", "Rcpp", "pROC", "data.table", 
                   "limma", "car", "fitdistrplus", "lars", "Hmisc", "magrittr", 
                   "methods", "xtable", "pls", "caret", "lattice", "igraph", 
                   "gplots", "KEGGgraph", "reshape", "RColorBrewer", "tibble", 
                   "siggenes", "plotly", "ropls", 'BridgeDbR')
    
    list_installed <- installed.packages()
    
    new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
    
    if(length(new_pkgs)!=0){
        
        # source("http://bioconductor.org/biocLite.R")
        # biocLite(new_pkgs, dependencies = TRUE, ask = FALSE)
        BiocManager::install(new_pkgs)
        print(c(new_pkgs, " packages added..."))
    }
    
    if((length(new_pkgs)<1)){
        print("No new packages added...")
    }
}

metanr_packages()

bioc_packages <- function(){
    
    metr_pkgs <- c('WGCNA', 'DGCA', 'mixOmics', 'ropls', 'ComplexHeatmap')
    
    list_installed <- installed.packages()
    
    new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
    
    if(length(new_pkgs)!=0){
        
        # source("http://bioconductor.org/biocLite.R")
        BiocManager::install(new_pkgs)
        print(c(new_pkgs, " packages added..."))
    }
    
    if((length(new_pkgs)<1)){
        print("No new packages added...")
    }
}

bioc_packages()

## metaboanalystr
metanr_packages <- function(){
    
    metr_pkgs <- c("Rserve", "ellipse", "scatterplot3d", "Cairo", 
                   "randomForest", "caTools", "e1071", "som", "impute", 
                   "pcaMethods", "RJSONIO", "ROCR", "globaltest", 
                   "GlobalAncova", "Rgraphviz", "preprocessCore", 
                   "genefilter", "pheatmap", "SSPA", "sva", "Rcpp", 
                   "pROC", "data.table", "limma", "car", "fitdistrplus", 
                   "lars", "Hmisc", "magrittr", "methods", "xtable", "pls", 
                   "caret", "lattice", "igraph", "gplots", "KEGGgraph", 
                   "reshape", "RColorBrewer", "tibble", "siggenes", "plotly", 
                   "xcms", "CAMERA", "fgsea", "MSnbase", "BiocParallel", 
                   "metap", "reshape2", "scales", "crmn")
    

    list_installed <- installed.packages()
    
    new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
    
    if(length(new_pkgs)!=0){
        
        if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
        BiocManager::install(new_pkgs)
        print(c(new_pkgs, " packages added..."))
    }
    
    if((length(new_pkgs)<1)){
        print("No new packages added...")
    }
}
metanr_packages()

BiocManager::install('ChemmineR')
devtools::install_git('F:/github/metabomapr/')