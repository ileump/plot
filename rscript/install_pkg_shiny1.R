#-------------------------------------------
# checkpoint
# replicate results using specific R package versions.
#-------------------------------------------
options(download.file.method = "auto")
# options(download.file.method = "auto", download.file.extra = c("--no-check-certificate"))


dir.create(file.path('~/R', ".checkpoint"), recursive = TRUE, showWarnings = FALSE)
options(install.packages.compile.from.source = "no")
oldLibPaths <- .libPaths()

# library(checkpoint)
# options(checkpoint.mranUrl = 'http://mran.microsoft.com/snapshot')
# checkpoint(
#     snapshotDate = '2019-01-22',
#     R.version = '3.5.1',
#     checkpointLocation = '~/R',
#     forceInstall = FALSE,
#     scanForPackages = FALSE,
#     project = getwd()
# )

options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
options()$BioC_mirror
# utils::setRepositories(ind=1)

cran_packages <- function(){
    
    metr_pkgs <- c('curl', 'httr', 'devtools', 'testthat', 'kableExtra',
                   'bookdown', 'knitr', 'pander', 'tableone', 'broom',
                   'ggplot2', 'ggdendro', 'ggalluvial', 'tidyverse', 
                   'XLConnect', 'data.table', 
                   'userfriendlyscience', 'FSA',
                   'pheatmap', 'factoextra', 'igraph', 'plotly', 'psych',
                   'DT', 'multcompView', 'bootstrap', 'Rtsne', 'VennDiagram',
                   'shinyjs')
    # WGCNA, DGCA
    list_installed <- installed.packages()
    
    new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
    
    if(length(new_pkgs)!=0){
        
        install.packages(new_pkgs, dependencies = TRUE)
        print(c(new_pkgs, " packages added..."))
    }
    
    if((length(new_pkgs)<1)){
        print("No new packages added...")
    }
}

cran_packages()

## bioconductor
# install.packages("BiocManager", dependencies = T)
# BiocManager::version()
# BiocManager::valid()
# BiocManager::install(version = '3.8')

source("http://bioconductor.org/biocLite.R")

metanr_packages <- function(){
    
    metr_pkgs <- c("Rserve", "ellipse", "scatterplot3d", "Cairo", 
                   "randomForest", "caTools", "e1071", "som", "impute", 
                   "pcaMethods", "RJSONIO", "ROCR", "globaltest", 
                   "GlobalAncova", "Rgraphviz", "preprocessCore", "genefilter", 
                   "pheatmap", "SSPA", "sva", "Rcpp", "pROC", "data.table", 
                   "limma", "car", "fitdistrplus", "lars", "Hmisc", "magrittr", 
                   "methods", "xtable", "pls", "caret", "lattice", "igraph", 
                   "gplots", "KEGGgraph", "reshape", "RColorBrewer", "tibble", 
                   "siggenes", "plotly", "ropls")
    
    list_installed <- installed.packages()
    
    new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
    
    if(length(new_pkgs)!=0){
        
        source("http://bioconductor.org/biocLite.R")
        biocLite(new_pkgs, dependencies = TRUE, ask = FALSE)
        print(c(new_pkgs, " packages added..."))
    }
    
    if((length(new_pkgs)<1)){
        print("No new packages added...")
    }
}

metanr_packages()

bioc_packages <- function(){
    
    metr_pkgs <- c('WGCNA', 'DGCA')
    
    list_installed <- installed.packages()
    
    new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
    
    if(length(new_pkgs)!=0){
        
        source("http://bioconductor.org/biocLite.R")
        biocLite(new_pkgs, dependencies = TRUE, ask = FALSE)
        print(c(new_pkgs, " packages added..."))
    }
    
    if((length(new_pkgs)<1)){
        print("No new packages added...")
    }
}

bioc_packages()

