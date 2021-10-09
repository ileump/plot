###############################
## Checkpoint
###############################
# dir.create(file.path('D:/R', "checkpoint"),
#            recursive = TRUE, showWarnings = FALSE)

getOption('repos')
local({
    r <- getOption("repos")
    r["CRAN"] <- "https://mran.microsoft.com/snapshot/2019-03-01"
    options(repos = r)
})

# getOption('repos')
# local({
#     r <- getOption("repos")
#     r["CRAN"] <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"
#     options(repos = r)
# })

options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

require('checkpoint')


if (require(checkpoint, quietly = T)) {
    options(install.packages.compile.from.source = "no")
    oldLibPaths <- .libPaths()

    options(checkpoint.mranUrl = 'http://mran.microsoft.com/snapshot')
    checkpoint(
        snapshotDate = '2019-03-01',
        R.version = '3.5.2',
        checkpointLocation = 'D:/R',
        forceInstall = F,
        scanForPackages = F,
        project = getwd()
    )
}

###############################
## Chinese language support
###############################
# sessionInfo()
Sys.setlocale('LC_ALL', 'Chinese')

###############################
## CRAN package
###############################

## DGCA:: differential correlation analysis

pkgs <- c('tidyverse', 'testthat', 'devtools',
          'gridExtra', 'ggrepel', 'ggsignif', 'ggthemes',
          'factoextra', 'DGCA',
          'dunn.test', 'multcomp', 'multcompView', 'psych', 'FSA', 'userfriendlyscience',
          'igraph', 'circlize', 'corrplot', 'cowplot', 'pheatmap',
          'tableone', 'broom', 'data.table',
          'ggpubr',
          'caret', 'pROC',
          'knitr', 'rmarkdown', 'bookdown', 'shiny', 'DT', 'plotly',
          'XLConnect', 'sas7bdat',
          'circlize',
          'BiocManager'
)

pkgs <- setdiff(pkgs, installed.packages())



install.packages(pkgs)

###############################
## Bioconductor package
###############################

bioc.pkgs <- c('mixOmics', 'ropls', 'WGCNA', 'MEGENA', 'DGCA')
pkgs <- setdiff(bioc.pkgs, installed.packages())


BiocManager::install(pkgs)

###############################
## MetaboAnalystR
###############################

metanr_packages <- function(){
  
  metr_pkgs <- c("Rserve", "ellipse", "scatterplot3d", "Cairo", 
                 "randomForest", "caTools", "e1071", "som", 
                 "impute", "pcaMethods", "RJSONIO", "ROCR", 
                 "globaltest", "GlobalAncova", "Rgraphviz", 
                 "preprocessCore", "genefilter", "pheatmap", 
                 "SSPA", "sva", "Rcpp", "pROC", "data.table", 
                 "limma", "car", "fitdistrplus", "lars", "Hmisc", 
                 "magrittr", "methods", "xtable", "pls", "caret", 
                 "lattice", "igraph", "gplots", "KEGGgraph", "reshape", 
                 "RColorBrewer", "tibble", "siggenes", "plotly", "xcms", 
                 "CAMERA", "fgsea", "MSnbase", "BiocParallel", "metap", 
                 "reshape2", "scales", 'spls')
  
  list_installed <- installed.packages()
  
  new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"]))
  
  if(length(new_pkgs)!=0){
    
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(new_pkgs, version = "3.8")
    print(c(new_pkgs, " packages added..."))
  }
  
  if((length(new_pkgs)<1)){
    print("No new packages added...")
  }
}

metanr_packages()

devtools::install_github(
    "xia-lab/MetaboAnalystR", build = TRUE, build_vignettes = FALSE,
    build_opts = c("--no-resave-data", "--no-manual", "--no-build-vignettes"))
# git clone https://github.com/xia-lab/MetaboAnalystR.git
# R CMD build MetaboAnalystR
# R CMD INSTALL MetaboAnalystR_*.tar.gz
