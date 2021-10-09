#-------------------------------------------
# checkpoint
# replicate results using specific R package versions.
#-------------------------------------------
# dir.create(file.path('~/R', ".checkpoint"), recursive = TRUE, showWarnings = FALSE)
# options(install.packages.compile.from.source = "no")
# oldLibPaths <- .libPaths()

# library(checkpoint)
# options(checkpoint.mranUrl = 'http://mran.microsoft.com/snapshot')
# checkpoint(
#     snapshotDate = '2019-03-01',
#     R.version = '3.5.2',
#     checkpointLocation = 'D:/R',
#     forceInstall = FALSE,
#     scanForPackages = FALSE,
#     project = getwd()
# )


knitr::opts_chunk$set(
    echo = FALSE,
    message = FALSE,
    include = FALSE,
    warning = FALSE,
    purl = FALSE,
    dpi = 300
)

options(java.parameters = '-Xmx4g')
options(stringsAsFactors = FALSE)

# library(data.table, quietly = TRUE)
# library(DT)

# color blind friendly color pallete
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
               "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

set.seed(123)