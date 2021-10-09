## ht.file : hypothesis test
## data.file : group mean

library(magrittr)

data.file <- 'F:/Lipidall/projects/Reports/2019/2019-121-C-01 metabolomic flux/data/20200203/data_group.csv'
ht.file <- 'F:/Lipidall/projects/Reports/2019/2019-121-C-01 metabolomic flux/report/20200203/report/05_Hypothesis_test 假设检验的结果/hypothesis_test.csv'
var.file <- 'F:/Lipidall/projects/Reports/2019/2019-121-C-01 metabolomic flux/data/20200203/var.csv'

d.data <- read.csv(data.file, check.names = F, row.names = 1)
d.var <- read.csv(var.file, check.names = F)
d.ht <- read.csv(ht.file, check.names = F, row.names = 1) %>% dplyr::mutate(
    sig = ifelse(`non-parametric pvalue` < 0.05, '*', ' ')
)

pdf('F:/Lipidall/projects/Reports/2019/2019-121-C-01 metabolomic flux/report/20200203/heatmap.pdf',
    w=4, h=4)

purrr::walk(unique(d.var$Class), function(class.i) {
    var.ind <- which(d.var$Class == class.i)
    
    d.ind <- d.data[, var.ind, drop = F] %>% t
    rownames(d.ind) <- paste(
        d.ht[var.ind, 'sig'],
        stringr::str_extract(rownames(d.ind), 'M\\+[0-9]+$')
    )
    
    pheatmap::pheatmap(d.ind, cellwidth = 20, cellheight = 10, fontsize = 10,
                       scale = 'none', main = class.i, angle_col = 0,
                       display_numbers = T,
                       legend = F,
                       cluster_rows = FALSE, cluster_cols = FALSE) %>%
        print
})

dev.off()