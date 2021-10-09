pheatmapByClass <- function(x, v, s, test, prefix) {
    ## x :: data.frame, data
    ## v :: data.frame, getAttrs(colnames(x))
    ## s :: sample annotation
    for (i in 1:length(unique(v[, 'class']))) {
        ind <- which(v[, 'class'] == unique(v[, 'class'])[i])
        if (length(ind) < 3) next
        d <- t(x[, ind])
        
        ann_row <- data.frame(ifelse(test[, 1] < 0.05, '<0.05', '>0.05')) %>% 
            `rownames<-`(rownames(test)) %>%
            `colnames<-`(colnames(test)[1])
        
        ann_cols <- vector('list', ncol(ann_row)) %>% `names<-`(colnames(ann_row))
        for (j in 1:ncol(ann_row))
            ann_cols[[j]] <- c('<0.05' = 'red', '>0.05' = 'white')
        
        sub_chunk <- paste0(
            '\n```{r ', prefix, i, ', echo=F, include=T,', "results='asis',",
            'fig.width=', 5.5, 
            ', fig.height=', 3+length(ind)/8,
            '}\n',
            "


g <- pheatmap::pheatmap(
    d,
    cluster_cols = F,
    scale = 'row',
    annotation_col = s,
    annotation_row = ann_row,
    annotation_colors = ann_cols,
    silent = T
)

print(g)
## to avoid adding multiple plots to the same device
invisible(dev.off())
```\n
")

cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = T))
    }
}
