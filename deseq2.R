d.gene.counts <- readr::read_delim(
    'BGI_F21FTSECKF0249_PEOaobdT/Expression/gene_expression.xls',
    delim = '\t'
) %>%
    select(
        gene_id, gene_symbol, starts_with('read')
    ) %>%
    mutate(
        across(starts_with('read'), floor)
    )

colData <- data.frame(
    row.names = colnames(d.gene.counts)[c(-1:-2, -11)],
    group = rep(c('Olig2OPC', 'VectorOPC'), c(6, 5)) %>% 
        factor(levels = c('VectorOPC', 'Olig2OPC'))
)

dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = d.gene.counts[, c(-1:-2, -14)],
    colData = colData,
    design = ~group
) %>% 
    DESeq2::DESeq()

dds.res <- DESeq2::results(
    dds, independentFiltering = F, contrast = c('group', 'Olig2OPC', 'VectorOPC'))

res <- dds.res@listData %>%
    data.frame(
        gene_id = d.gene.counts$gene_id,
        gene_symbol = d.gene.counts$gene_symbol,
        .
    ) %>%
    filter(
        !is.na(log2FoldChange) &
            padj < 0.05
    ) %>%
    arrange(
        as.integer(gene_id)
    ) 

write.csv(res, 'output/210/DESeq2_diff_gene.csv', row.names = F)