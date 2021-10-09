library(multiGSEA)

library(metaboliteIDmapping)

file.RData <- 'RData/210-limma_metabolome.RData'

if (file.exists(file.RData)) {
  load(file.RData)
  
  fit.single <- limma.fit.metabolome
}

testthat::expect_equal(
  colnames(data.metabolome$log2),
  data.metabolome$var$name
)

# 'kegg', 'panther', 'pathbank', 'pharmgkb', 'reactome', 'smpdb', 'wikipathways'

databases <- c('kegg')
organism <- 'hsapiens'

pathways.graphite <- lapply(databases, function(x) {
  graphite::pathways(organism, x)
})
names(pathways.graphite) <- databases

d.pathway <- lapply(pathways.graphite, function(pathway.db) {
  lapply(pathway.db, function(pathway) {
    c(pathway@id, pathway@title, pathway@database)
  }) %>% 
  do.call(what = rbind) %>%
    as.data.frame %>% 
    magrittr::set_colnames(c('pathway_id', 'pathway_title', 'database')) %>%
    dplyr::mutate(
      name = sprintf('(%s) %s', database, pathway_title)
    )
  }) %>%
  do.call(what = rbind) %>%
  magrittr::set_rownames(NULL)


d.metabolites <- purrr::map_dfr(pathways.graphite, function(pathway.db) {
  lapply(pathway.db, function(pathway) {
    native_id <- graphite::nodes(pathway, which = 'metabolites')
    
    if (length(native_id) > 0) {
      data.frame(
        database = pathway@database,
        pathway_id = pathway@id,
        pathway_title = pathway@title, 
        native_id = native_id
      )
    } else {
      NULL
    }
  }) %>% 
    do.call(what = rbind)
}) %>% 
  magrittr::set_rownames(NULL) %>% 
  dplyr::mutate(
    KEGG = stringi::stri_extract_first_regex(native_id, '(?<=KEGGCOMP:)C[0-9]{5}$')
  ) %>%
  dplyr::left_join(
    tidyr::drop_na(metabolitesMapping, KEGG), by = 'KEGG'
  )

data.metabolome <- read.csv('var.csv', check.names = F, row.names = 1)

pathways <- getMultiOmicsFeatures(
  dbs = databases, layer = 'metabolome', 
  # returnTranscriptome = 'SYMBOL', 
  # returnProteome = 'SYMBOL',
  returnMetabolome = 'HMDB'
)


res.fgsea <- purrr::map(
  colnames(fit2.single$contrasts) %>% purrr::set_names(), 
  ~ topTable(fit2.single, coef = ., number = Inf, adjust.method = 'BH')[
    rownames(data.metabolome), ]
    # colnames(data.metabolome$log2), ]
)[
  c('T2.N - T1.N', 
    'T3.N - T2.N',
    'T4.N - T3.N',
    'T5.N - T4.N',
    'T6.N - T5.N',
    'T7.N - T6.N',
    'T8.N - T7.N')
] %>%
  purrr::imap(function(d.limma, contrast.name) {
    testthat::expect_equal(
      rownames(d.limma), 
      rownames(data.metabolome)
    )
    
    t.stat <- d.limma %>% 
      dplyr::mutate(
        HMDB = data.metabolome$`HMDB ID`
      ) %>%
      tidyr::drop_na(HMDB) %>% 
      dplyr::pull(t, name = HMDB)
    
    set.seed(1)
    fgsea::fgseaMultilevel(pathways$metabolome, t.stat, eps = 0) %>%
      dplyr::mutate(
        ## leadingEdge is a list
        leadingEdge = purrr::map_chr(leadingEdge, ~ paste0(., collapse = ','))
      )
  })

output.dir.1 <- '汇报/'

dir.name <- file.path(output.dir.1, 'fgsea_matabolome')

if (!dir.exists(dir.name)) {
  dir.create(dir.name, recursive = T)
}


purrr::iwalk(res.fgsea, function(x, contrast.name) {
  file.csv <- file.path(dir.name, paste0('fgsea_', contrast.name, '.csv'))
  
  if (!file.exists(file.csv)) {
    write.csv(x, file.csv)
  }
})

file.pdf <- file.path(dir.name, 'fgsea_KEGG_summary.pdf')

if (!file.exists(file.pdf)) {
  cairo_pdf(file.pdf, onefile = T, width = 16, height = 12)
  
  purrr::walk(c('pval', 'padj'), function(p) {
    p_exp <- parse_expr(p)
    
    d1 <- purrr::imap_dfr(res.fgsea, function(x, y) {
      x %>% 
        dplyr::select(pathway, !!p_exp, size) %>%
        dplyr::mutate(
          contrast = y, 
          pathway = forcats::fct_reorder2(pathway, !!p_exp, -size)
        )
    })
    
    p1 <- ggplot(d1, aes(x = pathway, y = -log10(!!p_exp))) +
      geom_point(aes(size = size, color = !!p_exp < 0.05)) + 
      scale_color_manual(values = c('gray', 'red')) +
      geom_hline(yintercept = -log10(0.05)) +
      labs(x = '') +
      coord_flip() + 
      facet_wrap(~ contrast, nrow = 1) +
      theme_bw()
    
    print(p1)
  })
  
  dev.off()
}

















































































