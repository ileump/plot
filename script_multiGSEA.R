library(multiGSEA)
library(org.Hs.eg.db)
library(metaboliteIDmapping)
library(purrr)
library(rlang)

layers <- c('transcriptome', 'metabolome')
data(list = layers, package = 'multiGSEA')




odata <- initOmicsDataStructure(layer = layers)
odata <- layers %>%
  set_names() %>%
  map(function(x) {
    ome <- get(x)
    features <- rankFeatures(ome$logFC, ome$pValue)
    names(features) <- ome[, 1] %>% unlist %>% 
      gsub(pattern = 'HMDB', replacement = 'HMDB00')
    sort(features)
  })

metabo <- read.csv('/Users/leump/Desktop/metabo_blue.csv', check.names = F, row.names = 1)
rna <- read.csv('/Users/leump/Desktop/rna_blue.csv', check.names = F, row.names = 1)

metabo <- sort(metabo$`Fold: KO/WT`)
rna <- sort(rna$Fold..KO.WT)

l1 <- list(transcriptome = rna, 
           metabolome = metabo)

str(l1)





databases <- c('kegg')
pathways <- getMultiOmicsFeatures(
  dbs = databases,
  layer = layers,
  returnTranscriptome = 'SYMBOL',
  # returnProteome = 'SYMBOL',
  returnMetabolome = 'HMDB'
)

es <- multiGSEA(pathways, odata)

df <- extractPvalues(enrichmentScores = es,
                     pathwayNames = names(pathways$transcriptome))

## pathway column is not in df
df$combined_pval <- combinePvalues(df)
df$combined_padj <- p.adjust(df$combined_pval, method = 'BH')
df$pathway <- names(pathways$transcriptome)

head(df[, c('pathway', 'combined_pval', 'combined_padj')])

write.csv(df, '/Users/leump/Desktop/multiGSEA.csv', row.names = T)












