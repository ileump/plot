setwd(fs::path_norm(r"(F:\Lipidall\Shui G\2020\ACA\R)"))


##---------------------------------
## Parse json
##---------------------------------
library(tidyjson)
library(magrittr)

doc <- jsonlite::fromJSON('br08002.json')

getNodeData <- function(x, parent = '') {
  if (length(x) == 1 && names(x) == 'name') {
    data.frame(
      parent = parent,
      stringr::str_split(x$name, ' +', n = 2, simplify = T)
    )
  } else if (length(x) == 2 && names(x) == c('name', 'children')) {
    # root
    if (length(x$name) != length(x$children)) {
      return(getNodeData(x$children))
    } else {
      lapply(1:length(x$name), function(i) {
        getNodeData(x$children[[i]], 
                    parent = paste0(ifelse(parent == '', '', paste0(parent, '; ')), x$name[i]))
      }) %>% do.call(rbind, .)
    }
  }
}

doc.data <- getNodeData(doc)

colnames(doc.data) <- c('parent', 'compound_id', 'compound_name')

doc.data <- doc.data[, c('compound_id', 'compound_name', 'parent')]

parents <- doc.data$parent %>% stringr::str_split('; ', simplify = T) %>%
  `colnames<-`(paste('parent level', 1:ncol(.)))

doc.data <- cbind(doc.data, parents)
write.csv(doc.data, 'KEGG_lipids.csv')

##---------------------------------
## KEGG human pathway
##---------------------------------
library(KEGGREST, quietly = TRUE)
library(tidyverse, quietly = TRUE)

## retrieve all human pathway from KEGG
## N = 337, 20200816
hsaList <- keggList("pathway", "hsa")
IDList <- names(hsaList) %>% map_chr(pathwayID)
hsaPathway <- data.frame(pathway_id=IDList, pathway_name=hsaList)

if (!dir.exists('kegg')) {
  dir.create('kegg')
}

compoundSymbol <- function(x) {
  cmpList <- keggGet(x)[[1]]$COMPOUND
  if(!is.null(cmpList)){
    listLength <- length(cmpList)
    symbolList <- data.frame(compound_id = names(cmpList), compound = cmpList)
  } else {
    return(NULL)
  }
  symbolList
}

## write each pathway and its compounds to csv in kegg folder
## it often loses connection
## write each record to csv file, so it can start from where it was terminated last time
plyr::l_ply(hsaPathway$pathway_id, function(x) {
  file.out <- file.path('kegg', paste0(x, '.csv'))
  if (!file.exists(file.out)) {
    res <- compoundSymbol(x)
    if (!is.null(res)) {
      res <- data.frame(pathway_id = x, res)
      write.csv(res, file.out)
    }
  }
})

## read csv files
pathway_cmp <- lapply(list.files('kegg'), function(x) {
  file.in <- file.path('kegg', x)
  d <- read.csv(file.in, row.names = 1)
  d
}) %>% do.call(rbind, .)

## join hsaPathway to get pathway names
pathway_cmp <- dplyr::left_join(pathway_cmp, hsaPathway, by = 'pathway_id')
colnames(pathway_cmp)[2] <- 'compound_id'
pathway_cmp <- pathway_cmp %>% dplyr::select(
  pathway_id, pathway_name, compound_id, compound
)

pathway_cmp <- pathway_cmp %>%
  dplyr::mutate(
    kegg_lipid = ifelse(compound_id %in% doc.data$compound_id, 'Yes', '')
  )

## combine pathway and lipid data
pathway_cmp_join <- dplyr::left_join(pathway_cmp, doc.data, by = 'compound_id')

write.csv(pathway_cmp_join, 'KEGG_pathway_compound_hsa.csv')

## total compound and lipids in each pathway
d1 <- pathway_cmp %>% 
  dplyr::group_by(pathway_id, pathway_name) %>%
  dplyr::summarise(N_compound = n(), N_lipid = length(which(kegg_lipid == 'Yes'))) %>%
  dplyr::arrange(desc(N_lipid))

## by lipid parent class
d2 <- pathway_cmp_join %>% 
  dplyr::filter(!is.na(parent)) %>%
  dplyr::group_by(pathway_id, pathway_name, parent) %>%
  dplyr::summarise(N_parent = n())

## by level 1 parent class
d3 <- pathway_cmp_join %>% 
  dplyr::filter(!is.na(`parent level 1`)) %>%
  dplyr::group_by(pathway_id, pathway_name, `parent level 1`) %>%
  dplyr::summarise(N_parent = n())

file.out <- 'KEGG_pathway_lipid_summary.xlsx'

XLConnect::writeWorksheetToFile(
  file.out, d1, sheet = 'summary', rownames = F
)

XLConnect::writeWorksheetToFile(
  file.out, d2, sheet = 'summary_parent', rownames = F
)

XLConnect::writeWorksheetToFile(
  file.out, d3, sheet = 'summary_parent_L1', rownames = F
)


##------------------------
##
##------------------------
library(grid)
library(ggplot2)

p1 <- ggplot(d1, aes(x = N_lipid)) +
  geom_histogram() +
  stat_bin(aes(y=cumsum(..count..)),geom="line", color="grey") +
  labs(x = 'Number of lipids', y = 'Number of pathways',
       title = 'All lipid classes') +
  annotation_custom(
    textGrob(nrow(d1), x = unit(0.9, 'npc'), y = unit(0.9, 'npc')),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    ) +
  theme_bw()

p3 <- lapply(unique(d3$`parent level 1`), function(x) {
  dd <- subset(d3, `parent level 1` == x)
  ggplot(dd, aes(x = N_parent)) +
    geom_histogram() +
    stat_bin(aes(y=cumsum(..count..)),geom="line", color="grey") +
    labs(x = 'Number of lipids', y = 'Number of pathways',
         title = x) +
    annotation_custom(
      textGrob(nrow(dd), x = unit(0.9, 'npc'), y = unit(0.9, 'npc')),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
    ) +
    theme_bw()
})


p <- ggpubr::ggarrange(plotlist = c(list(p1), p3), nrow = 3, ncol = 3)

ggsave('Histogram_KEGG_pathway_lipids.pdf', p, width = 8, height = 8)
ggsave('Histogram_KEGG_pathway_lipids.tiff', p, width = 8, height = 8,
       dpi = 300, compress = 'lzw+p')
