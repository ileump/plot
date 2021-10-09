
#' Perform MSEA and MetPA
#'
#' @param ht.res a data.frame: hypothesis test results
#' @param d.var a data.frame: var.csv, must have HMDB ID, KEGG
#' @param test character: hypothesis test
#' @param analyses character vector: c('MSEA', MetPA
#'
#' @return
#' @export
#'
#' @examples
metabolomicsPathwayAnalysis <- function(
  ht.res, d.var, 
  test = c('parametric pvalue', 'non-parametric pvalue', 'Dunn', 'Games-Howell', 'T', 'Mann-Whitney-U', "Wilcox"), 
  analyses = c('MSEA', 'MetPA'), 
  output.dir, 
  msea.lib.type = 'smpdb_pathway',
  metpa.lib.type
) {
  
  msea.ref.vec <- unique(d.var[, 'Match'])
  msea.ref.type <- 'name'
  metpa.ref.vec <- unique(d.var[, 'KEGG']) %>% Filter(function(x) x != '', .)
  metpa.ref.type <- 'kegg'
  
  testthat::expect_true('HMDB ID' %in% colnames(d.var))
  testthat::expect_true('KEGG' %in% colnames(d.var))
  
  if (!dir.exists(output.dir)) {
    dir.create(output.dir, recursive = T)
  }
  
  ## test.name is used to match the columns in ht.res
  test <- match.arg(test)

  
  ## hypothesis test match test.name
  col.test <- startsWith(
    colnames(ht.res),
    test.name
  ) %>% which
  
  ## fold change columns
  col.fc <- startsWith(colnames(ht.res), 'Fold') %>% which
  
  testthat::expect_equal(
    length(col.test),
    length(col.fc)
  )
  
  if (length(col.test) > 1) {
    for (i in seq_along(col.test)) {
      ## hypothesis test and fold change are consistent for each pair-wise comparison
      name.diff <- Reduce(
        setdiff, 
        c(stringr::str_sub(colnames(ht.res)[col.test[i]], nchar(test.name) + 3, -1),
          stringr::str_sub(colnames(ht.res)[col.fc[i]], nchar('Fold') + 3, -1)) %>%
          stringr::str_split('')
      )
      
      testthat::expect_equal(name.diff, ':')
    }
  }
  
  
  id.hmdb <- d.var[, 'HMDB']
  id.kegg <- d.var[, 'KEGG']
  
  for (i in seq_along(col.test)) {
    ## P value and Fold change
    p <- ht.res[, col.test[i]]
    f <- ht.res[, col.fc[i]]
    
    
    cat('\n\n')
    cat(
      '#',
      stringr::str_sub(colnames(ht.res)[col.fc[i]], nchar('Fold') + 3, -1) %>% stringr::str_trim(),
      '\n\n'
    )
    
    if ('MSEA' %in% analyses) {
      file.msea <- file.path(
        output.dir,
        paste0('MSEA_', stringr::str_sub(colnames(ht.res)[col.fc[i]], nchar('Fold') + 3, -1) %>% stringr::str_trim(), '.csv') %>% fs::path_sanitize(replacement = '_')
      )
      
      ## MSEA
      cat('\n\n')
      cat('## Metabolite Set Enrichment Analysis', '\n\n')
      ids <- which((p < 0.05) & !is.na(id.hmdb) & id.hmdb != 'NA' & nchar(as.character(id.hmdb)) > 1)
      
      if(length(ids) == 0) {
        cat('\n\n')
        cat(sprintf('��????֮?????????????Ĵ?л?P<0.05??????Ϊ%d, ??????ƥ?䵽HMDB ID?Ĵ?л??????Ϊ%d??\n\n', length(which(p < 0.05)), length(ids)))
      } else if (length(ids) > 0) {
        print(
          knitr::kable(
            data.frame(
              name = d.var[ids, 'Metabolite Name'],
              HMDB = d.var[ids, 'HMDB'],
              Class = d.var[ids, 'Class'],
              KEGG = d.var[ids, 'KEGG'],
              P = p[ids] %>% unname,
              FoldChange = round(f[ids], 3) %>% unname
            ),
            caption = 'Table: Metabolites (P < 0.05) used for metabolite set enrichment analysis.'
          ) %>% kable_styling(bootstrap_options = c("striped", "hover"))
        )
        
        tmp <- capture.output(
          msea.res <- DoMSEA(id.hmdb, p , q.type = 'hmdb', 
                             lib.type = msea.lib.type, 
                             ref.vec = msea.ref.vec, 
                             ref.type = msea.ref.type)
        )
        
        ## print table of metabolites passed to MSEA
        if (!is.list(msea.res) && is.numeric(msea.res) && msea.res == 0) {
          cat('\n\n')
          # cat("No match was found to the selected metabolite set library!", '\n\n')
          cat('??л?Ｏ?Ͽ???û???ҵ?ƥ???Ĵ?л?')
        } else {
          msea1 <- msea.res$table
          
          if (!is.null(msea.ref.vec)) {
            cat(sprintf('??л?Ｏ?Ͽ?????%d????л?Ｏ?ϣ??ܹ?????%d????л????????μ??⵽????%d??????????%d???????????졣\n\n',
                        msea.res$mSet.size, msea.res$uniq.count, 
                        length(msea.res$filtered.mset), length(msea.res$ora.hits)
                        ))
          }
         
        
        
        # assign metabolite names in hits column to up-regulated and down-regulated
        msea1 <- msea1 %>% mutate(
          `Up-regulated` = stringr::str_split(Hits, ',') %>%
            sapply(function(x) {
              x <- stringr::str_trim(x)
              paste(x[x %in% d.var[f>1, 'Match']], collapse = ',')
            }),
          `Down-regulated` = stringr::str_split(Hits, ',') %>%
            sapply(function(x) {
              x <- stringr::str_trim(x)
              paste(x[x %in% d.var[f<1, 'Match']], collapse = ',')
            }),
          Hits = NULL
        )
        
        print(
          knitr::kable(
            subset(msea1, fold > 1) %>% dplyr::arrange(`Raw p`),
            row.names = F,
            caption = sprintf('Table: Metabolite sets with fold > 1. Full table can be found in \"%s\"', fs::path_file(file.msea)),
          ) %>% kable_styling(bootstrap_options = c("striped", "hover"))
        )
        
        if(!file.exists(file.msea))
          write.csv(msea1, file.msea, row.names = F)
        }
      }
    }
    
    if ('MetPA' %in% analyses) {
      file.metpa <- file.path(
        output.dir,
        paste0('MetPA_', stringr::str_sub(colnames(ht.res)[col.fc[i]], nchar('Fold') + 3, -1) %>% stringr::str_trim(), '.csv') %>% fs::path_sanitize(replacement = '_')
      )
      
      ## MetPA
      cat('\n\n')
      cat('## Metabolic Pathway Analysis', '\n\n')
      ids <- which((p < 0.05) & !is.na(id.kegg) & id.kegg != 'NA' & nchar(as.character(id.kegg)) > 1)
      
      if(length(ids)==0) {
        cat('\n\n')
        # cat('0 significant metabolites (P < 0.05)', '\n\n')
        cat(sprintf('��????֮?????????????Ĵ?л?P<0.05??????Ϊ%d, ??????ƥ?䵽KEGG ID?Ĵ?л??????Ϊ%d??\n\n', length(which(p < 0.05)), length(ids)))
      } else if (length(ids) > 0) {
        print(
          knitr::kable(
            data.frame(
              name = d.var[ids, 'Metabolite Name'],
              HMDB = d.var[ids, 'HMDB'],
              Class = d.var[ids, 'Class'],
              KEGG = d.var[ids, 'KEGG'],
              P = p[ids] %>% unname,
              FoldChange = round(f[ids], 3) %>% unname
            ),
            caption = 'Table: Metabolites (P < 0.05) used for metabolic pathway analysis.'
          ) %>% kable_styling(bootstrap_options = c("striped", "hover"))
        )
        
        cat('\n\n')
        
        ## yeast sce
        tmp <- capture.output(
          metpa.res <- DoPathway(id.kegg, 
                              p, f,
                              q.type = 'kegg', 
                              ref.vec = metpa.ref.vec,
                              ref.type = metpa.ref.type,
                              lib.type = metpa.lib.type)
        )
        
        
        
        if (nrow(metpa.res$table) == 0) {
          cat('\n\n')
          # cat('Metabolites were not found in any pathway.', '\n\n')
          cat('??лͨ·??û???ҵ?ƥ???Ĵ?л?\n\n')
        } else {
          metpa1 <- metpa.res$table
          
          if (!is.null(metpa.ref.vec)) {
            cat(sprintf('??лͨ·???ݿ?????%d????л??ͨ·???ܹ?????%d????л????????μ??⵽????%d??????????%d???????????졣\n\n',
                        metpa.res$mSet.size, metpa.res$uniq.count, 
                        length(metpa.res$filtered.mset), length(metpa.res$ora.hits)
            ))
          }
          
          metpa1.sig <- subset(metpa1, `Raw p` < 0.05)
          if (nrow(metpa1.sig) > 0) {
            print(
              knitr::kable(
                subset(metpa1, `Raw p` < 0.05),
                row.names = F,
                caption = sprintf('Table: Over-represented pathways with raw P < 0.05. Full table can be found in \"%s\"', fs::path_file(file.metpa)),
              ) %>% kable_styling(bootstrap_options = c("striped", "hover"))
            )
          } else {
            cat('\n\n')
            # cat('\n\n', 'No pathway was found to be significantly over-represented.')
            cat('û???ҵ????????????ߵĴ?лͨ·??\n\n')
          }
          
          cat('\n\n')
          if(!file.exists(file.metpa))
            write.csv(metpa1, file.metpa, row.names = F)
        }
        
      }
    }
  }
}

##------------------------------------------------------
## MSEA
##------------------------------------------------------

# debugonce(download.file)
# http://www.metaboanalyst.ca/resources/libs/compound_db.rds
# http://www.metaboanalyst.ca/resources/libs/syn_nms.rds
# https://www.metaboanalyst.ca/resources/libs/msets/kegg_pathway.rda
# https://www.metaboanalyst.ca/resources/libs/msets/smpdb_pathway.rda
# https://github.com/xia-lab/MetaboAnalystR/issues/69
# debugonce(SetKEGG.PathLib)

#' Calculate hypergeometric test score (modified)
#'
#' The original version did not get q.size right. q.size is the number of metabolites 
#' provided by user that are also in cpmd.db (total 19024 metabolites from compound_db.rds) 
#' in the original code, but there are only 1024 unique metabolites in metabolite set.
#' 
#' @param mSetObj metaboanalyst mSet data object
#'
#' @return metaboanalyst mSet data object with hypergeometric test results
#' @export
#'
#' @examples
CalculateHyperScore2 <- function (mSetObj = NA) 
{
  mSetObj <- MetaboAnalystR:::.get.mSet(mSetObj)
  nm.map <- GetFinalNameMap(mSetObj)
  valid.inx <- !(is.na(nm.map$hmdb) | duplicated(nm.map$hmdb))
  ora.vec <- nm.map$hmdb[valid.inx]
  q.size <- length(ora.vec)
  if (is.na(ora.vec) || q.size == 0) {
    AddErrMsg("No valid HMDB compound names found!")
    return(0)
  }
  # list of mSet name and metabolites
  current.mset <- current.msetlib$member
  if (mSetObj$dataSet$use.metabo.filter && !is.null(mSetObj$dataSet$metabo.filter.hmdb)) {
    current.mset <- lapply(current.mset, function(x) {
      x[x %in% mSetObj$dataSet$metabo.filter.hmdb]
    })
    mSetObj$dataSet$filtered.mset <- current.mset
  }
  uniq.count <- length(unique(unlist(current.mset, use.names = FALSE)))
  set.size <- length(current.mset)
  if (set.size == 1) {
    AddErrMsg("Cannot perform enrichment analysis on a single metabolite set!")
    return(0)
  }
  hits <- lapply(current.mset, function(x) {
    x[x %in% ora.vec]
  })
  hit.num <- unlist(lapply(hits, function(x) length(x)), use.names = FALSE)
  if (sum(hit.num > 0) == 0) {
    AddErrMsg("No match was found to the selected metabolite set library!")
    return(0)
  }
  set.num <- unlist(lapply(current.mset, length), use.names = FALSE)
  res.mat <- matrix(NA, nrow = set.size, ncol = 6)
  rownames(res.mat) <- names(current.mset)
  colnames(res.mat) <- c("total", "expected", "hits", "Raw p", 
                         "Holm p", "FDR")
  ## make sure q.size if the number of metabolites in the database
  q.size <- length(which(ora.vec %in% unique(unlist(current.mset))))
  
  for (i in 1:set.size) {
    res.mat[i, 1] <- set.num[i]
    res.mat[i, 2] <- q.size * (set.num[i]/uniq.count)
    res.mat[i, 3] <- hit.num[i]
    res.mat[i, 4] <- phyper(hit.num[i] - 1, set.num[i], 
                            uniq.count - set.num[i], q.size, lower.tail = F)
  }
  res.mat[, 5] <- p.adjust(res.mat[, 4], "holm")
  res.mat[, 6] <- p.adjust(res.mat[, 4], "fdr")
  res.mat <- res.mat[hit.num > 0,  ,drop = F]
  ord.inx <- order(res.mat[, 4])
  mSetObj$analSet$ora.mat <- signif(res.mat[ord.inx, ,drop = F], 3)
  mSetObj$analSet$ora.hits <- hits
  write.csv(mSetObj$analSet$ora.mat, file = "msea_ora_result.csv")
  if (MetaboAnalystR:::.on.public.web) {
    .set.mSet(mSetObj)
    return(1)
  }
  return(MetaboAnalystR:::.set.mSet(mSetObj))
}

#' Metabolite set enrichment analysis - over-representation analysis.
#'
#' @param cmpd.vec vector: compound list.
#' @param lib name of metabolite set library
#' @param q.type type of metabolite ID
#' @param ref.vec compound list of reference metabolome
#' @param ref.type reference metabolome metabolite ID type
#' @param plot a boolean: whether a plot is produced
#' @param use.CalculateHyperScore2 a boolean: whether to use modified CalculateHyperScore
#'
#' @return mSet object
#' @export
#'
#' @examples
msetora <- function(
  cmpd.vec, 
  q.type = c('name', 'hmdb', 'kegg', 'pubchem', 'chebi', 'metlin', 'hmdb_kegg'),
  lib.type = c('smpdb_pathway', 'kegg_pathway_new', 'kegg_pathway_old', 'predicted'), 
  ref.vec = NULL, 
  ref.type = c('name', 'hmdb', 'kegg', 'pubchem', 'chebi', 'metlin'),
  plot = T,
  use.CalculateHyperScore2 = T,
  use.shortName = F
) {
  
  q.type <- match.arg(q.type)
  lib.type <- match.arg(lib.type)
  ref.type <- match.arg(ref.type)
  
  # Create mSetObj
  if (exists('mSet')) rm(mSet)
  mSet<-InitDataObjects("conc", "msetora", FALSE)
  
  #Set up mSetObj with the list of compounds
  mSet<-Setup.MapData(mSet, cmpd.vec)
  
  # Cross reference list of compounds against libraries (hmdb, pubchem, chebi, kegg, metlin)
  mSet<-CrossReferencing(mSet, q.type)
  ## mSet$name.map
  
  # Create the mapping results table
  mSet<-CreateMappingResultTable(mSet)
  
  if (!is.null(ref.vec)) {
    cmpd.db <- MetaboAnalystR:::.read.metaboanalyst.lib("compound_db.rds")

    if (ref.type == 'kegg') {
      hits <- tolower(ref.vec) %in% tolower(cmpd.db$kegg_id)
      mSet$dataSet$metabo.filter.kegg <- ref.vec[hits]
    } else if (ref.type == 'name') {
      hits <- tolower(ref.vec) %in% tolower(cmpd.db$name)
      ## Setup.HMDBReferenceMetabolome
      ## assign compound name, not HMDB ID to mSet$dataSet$metabo.filter.hmdb
      mSet$dataSet$metabo.filter.hmdb <- ref.vec[hits]
    } else {
      ## match ref.type againt cmpd.db column name
      hits <- match(tolower(ref.vec),  cmpd.db[, which(startsWith(colnames(cmpd.db), ref.type))[1]])
      ## KEGG IDs are used here because cmpd.db has version 4 HMDB ID
      mSet$dataSet$metabo.filter.kegg <- cmpd.db[hits, 'kegg_id']
    }
  }
  
  # Set the metabolite filter
  mSet<-SetMetabolomeFilter(mSet, !is.null(ref.vec))
  
  # Select metabolite set library, refer to 
  # Only use metabolite sets containing at least
  mSet<-SetCurrentMsetLib(mSet, lib.type, 2)
  
  # Calculate hypergeometric score, results table generated in your working directory
  if (use.CalculateHyperScore2) {
    mSet<-CalculateHyperScore2(mSet)
  } else {
    mSet <- CalculateHyperScore(mSet)
  }
  
  ## current.mset intersects reference metabolome
  # length(unique(unlist(mSet$dataSet$filtered.mset)))
  ## current.mset intersects compound list
  # length(unique(unlist(mSet$analSet$ora.hits)))

  
  if(is.numeric(mSet)) {
    if (mSet == 0)
      return(0)
  }   
  
  # Plot the ORA, bar-graph
  # mSet<-PlotORA(mSet, "ora_0_", "bar", "png", 72, width=NA)
  
  folds <- mSet$analSet$ora.mat[, 3]/mSet$analSet$ora.mat[, 2]
  
  ## use short metabolite names
  if (use.shortName) {
    names(folds) <- MetaboAnalystR:::GetShortNames(rownames(mSet$analSet$ora.mat))
  } else {
    names(folds) <- rownames(mSet$analSet$ora.mat)
  }
  
  pvals <- mSet$analSet$ora.mat[, 4]
  
  folds <- folds[folds > 1]
  pvals <- pvals[folds > 1]
  
  if (length(folds) > 0 & plot)
    PlotMSEA.Overview(folds, pvals)
  
  mSet
}


#' A wrapper for MSEA analysis
#'
#' @param cmpd.vec Character vector of compound list
#' @param p Numeric vector of P values
#' @param q.type Metabolite ID type
#' @param lib.type Metabolite set library type
#' @param ref.vec Character vector: a list of metabolite IDs as reference metabolome.
#' @param ref.type Character: ID type of metabolites in reference metabolome
#' @param quiet Print MSEA result table
#'
#' @return mSet object of MSEA results
#' @export
#'
#' @examples
DoMSEA <- function(cmpd.vec, p, q.type, lib.type, ref.vec, ref.type, quiet = T) {
  # cat('\n\n##### ', colnames(p)[i], ' ', lib.i, '\n\n')
  ## remove NA and empty string from metabolite ID
  cmpd.vec <- setdiff(cmpd.vec[p < 0.05] %>% unique, c('NA', '', NA))
  
  mSet <- msetora(cmpd.vec, q.type = q.type, lib.type = lib.type, 
               ref.vec = ref.vec, ref.type = ref.type)
  
  if(is.numeric(mSet) && mSet == 0) {
    return(0)
  }
  cat('\n\nMap metabolite names to reference names \n\n')
  
  out.table <- table(mSet[['name.map']][['match.state']])
  names(out.table) <- plyr::mapvalues(names(out.table), 
                                      from = c('1', '0'),
                                      to = c('Matched', 'Unmatched'))
  knitr::kable(
    out.table ,
    col.names = c('', 'Freq')
  ) %>% print
  
  if (any(mSet[['name.map']][['match.state']] == 0))
    knitr::kable(
      mSet[['name.map']][['query.vec']][mSet[['name.map']][['match.state']] == 0],
      caption = 'Not matched'
    ) %>% print
  
  cat('\n\nEnriched metabolite sets with fold > 1 \n\n')
  
  out <- mSet$analSet$ora.mat %>%
    as.data.frame %>%
    `rownames<-`(rownames(mSet$analSet$ora.mat))
  
  out$fold <- round(out$hits / out$expected, 2)
  
  out <- out[, c('total', 'expected', 'hits', 'fold', 'Raw p', 'Holm p', 'FDR')]
  
  ## order by significance
  out <- out[order(out[, 'Raw p'], decreasing = F), ]
  
  if(!quiet) {
    knitr::kable(
      out
    ) %>% print
  }
  
  
  if (!is.null(ref.vec)) {
    db1 <- plyr::ldply(mSet$dataSet$filtered.mset, function(x) {
      paste0(x, collapse = ', ')
    }) %>% `colnames<-`(c('id', 'Total'))
  }

  
  db2 <- plyr::ldply(mSet$analSet$ora.hits, function(x) {
    paste0(x, collapse = ', ')
  }) %>% `colnames<-`(c('id', 'Hits'))
  
  db <- merge(db1, db2, by = 'id')
  
  out[, 'id'] <- rownames(out)
  
  out <- merge(out, db, by = 'id', all.x = T)
  
  colnames(out)[colnames(out) == 'id'] <- 'Metabolite Set'
  # out[, 'id'] <- NULL
  rownames(out) <- NULL
  
  list(table = out,
       mSet.size = length(current.msetlib$member),
       uniq.count = mSet$dataSet$uniq.count,
       filtered.mset = unique(unlist(mSet$dataSet$filtered.mset)),
       ora.hits = unique(unlist(mSet$analSet$ora.hits)))
}




##------------------------------------------------------
## MetPA
##------------------------------------------------------

#' Make bubble plot for MetPA analysis 
#'
#' @param mSetObj 
#'
#' @return
#' @export
#'
#' @examples
plotPathSummary.2 <- function(mSetObj = NA) {
  mSetObj <- MetaboAnalystR:::.get.mSet(mSetObj)
  if (mSetObj$analSet$type == "pathora") {
    x <- mSetObj$analSet$ora.mat[, 8]
    y <- mSetObj$analSet$ora.mat[, 4]
  } else {
    x <- mSetObj$analSet$qea.mat[, 7]
    y <- mSetObj$analSet$qea.mat[, 3]
  }
  if (length(y) < 2)
    return(0)
  y = -log10(y)
  inx <- order(y, decreasing = T)
  x <- x[inx]
  y <- y[inx]
  sqx <- sqrt(x)
  min.x <- min(sqx)
  max.x <- max(sqx)
  maxR <- (max.x - min.x)/40
  minR <- (max.x - min.x)/160
  radi.vec <- minR + (maxR - minR) * (sqx - min.x)/(max.x - 
                                                      min.x)
  bg.vec <- heat.colors(length(y))
  
  d <- data.frame(
    x = x,
    y = y,
    radi = radi.vec,
    bg = bg.vec,
    id = names(x),
    label = names(metpa$path.ids)[match(names(x), metpa$path.ids)]
  )
  
  ## only label top 6 in terms of P value
  d[y <= -log10(0.05), 'label'] <- ''
  g <- ggplot(d, aes(x = x, y = y, fill = I(bg), size = I(radi * 900))) +
    geom_point(shape=21, color = 'black') +
    ggrepel::geom_text_repel(aes(label = label), size = 3, color = 'black') +
    labs(x = 'Pathway Impact', y = '-log10(p)') +
    theme_bw()
  
  print(g)
}

#' Title
#'
#' @param cmpd.vec
#' @param organism SetKEGG.PathLib organism code (KEGG)
#' @param q.type Character: metabolite ID type
#' @param nodeImp CalculateOraScore, Indicate the pathway topology analysis, "rbc" for relative-betweeness centrality, and "dgr" for out-degree centrality.
#' @param method CalculateOraScore, is "fisher" or "hyperg"
#' @param ref.vec Character vector: list of metabolites as reference metabolome
#' @param ref.type Character: type of metabolite ID in reference metabolome
#'
#' @return
#' @export
#'
#' @examples
pathora <- function(cmpd.vec, lib.type, 
                    q.type = c('name', 'hmdb', 'kegg', 'pubchem', 'chebi', 'metlin', 'hmdb_kegg'), 
                    nodeImp = c('rbc', 'dgr'), method = c('fisher', 'hyperg'),
                    ref.vec = NULL, ref.type = c('kegg', 'hmdb')) {
  ## organism: -SetKEGG.PathLib organism code (KEGG)
  ## q.type: CrossReferencing, Input the query type, "name" for compound names, 
  ##        "hmdb" for HMDB IDs, "kegg" for KEGG IDs, 
  ##        "pubchem" for PubChem CIDs, "chebi" for ChEBI IDs, 
  ##        "metlin" for METLIN IDs, and "hmdb_kegg" for a both KEGG and HMDB IDs.
  ## nodeImp: CalculateOraScore, Indicate the pathway topology analysis, 
  ##        "rbc" for relative-betweeness centrality, and "dgr" for out-degree centrality.
  ## method: CalculateOraScore, is "fisher" or "hyperg"
  ## ref.vec: refer to Setup.KEGGReferenceMetabolome, Setup.HMDBReferenceMetabolome
  ##        a vector of KEGG/HMDB ID as reference metabolome
  ## ref.type: 'hmdb' or 'kegg'
  
  q.type <- match.arg(q.type)
  nodeImp <- match.arg(nodeImp)
  method <- match.arg(method)
  ref.type <- match.arg(ref.type)
  
  # Create mSetObj for storing objects created during your analysis
  mSet<-InitDataObjects("conc", "pathora", FALSE)
  
  # Set up mSetObj with the list of compounds
  mSet<-Setup.MapData(mSet, cmpd.vec);
  
  
  
  if (!is.null(ref.vec)) {
    cmpd.db <- MetaboAnalystR:::.read.metaboanalyst.lib("compound_db.rds")

    if (ref.type == 'kegg') {
      hits <- tolower(ref.vec) %in% tolower(cmpd.db$kegg_id)
      mSet$dataSet$metabo.filter.kegg <- ref.vec[hits]
    } else if (ref.type == 'hmdb') {
      hits <- tolower(ref.vec) %in% tolower(cmpd.db$name)
      mSet$dataSet$metabo.filter.hmdb <- ref.vec[hits]
    }
    
  }
  
  # Cross reference list of compounds against libraries (hmdb, pubchem, chebi, kegg, metlin)
  mSet<-CrossReferencing(mSet, q.type);
  
  mSet<-CreateMappingResultTable(mSet);
  
  # Select the pathway library, ranging from mammals to prokaryotes
  mSet<-SetKEGG.PathLib(mSet, lib.type, lib.version = 'current')
  
  # Set the metabolite filter
  # F when ref.vec is null
  mSet<-SetMetabolomeFilter(mSet, !is.null(ref.vec));
  
  # Calculate the over representation analysis score, here we selected to use the hypergeometric test (alternative is Fisher's exact test)
  # A results table "pathway_results.csv" will be created and found within your working directory
  mSet<-CalculateOraScore2(mSet, nodeImp, method)
  
  # Plot of the Pathway Analysis Overview 
  # mSet<-PlotPathSummary.2(mSet, "path_view_0_", "png", 72, width=NA)
  plotPathSummary.2(mSet)
  
  # Plot a specific metabolic pathway, in this case "Glycine, serine and threonine metabolism"
  # mSet<-PlotMetPath(mSet, "Glycine, serine and threonine metabolism", 480, 420)
  
  mSet
}

DoPathway <- function(cmpd.vec, p, f, q.type = 'kegg', ref.vec, ref.type, lib.type) {
  # cmpd.vec: vector of KEGG ID
  # p: vector of P values
  # f: vector of fold change
  
  # filter metabolite names/ID which have P value less than 0.05
  # remove metabolite names that are NA or blank
  x <- setdiff(cmpd.vec[p < 0.05] %>% unique, c('NA', '', NA))
  
  # run MetaboAnalyst pathway overrepresentation analysis
  # mSet <- pathway(d, organism, q.type)
  ref.vec <- ref.vec[!is.na(ref.vec) & ref.vec != '']
  mSet <- pathora(x, lib.type, q.type = q.type, 
    nodeImp = 'rbc', method = 'hyperg',
    ref.vec = ref.vec, ref.type = ref.type)
  
  cat('\n\nMap metabolite names to reference names \n\n')
  
  out.table <- table(mSet[['name.map']][['match.state']])
  names(out.table) <- plyr::mapvalues(names(out.table), 
                                      from = c('1', '0'),
                                      to = c('Matched', 'Unmatched'))
  knitr::kable(
    out.table ,
    col.names = c('', 'Freq')
  ) %>% print
  
  cat('\n\nEnriched pathways with raw p < 0.05 \n\n')
  
  out <- data.frame(mSet$analSet$ora.mat, check.names = F)
  path.ids <- data.frame(
    id = metpa$path.ids,
    pathway = names(metpa$path.ids)
  ) %>% `rownames<-`(metpa$path.ids)
  
  out$pathway <- path.ids[rownames(out), 'pathway']
  
  colnames(out)[5] <- '-log10(p)'
  out[, '-log10(p)'] <- -log10(out[, 'Raw p'])
  
  db1 <- plyr::ldply(mSet$analSet$ora.hits, function(x) {
    if (length(x) == 0) {
      return(c('', ''))
    } else {
      ## upregulated with fc > 1
      x.up <- x[x %in% cmpd.vec[f > 1]]
      x.down <- x[x %in% cmpd.vec[f < 1]]
      return(c(
        ifelse(length(x.up) == 0,
               '',
               paste0(paste0(names(x.up), ' (', x.up, ')'), collapse = ', ')),
        ifelse(length(x.down) == 0,
               '',
               paste0(paste0(names(x.down), ' (', x.down, ')'), collapse = ', '))
      ))
    }
  }) %>% `colnames<-`(c('id', 'Up-regulated', 'Down-regulated'))
  
  out[, 'id'] <- rownames(out)
  
  out <- merge(out, db1, by = 'id', all.x = T)
  
  
  list(table = out,
       mSet.size = length(metpa$mset.list),
       uniq.count = length(unique(unlist(metpa$mset.list))),
       filtered.mset = unique(unlist(mSet$dataSet$metabo.filter.kegg)),
       ora.hits = unique(unlist(mSet$analSet$ora.hits)))
}

CalculateOraScore2 <- function (mSetObj = NA, nodeImp, method) 
{
  mSetObj <- MetaboAnalystR:::.get.mSet(mSetObj)
  nm.map <- GetFinalNameMap(mSetObj)
  if (mSetObj$pathwaylibtype == "KEGG") {
    valid.inx <- !(is.na(nm.map$kegg) | duplicated(nm.map$kegg))
    ora.vec <- nm.map$kegg[valid.inx]
  }
  else if (mSetObj$pathwaylibtype == "SMPDB") {
    valid.inx <- !(is.na(nm.map$hmdbid) | duplicated(nm.map$hmdbid))
    ora.vec <- nm.map$hmdbid[valid.inx]
  }
  q.size <- length(ora.vec)
  if (is.na(ora.vec) || q.size == 0) {
    if (mSetObj$pathwaylibtype == "KEGG") {
      AddErrMsg("No valid KEGG compounds found!")
    }
    else if (mSetObj$pathwaylibtype == "SMPDB") {
      AddErrMsg("No valid SMPDB compounds found!")
    }
    return(0)
  }
  current.mset <- metpa$mset.list
  uniq.count <- metpa$uniq.count
  if (mSetObj$dataSet$use.metabo.filter && !is.null(mSetObj$dataSet$metabo.filter.kegg)) {
    current.mset <- lapply(current.mset, function(x) {
      x[x %in% mSetObj$dataSet$metabo.filter.kegg]
    })
    mSetObj$analSet$ora.filtered.mset <- current.mset
    uniq.count <- length(unique(unlist(current.mset, use.names = FALSE)))
  }
  hits <- lapply(current.mset, function(x) {
    x[x %in% ora.vec]
  })
  hit.num <- unlist(lapply(hits, function(x) {
    length(x)
  }), use.names = FALSE)
  set.size <- length(current.mset)
  set.num <- unlist(lapply(current.mset, length), use.names = FALSE)
  res.mat <- matrix(0, nrow = set.size, ncol = 8)
  rownames(res.mat) <- names(current.mset)
  colnames(res.mat) <- c("Total", "Expected", "Hits", "Raw p", 
                         "-log(p)", "Holm adjust", "FDR", "Impact")
  if (nodeImp == "rbc") {
    imp.list <- metpa$rbc
    mSetObj$msgSet$topo.msg <- "Your selected node importance measure for topological analysis is \\textbf{relative betweenness centrality}."
  }
  else {
    imp.list <- metpa$dgr
    mSetObj$msgSet$topo.msg <- "Your selected node importance measure for topological analysis is \\textbf{out degree centrality}."
  }
  
  ## make sure metabolites exist in the database
  q.size <- length(which(ora.vec %in% unique(unlist(current.mset))))
  
  res.mat[, 1] <- set.num
  res.mat[, 2] <- q.size * (set.num/uniq.count)
  res.mat[, 3] <- hit.num
  if (method == "fisher") {
    res.mat[, 4] <- GetFisherPvalue(hit.num, q.size, set.num, 
                                    uniq.count)
    mSetObj$msgSet$rich.msg <- "The selected over-representation analysis method is \\textbf{Fishers' exact test}."
  }
  else {
    res.mat[, 4] <- phyper(hit.num - 1, set.num, uniq.count - 
                             set.num, q.size, lower.tail = F)
    mSetObj$msgSet$rich.msg <- "The selected over-representation analysis method is \\textbf{Hypergeometric test}."
  }
  res.mat[, 5] <- -log(res.mat[, 4])
  res.mat[, 6] <- p.adjust(res.mat[, 4], "holm")
  res.mat[, 7] <- p.adjust(res.mat[, 4], "fdr")
  res.mat[, 8] <- mapply(function(x, y) {
    sum(x[y])
  }, imp.list, hits)
  res.mat <- res.mat[hit.num > 0, , drop = FALSE]
  res.mat <- res.mat[!is.na(res.mat[, 8]), , drop = FALSE]
  if (nrow(res.mat) > 1) {
    ord.inx <- order(res.mat[, 4], res.mat[, 8])
    res.mat <- res.mat[ord.inx, ]
  }
  mSetObj$analSet$ora.mat <- signif(res.mat, 5)
  mSetObj$analSet$ora.hits <- hits
  mSetObj$analSet$node.imp <- nodeImp
  MetaboAnalystR:::.set.mSet(mSetObj)
  save.mat <- mSetObj$analSet$ora.mat
  rownames(save.mat) <- GetORA.pathNames(mSetObj)
  write.csv(save.mat, file = "pathway_results.csv")
  if (MetaboAnalystR:::.on.public.web) {
    return(1)
  }
  else {
    return(MetaboAnalystR:::.set.mSet(mSetObj))
  }
}

##------------------------------------------------------
## utilities
##------------------------------------------------------
idMap <- function(ids) {
  ## Remove NA and empty
  ids <- ids[!is.na(ids) & ids != '']
  out <- NULL
  
  id_type <- lapply(ids, function(x) {
    if (stringr::str_detect(x, '^HMDB')) {
      return(c(x, 'hmdb'))
    }
    
    if (stringr::str_detect(x, '^KEGG')) {
      return(c(stringr::str_split_fixed(x, ' ', 2)[2], 'kegg'))
    }
    
    if (stringr::str_detect(x, '^METLIN')) {
      return(c(stringr::str_split_fixed(x, ' ', 2)[2], 'metlin'))
    }
    
    warning(x, ' does not match any known database ID pattern.')
    
    return(c(x, ''))
  }) %>% do.call(rbind, .)
  
  ids <- id_type[, 1]
  id.type <- id_type[, 2]
  # if (length(id.type) == 1) {
  #   ## make use of metaboanalyst CrossReferencing function to map ids
  #   mSet <- InitDataObjects("conc", "msetora", FALSE)
  #   
  #   #Set up mSetObj with the list of compounds
  #   mSet<-Setup.MapData(mSet, ids)
  #   
  #   # Cross reference list of compounds against libraries (hmdb, pubchem, chebi, kegg, metlin)
  #   mSet<-CrossReferencing(mSet, id.type)
  #   ## mSet$name.map
  #   
  #   # Create the mapping results table
  #   mSet <- CreateMappingResultTable(mSet)
  #   
  #   out <- mSet$dataSet$map.table
  #   
  #   out[which(out[, 'KEGG'] %in% c('NA')), 'KEGG'] <- ''
  #   return(out)
  # }
  
  # if (length(ids) != length(id.type))
  #   stop('Lengths are not equal')
  
  for (type.i in unique(id.type)) {
    if (type.i %in% c('name', 'hmdb', 'pubchem', 'chebi', 'kegg', 'metlin')) {
      ind <- which(id.type == type.i)
      
      ## make use of metaboanalyst CrossReferencing function to map ids
      mSet <- InitDataObjects("conc", "msetora", FALSE)
      
      #Set up mSetObj with the list of compounds
      mSet<-Setup.MapData(mSet, ids[ind])
      
      # Cross reference list of compounds against libraries (hmdb, pubchem, chebi, kegg, metlin)
      mSet<-CrossReferencing(mSet, type.i)
      ## mSet$name.map
      
      # Create the mapping results table
      mSet <- CreateMappingResultTable(mSet)
      
      if (is.null(out))
        out <- mSet$dataSet$map.table
      else
        out <- rbind(out, mSet$dataSet$map.table)
    } else {
      warning(type.i, 'is not supported. Use hmdb, pubchem, chebi, kegg, metlin.')
    }
  }
  
  out <- out[match(ids, out[, 'Query']), ]
  
  ## replace NA and 'NA' to ''
  for (col.i in c('Match', 'HMDB', 'PubChem', 'KEGG', 'SMILES', 'Comment')) {
    ind <- out[, col.i] == 'NA' | is.na(out[, col.i])
    out[ind, col.i] <- ''
  }
  out
}

## Get KEGG pathway name for KEGG IDs
getPathwayName <- function(x, organism, q.type) {
  # Create mSetObj for storing objects created during your analysis
  mSet<-InitDataObjects("conc", "pathora", FALSE)
  
  # Set up mSetObj with the list of compounds
  mSet<-Setup.MapData(mSet, x);
  
  # Cross reference list of compounds against libraries (hmdb, pubchem, chebi, kegg, metlin)
  mSet<-CrossReferencing(mSet, q.type, kegg = T, hmdb = T);
  
  mSet<-CreateMappingResultTable(mSet);
  
  # Select the pathway library, ranging from mammals to prokaryotes
  mSet<-SetKEGG.PathLib(mSet, organism)
  
  out <- data.frame(original_name = x, mapped_name = mSet$name.map$hit.values)
  
  ## pathway id, metabolite
  ## long format
  y <- do.call(c, metpa$mset.list) %>% 
    names %>% 
    stringr::str_split('\\.') %>%
    do.call('rbind', .)
  rownames(y) <- y[, 2]
  colnames(y) <- c('pathway_id', 'mapped_name')
  
  ## pathway id, pathway name
  z <- data.frame(metpa$path.ids, names(metpa$path.ids))
  colnames(z) <- c('pathway_id', 'pathway_name')
  
  path <- merge(y, z, by = 'pathway_id')
  
  out <- merge(out, path, by = 'mapped_name', all.x = T)
  out
}


##------------------------------------------------------
## Pathview
##------------------------------------------------------
pathviewAnnotatePValue <- function(cpd.sig = NULL, pathway.id, species = 'hsa',
                                   out.suffix = '', in.dir, out.dir, kegg.dir) {
  ## Annotate significant metabolite with *
  ## cpd.sig is a named boolean vector 
  ## in.dir is the directory of pathway maps generated by pathview 
  ## out.dir is the directory to save annotated pathway map
  ## kegg.dir is the directory of kegg pathway xml files
  
  img.file <- file.path(in.dir, 
                        paste0(species, pathway.id, '.', out.suffix, '.png'))
  kegg.file <- file.path(kegg.dir, paste0(species, pathway.id, '.xml'))
  out.file <- file.path(out.dir, 
                        paste0(species, pathway.id, '.', out.suffix, '.png'))
  
  if (!file.exists(img.file))
    stop(img.file, ' does not exist!')
  if (!file.exists(kegg.file))
    stop(kegg.file, 'does not exist!')
  
  ## read png file
  img <- png::readPNG(img.file)
  ## parse kegg pathway xml
  node <- pathview:::node.info(kegg.file)
  
  cpd.names <- node$labels[node$type == 'compound']
  
  if (any(cpd.sig))
    cpd.sig <- cpd.sig[cpd.sig]
  ## match compound names with pathway node labels
  ## keep none-NA values
  ind <- match(names(cpd.sig), node$labels) %>% Filter(Negate(is.na), .)
  
  ## output image parameters
  width <- ncol(img)
  height <- nrow(img)
  res <- 300
  
  ## print the original image first
  png(out.file, width = width, height = height, res = res)
  par(mar = c(0, 0, 0, 0))
  plot(c(0, width), c(0, height), type = "n", xlab = "", 
       ylab = "", xaxs = "i", yaxs = "i")
  rasterImage(img, 0, 0, width, height, interpolate = F)
  
  ## add annotation
  if(length(ind) > 0)
    for (i in 1:length(ind)) {
      points(x = node$x[ind] + 8, y = height - node$y[ind] + 8, 
             pch = '*', col = 'red', cex = 0.8)
    }
  
  dev.off()
}


pathviewWrap <- function(cpd.data, cpd.name, p, fold, pathways, organism,
                         out.suffix, kegg.dir, out.dir.1, out.dir.2) {
  ## kegg.dir is relative path from out.dir
  ## out.dir.1 is the relative path to store pathview map
  ## out.dir.2 is the relative path to store annotated pathway map
  ## pathview-generated pathway maps are saved in current directory
  # library(parallel)
  
  # Calculate the number of cores
  # no_cores <- detectCores() - 1
  
  # Initiate cluster
  # cl <- makeCluster(no_cores, type="FORK")
  
  if (!(length(cpd.name) == length(p) & 
        length(cpd.data) == length(cpd.name) &
        length(p) == length(fold)))
    stop('Length are not all equal')
  
  
  ## keep compound data with KEGG name
  cpd.data <- structure(fold, names = cpd.name)
  cpd.sig <- structure(p, names = cpd.name)
  
  if (!dir.exists(out.dir.1))
    dir.create(out.dir.1, recursive = T)
  if (!dir.exists(out.dir.2))
    dir.create(out.dir.2, recursive = T)
  
  ## names(table) remove NA
  lapply(
    pathways, 
    function(pw.i) {
      tryCatch(
        pv.out <- pathview::pathview(
          cpd.data = log2(cpd.data),
          pathway.id = pw.i, species = organism, 
          out.suffix = out.suffix,
          keys.align = "y", kegg.native = T, 
          kegg.dir = kegg.dir,
          key.pos = 'topright'
        ),
        error = function(e) {print(e); print('pw.i')}
      )
      
      
      ## move pathway png file to out.dir.1
      f <- paste0(organism, pw.i, '.', out.suffix, '.png')
      if (file.exists(f)) {
        if(file.rename(f, file.path(out.dir.1, f)))
          pathviewAnnotatePValue(
            cpd.sig, pathway.id = pw.i, species = organism, 
            out.suffix = out.suffix, in.dir = out.dir.1,
            out.dir = out.dir.2,
            kegg.dir = kegg.dir
          )
      }
    })
  
  
  # stopCluster(cl)
  invisible()
}

