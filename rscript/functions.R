# coln is colnames of data.frame
# return a data.frame with name and class
getAttrs <- function(coln) {
    # name, class, #carbon, #double bond, color of class
    out <- cbind(
        name = as.character(coln),
        class = sapply(c(
            '^PA(?![A-z])', '^LPA(?![A-z])',
            '^PC(?![A-z])', '^LPC(?![A-z])', '^LysoPC(?![A-z])', 'PC$',
            '^PE(?![A-z])', '^LPE(?![A-z])', 'HETE-PE$', '^OX-PE$',
            '^PG(?![A-z])',
            '^PI(?![A-z])', '^LPI(?![A-z])', 
            '^PS(?![A-z])', '^LPS(?![A-z])',
            '^CL(?![A-z])',
            '^Cer(?!1P)', '^Cer1P', 
            '^GluCer(?![A-z])',
            '^GalCer(?![A-z])',
            '^LacCer(?![A-z])',
            '^PhytoCer(?![A-z])',
            '^SM(?![A-z])', '^Lyso-SM(?![A-z])', '^LSM$',
            '^Sph(?![A-z])', '^S1P(?![A-z])',
            '^Gb3(?![A-z])', '^SL(?![A-z])',
            '^Cho$', '^CE(?![A-z])',
            '^DAG(?![A-z])', '^TAG(?![A-z])',  
            '^FFA(?![A-z])', 
            '^PAHSA(?![A-z])', '^FAHFA$',
            'carnitine$',
            '^MGDG(?![A-z])', '^DGDG(?![A-z])',
            '^GM1(?![A-z])', '^GM2(?![A-z])', '^GM3(?![A-z])', 
            '^LBPA(?![A-z])', 
            '^WE(?![A-z])', 
            '^BMP(?![A-z])'),
            function(x) {
                res <- stringr::str_extract(coln, x)
                ## assign Cho to CE
                if (x == '^Cho$')
                    res[which(res == 'Cho')] <- 'CE'
                if (x == '^LysoPC(?![A-z])')
                    res[which(res == 'LysoPC')] <- 'LPC'
                if (x == '^Lyso-SM(?![A-z])')
                    res[which(res == 'Lyso-SM')] <- 'LSM'
                if (x == 'carnitine$')
                    res[which(res == 'carnitine')] <- 'Acylcarnitine'
                if (x == '^PE(?![A-z])')
                    res[res == 'PE' & stringr::str_detect(coln, '\\[O\\]')] <- 'OX-PE'
                if (x == 'HETE-PE$')
                    res[which(res == 'HETE-PE')] <- 'OX-PE'
                if (x == '^PAHSA(?![A-z])')
                    res[which(res == 'PAHSA')] <- 'FAHFA'
                if (x == 'PC$')
                    res[which(res == 'PC')] <- 'OXPC'
                res
            }) %>% apply(1, function(x) {
                x[!is.na(x)][1]
            })
    ) %>% data.frame %>% dplyr::mutate(
        ## sphingolipid and phospholipids are named differently
        carbon = ifelse(class %in% c('SL', 'Cer', 'GM1', 'GM2', 'GM3',
                                     'SM', 'GluCer', 'GalCer',
                                     'LacCer', 'Gb3'),
                        stringr::str_extract(
                            stringr::str_extract(name, '[0-9]{1,3}:[0-9][h)]?$'),
                            '[0-9]+(?=:)'),
                        stringr::str_extract(name, '[0-9]+(?=:)')) %>% as.numeric,
        dbond = ifelse(class %in% c('SL', 'Cer', 'GM1', 'GM2', 'GM3', 'SM', 'GluCer',
                                    'LacCer', 'Gb3'),
                       stringr::str_extract(
                           stringr::str_extract(name, '[0-9]{1,3}:[0-9][h)]?$'),
                           '(?<=:)[0-9]+'),
                       stringr::str_extract(coln, '(?<=:)[0-9]+')) %>% as.numeric
    )  %>% dplyr::mutate(
        carbon_cat = ifelse(
            class %in% c('PC', 'PE', 'PI', 'PS', 'PA', 'PG', 'DAG', 'LBPA'),
            ifelse(carbon >= 44, 'Very long', ifelse(carbon >= 32, 'Long', 'Short')),
            ifelse(class %in% c('TAG'),
                   ifelse(carbon >= 66, 'Very long',
                          ifelse(carbon >= 48, 'Long', 'Short')),
                   ifelse(carbon >= 22, 'Very long',
                          ifelse(carbon >= 16, 'Long', 'Short'))
            )),
        dbond_cat = ifelse(
            class %in% c('PC', 'PE', 'PI', 'PS', 'PA', 'PG', 'DAG', 'LBPA'),
            ifelse(dbond > 4, 'Poly-unsaturated',
                   ifelse(dbond > 0, 'Mono/di-unsaturated', 'Saturated')),
            ifelse(class %in% c('TAG'),
                   ifelse(dbond > 6, 'Poly-unsaturated',
                          ifelse(dbond > 1, 'Mono/di-unsaturated', 'Saturated')),
                   ifelse(dbond > 2, 'Poly-unsaturated',
                          ifelse(dbond > 0, 'Mono/di-unsaturated', 'Saturated')))
        )
    ) %>% dplyr::mutate(
        carbon_cat = factor(
            carbon_cat,
            levels = c('Very long', 'Long', 'Short')
        ),
        dbond_cat = factor(
            dbond_cat,
            levels = c('Poly-unsaturated', 'Mono/di-unsaturated', 'Saturated'))
    )
    # color
    # remove black from color sequence
    out$color <- WGCNA::labels2colors(
        as.numeric(as.factor(out$class)),
        # colorSeq = colorRamps::primary.colors(50, step = 5)[-1])
        colorSeq = WGCNA::standardColors(50)[-7])
    # rowname
    rownames(out) <- out$name
    out
}

# In lipid classes which contain 1 carbon chain
# Short: [0, 16)
# Long:  [16, 22)
# Very long: [22, inf)
#
# Saturated: 0
# Mono/di-unsaturated: 1-2
# Poly-unsaturated: > 2
#
# In lipid classes which contain 2 carbon chains: 'PC', 'PE', 'PI', 'PS', 'PA', 'PG', 'DAG', 'LBPA'
# Short: [0, 32)
# Long:  [32, 44)
# Very long: [44, inf)
#
# Saturated: 0
# Mono/di-unsaturated: [1, 4]
# Poly-unsaturated: [5, Inf)
#
# In lipid classes which contain 3 carbon chains: TAG
# Short: [0, 48)
# Long:  [48, 66)
# Very long: [66, inf)
#
# Saturated: 0
# Mono-unsaturated: [1, 6]
# Poly-unsaturated: [7, Inf)
