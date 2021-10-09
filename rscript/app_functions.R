#------------------------------------
# Games-Howell Post Hoc Test in R
# https://stats.stackexchange.com/questions/83941/games-howell-post-hoc-test-in-r
#-------------------------------------

## obselete
posthoc.tgh <- function(y, x, method=c("games-howell", "tukey"), digits=2) {
    ### Based on http://www.psych.yorku.ca/cribbie/6130/games_howell.R
    method <- tolower(method);
    tryCatch(method <- match.arg(method), error=function(err) {
        stop("Argument for 'method' not valid!");
    });
    
    res <- list(input = list(x=x, y=y, method=method, digits=digits));
    
    res$intermediate <- list(x = factor(x[complete.cases(x,y)]),
                             y = y[complete.cases(x,y)]);
    res$intermediate$n <- tapply(y, x, length);
    res$intermediate$groups <- length(res$intermediate$n);
    res$intermediate$df <- sum(res$intermediate$n) - res$intermediate$groups;
    res$intermediate$means <- tapply(y, x, mean);
    res$intermediate$variances <- tapply(y, x, var);
    
    res$intermediate$pairNames <- combn(levels(res$intermediate$x),
                                        2, paste0, collapse=":");
    
    res$intermediate$descriptives <- cbind(res$intermediate$n,
                                           res$intermediate$means,
                                           res$intermediate$variances);
    rownames(res$intermediate$descriptives) <- levels(res$intermediate$x);
    colnames(res$intermediate$descriptives) <- c('n', 'means', 'variances');
    
    ### Start on Tukey
    res$intermediate$errorVariance <-
        sum((res$intermediate$n-1) * res$intermediate$variances) /
        res$intermediate$df;
    res$intermediate$t <- combn(res$intermediate$groups, 2, function(ij) {
        abs(diff(res$intermediate$means[ij]))/
            sqrt(res$intermediate$errorVariance*sum(1/res$intermediate$n[ij]));
    } );
    res$intermediate$p.tukey <- ptukey(res$intermediate$t*sqrt(2),
                                       res$intermediate$groups,
                                       res$intermediate$df,
                                       lower.tail=FALSE);
    res$output <- list();
    res$output$tukey <- cbind(res$intermediate$t,
                              res$intermediate$df,
                              res$intermediate$p.tukey)                                     
    rownames(res$output$tukey) <- res$intermediate$pairNames;
    colnames(res$output$tukey) <- c('t', 'df', 'p');
    
    ### Start on Games-Howell
    res$intermediate$df.corrected <- combn(res$intermediate$groups, 2, function(ij) {               
        sum(res$intermediate$variances[ij] /
                res$intermediate$n[ij])^2 / 
            sum((res$intermediate$variances[ij] /
                     res$intermediate$n[ij])^2 / 
                    (res$intermediate$n[ij]-1));
    } );
    res$intermediate$t.corrected <- combn(res$intermediate$groups, 2, function(ij) {               
        abs(diff(res$intermediate$means[ij]))/
            sqrt(sum(res$intermediate$variances[ij] /
                         res$intermediate$n[ij]));
    } );    
    res$intermediate$p.gameshowell <- ptukey(res$intermediate$t.corrected*sqrt(2),
                                             res$intermediate$groups,
                                             res$intermediate$df.corrected,
                                             lower.tail=FALSE)  
    res$output$games.howell <- cbind(res$intermediate$t.corrected,
                                     res$intermediate$df.corrected,
                                     res$intermediate$p.gameshowell);
    rownames(res$output$games.howell) <- res$intermediate$pairNames;
    colnames(res$output$games.howell) <- c('t', 'df', 'p');
    
    ### Set class and return object
    class(res) <- 'posthocTukeyGamesHowell';
    return(res);
    
}

print.posthocTukeyGamesHowell <- function(x, digits=x$input$digits, ...) {
    print(x$intermediate$descriptives, digits=digits);
    cat('\n');
    if (x$input$method == 'tukey') {
        print(x$output$tukey);
    }
    else if (x$input$method == 'games-howell') {
        print(x$output$games.howell, digits=digits);
    }
}

##--------------------------------------------------
## Hypothesis test
##--------------------------------------------------
#' Title
#'
#' @param x data.frame. data in the form of sample x variable
#' @param group factor of group
#' @param pair vector or NA. NA for non-paired analysis. Vector of pair id.
#' @param p.adjust.method character. P value adjust method for post-hoc analysis.
#'
#' @return data.frame of hypothesis test results and fold change
#' @export
#'
#' @examples
hypothesisTest <- function(x, group, pair = NA, p.adjust.method = union('none', p.adjust.methods)) {
    ## d1: data
    ## g: group
    p.adjust.method <- match.arg(p.adjust.method)
    
    if (!is.na(pair)) pair <- as.factor(pair)
    if (!is.data.frame(x)) x <- as.data.frame(x)
    
    if (nlevels(group) == 1) return(NULL)
    if (is.na(pair)) {
        ## not paired design
        test_result <- sapply(x, function(x.i) {
            # g_mean <- aggregate(x, by = list(Group = g), mean)
            
            ## combination of group pairs
            g_pair <- combn(levels(group), 2)
            
            g_fold <- apply(g_pair, 2, function(pair.i) {
                mean(x.i[group == pair.i[2]], na.rm = T) /
                    mean(x.i[group == pair.i[1]], na.rm = T)
            })
            
            names(g_fold) <- apply(g_pair, 2, function(x) {
                paste0('Fold: ', paste0(rev(x), collapse = '/'))
            })
            
            ## welch's ANOVA
            ## equivalent to t-test when there are 2 groups
            aov.out <- oneway.test(x.i ~ group)
            
            ## non-parametric
            ## wilcox when #groups == 2
            ## kruskal-wallis when #groups > 2
            if (nlevels(group) > 2)
                nonparam.out <- kruskal.test(x.i ~ group)
            else
                nonparam.out <- wilcox.test(x.i ~ group)
            
            if (nlevels(group) > 2) {
                ## posthoc.tgh in functions.R
                # posthoc.out1 <- posthoc.tgh(x, g)
                posthoc.out <- rosetta::posthocTGH(
                    x.i, group, 
                    p.adjust = ifelse(p.adjust.method %in% c('bh', 'by'), toupper(p.adjust.method), p.adjust.method), 
                    formatPvalue = F
                )
                
                ## replace - by :
                posthoc.names <- apply(
                    combn(levels(group),2), 2, 
                    function(x) {paste0(rev(x), collapse=':')})
                
                posthoc.out$output$tukey <- as.matrix(posthoc.out$output$tukey)
                rownames(posthoc.out$output$tukey) <- posthoc.names
                posthoc.out$output$games.howell <- as.matrix(posthoc.out$output$games.howell)
                rownames(posthoc.out$output$games.howell) <- posthoc.names
                
                ## dunn's test
                ## p value adjust
                
                dunn.out <- FSA::dunnTest(x.i ~ group, method = p.adjust.method)
                dunn.out$res[, 'Comparison'] <- stringr::str_replace(
                    dunn.out$res[, 'Comparison'],
                    ' - ', ':'
                ) ## make it same as Tukey
                ## dunnTest does not honor level order
                dunn.out$res <- dunn.out$res[
                    match(
                        combn(levels(group), 2, function(x) { paste0(sort(x), collapse=":") }),
                        dunn.out$res[, 'Comparison']
                    ),
                    ]
                ## reverse the order
                paste0Rev <- function(x) {paste0(rev(x), collapse = ':')}
                dunn.out$res$Comparison <- combn(levels(group), 2, paste0Rev)
                
                res <- c(`parametric pvalue` = aov.out$p.value,
                         posthoc.out$output$tukey[, ifelse(p.adjust.method == 'none', 'p', 'p.adjusted')],
                         posthoc.out$output$games.howell[, ifelse(p.adjust.method == 'none', 'p', 'p.adjusted')],
                         `non-parametric pvalue` = nonparam.out$p.value,
                         dunn.out$res[, 'P.adj'],
                         g_fold)
                
                ## Need to set names because column of data.frame
                ## is not named
                # names(res)[2:(1+nrow(posthoc.out$output$tukey))] <-
                #     stringr::rownames(posthoc.out$output$tukey)
                
                ## number of pair-wise comparison
                n_comb <- nrow(posthoc.out$output$tukey)
                names(res)[2:(1 + 2 * n_comb)] <- paste0(
                    rep(c('TukeyHSD: ', 'Games-Howell: '), each = n_comb),
                    names(res)[2:(1 + 2 * n_comb)]
                )
                
                names(res)[(3 + 2 * n_comb):(2 + 3 * n_comb)] <- paste0(
                    rep('Dunn: ', n_comb),
                    dunn.out$res$Comparison
                )
                
            } else {
                res <- c(`parametric pvalue` = aov.out$p.value, 
                         `non-parametric pvalue` = nonparam.out$p.value,
                         g_fold)
            }
            
            res 
        }) %>% t
    } else {
        ## paired design
        ## paired t-test
        test_result <- sapply(x, function(x.i) {
            order_pair <- order(pair)
            x_order <- x.i[order_pair]
            g_order <- group[order_pair]
            g1 <- which(g_order == levels(group)[1])
            g2 <- which(g_order == levels(group)[2])
            
            g_fold <- median(x_order[g2] / x_order[g1], na.rm = T)
            
            c(t.test(x_order[g1], x_order[g2], paired = T)$p.value, 
              wilcox.test(x_order[g1], x_order[g2], paired = T)$p.value,
              g_fold)
        }) %>% t
        
        colnames(test_result) <- c('parametric pvalue', 'non-parametric pvalue', paste0('Fold: ', levels(group)[2], '/', levels(group)[1]))
        
    }
    return(test_result)
}


hypothesisTest2 <- function(x, group, pair = NA, 
                            p.adjust.method = c("none", "holm", "hochberg", "bonferroni", 'BH', 'BY'),
                            test = c('T', 'Wilcox', 
                                     'Pairwise T', 'Pairwise Wilcox',
                                     'ANOVA', 
                                     'TukeyHSD', 'Games-Howell', 
                                     'Kruskal-Wallis', 'Dunn', 
                                     'Fold change')) {
    ## p.adjust.method is the intersection of p.adjust.methods and dunn.test::p.adjustment.methods[c(4, 2:3, 5:8, 1)]
    ## dunn.test::p.adjustment.methods[c(4, 2:3, 5:8, 1)] has bh and by, instead of BH and BY.
    p.adjust.method <- match.arg(p.adjust.method)
    
    if (!is.na(pair)) pair <- as.factor(pair)
    if (!is.data.frame(x)) x <- as.data.frame(x)
    
    if (nlevels(group) == 1) {
        warning('There is only 1 group.')
        return(NULL)
    }
    
    if (nlevels(group) > 2 & !is.na(pair)) {
        warning('Paired test is only applicable when number of groups is 2.')
        return(NULL)
    }
    
    if (!is.na(pair)) {
        ## paired test and nlevels(group) is 2
        test_result <- sapply(x, function(x.i) {
            
            test <- intersect(c('T', 'Wilcox', 'Fold change'), test)
            
            if (length(test) > 0) {
                order_pair <- order(pair)
                x_order <- x.i[order_pair]
                g_order <- group[order_pair]
                g1 <- which(g_order == levels(group)[1])
                g2 <- which(g_order == levels(group)[2])
                
                res <- lapply(test, function(test.i) {
                    if (test.i == 'T') {
                        return(structure(t.test(x_order[g1], x_order[g2], paired = T)$p.value, names = 'parametric pvalue'))
                    }
                    
                    if (test.i == 'Wilcox') {
                        return(structure(wilcox.test(x_order[g1], x_order[g2], paired = T)$p.value, names = 'non-parametric pvalue'))
                    }
                    
                    if (test.i == 'Fold change') {
                        return(structure(median(x_order[g2] / x_order[g1], na.rm = T), 
                                         names = paste0('Fold: ', levels(group)[2], '/', levels(group)[1])))
                    }
                })
            }
            
            do.call(c, res)
        })
    } else{
        if (nlevels(group) == 2) {
            test <- intersect(c('T', 'Wilcox', 'Fold change'), test)
        } else { ## more than 2 groups
            test <- intersect(c('ANOVA', 
                                'Pairwise T',
                                'TukeyHSD', 'Games-Howell', 
                                'Kruskal-Wallis', 
                                'Pairwise Wilcox', 
                                'Dunn', 'Fold change'), test)
        }
        
        if (length(test) > 0) {
            test_result <- sapply(x, function(x.i) {
                
                res <- lapply(test, function(test.i) {
                    if (test.i == 'ANOVA' | test.i == 'T') {
                        ## oneway.test is equivalent to t.test when number of groups is 2
                        return(structure(oneway.test(x.i ~ group)$p.value, names = 'parametric pvalue'))
                    }
                    
                    if (test.i == 'Wilcox') {
                        return(structure(wilcox.test(x.i ~ group)$p.value, names = 'non-parametric pvalue'))
                    }
                    
                    if (test.i == 'Fold change') {
                        ## combination of group pairs
                        g_pair <- combn(levels(group), 2)
                        
                        fold_change <- apply(g_pair, 2, function(pair.i) {
                            mean(x.i[group == pair.i[2]], na.rm = T) /
                                mean(x.i[group == pair.i[1]], na.rm = T)
                        })
                        
                        names(fold_change) <- apply(g_pair, 2, function(pair.i) {
                            paste0('Fold: ', paste0(rev(pair.i), collapse = '/'))
                        })
                        
                        return(fold_change)
                    }
                    
                    if (test.i == 'Kruskal-Wallis') {
                        return(structure(kruskal.test(x.i ~ group)$p.value, names = 'non-parametric pvalue'))
                    }
                    
                    if (test.i %in% c('TukeyHSD', 'Games-Howell')) {
                        posthoc.out <- rosetta::posthocTGH(
                            x.i, group, 
                            p.adjust = p.adjust.method, 
                            formatPvalue = F
                        )
                        
                        ## replace - by :
                        posthoc.names <- apply(
                            combn(levels(group),2), 2, 
                            function(x) {paste0(rev(x), collapse=':')})
                       
                        if (test.i == 'TukeyHSD') {
                            
                            posthoc.out$output$tukey <- as.matrix(posthoc.out$output$tukey)
                            rownames(posthoc.out$output$tukey) <- posthoc.names
                            res.TukeyHSD <- posthoc.out$output$tukey[, ifelse(p.adjust.method == 'none', 'p', 'p.adjusted')]
                            names(res.TukeyHSD) <- paste('TukeyHSD:', names(res.TukeyHSD))
                            return(res.TukeyHSD)
                            
                        } else if (test.i == 'Games-Howell') {
                            
                            posthoc.out$output$games.howell <- as.matrix(posthoc.out$output$games.howell)
                            rownames(posthoc.out$output$games.howell) <- posthoc.names
                            res.Games_Howell <- posthoc.out$output$games.howell[, ifelse(p.adjust.method == 'none', 'p', 'p.adjusted')]
                            names(res.Games_Howell) <- paste('Games-Howell:', names(res.Games_Howell))
                            return(res.Games_Howell)
                        }
                        
                        
                    } 
                    
                    if (test.i == 'Dunn') {
                        
                        dunn.out <- FSA::dunnTest(x.i ~ group, method = tolower(p.adjust.method))
                        dunn.out$res[, 'Comparison'] <- stringr::str_replace(
                            dunn.out$res[, 'Comparison'],
                            ' - ', ':'
                        ) ## make it same as Tukey
                        ## dunnTest does not honor level order
                        dunn.out$res <- dunn.out$res[
                            match(
                                combn(levels(group), 2, function(x) { paste0(sort(x), collapse=":") }),
                                dunn.out$res[, 'Comparison']
                            ),
                            ]
                        ## reverse the order
                        paste0Rev <- function(x) {paste0(rev(x), collapse = ':')}
                        dunn.out$res$Comparison <- combn(levels(group), 2, paste0Rev)
                        res.dunn <- dunn.out$res[, 'P.adj']
                        names(res.dunn) <- paste('Dunn:', dunn.out$res[, 'Comparison'])
                        return(res.dunn)
                    }
                    
                    if (test.i == 'Pairwise Wilcox') {
                        g_pair <- combn(levels(group), 2)
                        
                        res.wilcox <- apply(g_pair, 2, function(pair.i) {
                            wilcox.test(x.i[group == pair.i[1]], x.i[group == pair.i[2]])$p.value
                        })
                        
                        names(res.wilcox) <- apply(g_pair, 2, function(pair.i) {
                            paste0('Wilcox: ', paste0(rev(pair.i), collapse = ':'))
                        })
                        
                        return(res.wilcox)
                    }
                    
                    if (test.i == 'Pairwise T') {
                        g_pair <- combn(levels(group), 2)
                        
                        res.t <- apply(g_pair, 2, function(pair.i) {
                            t.test(x.i[group == pair.i[1]], x.i[group == pair.i[2]])$p.value
                        })
                        
                        names(res.t) <- apply(g_pair, 2, function(pair.i) {
                            paste0('T: ', paste0(rev(pair.i), collapse = ':'))
                        })
                        
                        return(res.t)
                    }
                })
                do.call(c, res)
            })
        }
        
    }
    
    return(t(test_result))
}

##--------------------------------------------------
## Volcano plot
##--------------------------------------------------
plot_volcano <- function(d, sig_lvl, fold_cutoff, max_label, title = '', 
                         ggtheme = theme_bw(), font_family = 'Arial',
                         force = 1, expand = 0.05) {
    ## d is data.frame with 3 columns: pvalue, fold, label
    ## sig_lvl is the significance level
    ## fold_cutoff is fold cutoff
    assertthat::assert_that(ncol(d) == 3)
    
    colnames(d) <- c('pvalue', 'fold', 'label')
    
    d <- d %>% dplyr::mutate(
        # pvalue and fold change both must pass threshold
        Significant = ifelse(pvalue < sig_lvl & 
                                 abs(log(fold)) > log(fold_cutoff), 
                             paste0('P<', sig_lvl, '&fold>', fold_cutoff), 'Not sig'),
        label = ifelse(Significant != 'Not sig', label, '')
    ) %>% dplyr::filter(
        !is.na(pvalue) & !is.na(fold)
    )
    
    
    if (max_label > 0) {
        ## number of significant variables
        n_sig <- length(which(d$label != ''))
        if (n_sig > max_label) {
            ## ids of top max_label smallest P value
            id_label <- order(d[, 'pvalue'])[1:max_label]
            d[-id_label, 'label'] <- ''
        }
    }
    # if (max_label > 0) {
    #     ## cut_off is the max_label largest P value or the largest P value
    #     ## of significant variables depending on if max_label is larger than
    #     ## the number of significant variables
    #     if (length(which(d$label != '')) > max_label) {
    #         p_cut <- sort(subset(d, label != '')[, 'pvalue'])[max_label]
    #     } else {
    #         p_cut <- 1.1
    #     }
    #     d <- d %>% dplyr::mutate(
    #         label = ifelse(
    #             pvalue >= p_cut,
    #             '', label)
    #     )
    # }
    
    if (fold_cutoff == 1) {
        ## change significant to P > 0.05 and P < 0.05
        ## change color_values accordingly
        d <- d %>% dplyr::mutate(
            Significant = ifelse(
                Significant == 'Not sig',
                paste0('P\u2265', sig_lvl),
                paste0('P<', sig_lvl)
            )
        )
        color_values <- c('gray40', 'red') %>%
            `names<-`(c(paste0('P\u2265', sig_lvl),
                        paste0('P<', sig_lvl)))
    } else {
        color_values <- c('gray40', 'red') %>%
            `names<-`(c('Not sig', paste0('P<', sig_lvl, '&fold>', fold_cutoff)))
    }
    
    p <- ggplot(d, aes(x = log2(fold), y = -log10(pvalue))) +
        geom_point(aes(color = Significant)) +
        scale_color_manual(values = color_values) +
        ggrepel::geom_text_repel(aes(label = label), family = font_family,
                                 force = force, max.iter = 2000) +
        geom_vline(xintercept = 0, linetype = 'dashed') +
        geom_hline(yintercept = -log10(sig_lvl), linetype = 'dashed') +
        scale_x_continuous(expand = expand_scale(mult = c(0, expand), add = c(0, 0))) +
        scale_y_continuous(expand = expand_scale(mult = c(0, expand), add = c(0, 0))) +
        labs(title = title) +
        ggtheme
    
    p
}

##--------------------------------------------------
## subchunkify
## set figure width and height dynamically within a chunk
## https://stackoverflow.com/questions/15365829/dynamic-height-and-width-for-knitr-plots
##--------------------------------------------------
# 
# subchunkify <- function(g, fig_height=7, fig_width=5) {
#     g_deparsed <- paste0(deparse(
#         function() {g}
#     ), collapse = '')
#     
#     sub_chunk <- paste0("
#                         `","``{r sub_chunk_", floor(runif(1) * 10000), ", fig.height=",
#                         fig_height, ", fig.width=", fig_width, ", echo=FALSE}",
#                         "\n(", 
#                         g_deparsed
#                         , ")()",
#                         "\n`","``
#                         ")
#     
#     cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
                        # }

# And use the function like this, defining your own figure sizes:
#     
# ```{r echo=FALSE, results='asis'}
# g <- ggplot(economics, aes(date, unemploy)) + 
#     geom_line()
# subchunkify(g, 10, 3)
# subchunkify(g, 7, 7)
# ```

# Or let the data define the sizes:
#     
#     ```{r echo=FALSE, results='asis'}
# g <- ggplot(economics, aes(date, unemploy)) + 
#     geom_line()
# for (i in seq(2, 5)) {
#     subchunkify(g, i / 2, i)
# }

# g <- ggplot(economics, aes(date, unemploy)) + 
#     geom_line()
# 
# cat('<h2>A Small Square Plot</h2>')
# subchunkify(g, 3, 3)


##--------------------------------------------------
## boxplot 
##--------------------------------------------------
batch_boxplot <- function(d, g, dh, add_point, add_sig, add_violin, notch, 
                          x_angle, h_just, v_just, posthoc, unit, 
                          ggtheme = theme_bw(), font_family = 'Arial') {
    # d: data.frame
    # g: group, length is the same as row number of d
    # dh: hypothesis test resutls
    d.colnames <- colnames(d)
    d <- cbind(d, group = g)
    
    ## check dh is numeric matrix
    if (is.data.frame(dh)) {
        rn <- rownames(dh)
        x <- as.matrix(dh) %>% 
            apply(2, as.numeric) %>%
            `rownames<-`(rn)
    } else if (is.matrix(dh)) {
        if (!is.numeric(dh)) {
            rn <- rownames(dh)
            dh <- apply(dh, 2, as.numeric) %>%
                `rownames<-`(rn)
        }
    } else {
        stop('dh is not data.frame not matrix.')
    }
    
    p_list <- lapply(d.colnames, function(var.i) {
        ## add violin
        if (add_violin) {
            p <- ggplot(d, aes(x = group, y = eval(parse(text = paste0('`', var.i, '`'))))) +
                geom_violin() +
                geom_boxplot(width = 0.2, notch = notch)
            
            ## add violin and jitter
            if (add_point) {
                p <- ggplot(d, aes(x = group, y = eval(parse(text = paste0('`', var.i, '`'))))) +
                    geom_violin() + 
                    geom_boxplot(width = 0.2, outlier.size = 0, notch = notch) +
                    geom_jitter(width = 0.1)
            }
        } else if (add_point) {
            p <- ggplot(d, aes(x = group, y = eval(parse(text = paste0('`', var.i, '`'))))) +
                geom_boxplot(outlier.size = 0, notch = notch) +
                geom_jitter(width = 0.1)
        } else {
            ## no voilin no jitter
            p <- ggplot(d, aes(x = group, y = eval(parse(text = paste0('`', var.i, '`'))))) +
                geom_boxplot(notch = notch) 
        }
        
        if (add_sig) {
            if (nlevels(g) == 1) {
                
            } else if(nlevels(g) > 2) {
                
                p_matrix <- matrix(
                    1, nrow = nlevels(g), ncol = nlevels(g),
                    dimnames = list(levels(g), levels(g))    
                )
                ## multcompLetters work on lower.tri
                p_matrix[lower.tri(p_matrix)] <- dh[var.i, stringr::str_detect(colnames(dh), posthoc)]
                
                sig_letters <- multcompView::multcompLetters(
                    p_matrix, threshold = 0.05
                )$Letters
                
                d_sig <- data.frame(group = names(sig_letters), 
                                    y = 1.05 * max(d[, var.i]) - 0.05 * min(d[, var.i]),
                                    label = sig_letters)                    
                
                suppressWarnings(
                    p <- p + geom_text(
                        data = d_sig,
                        aes(x = group, y = y, label = label), family = font_family
                    ) + scale_y_continuous(
                        expand = c(0.1, 0)
                    )
                )
            } else {
                p_col <- ifelse(posthoc == 'Parametric',
                                'parametric pvalue', 'non-parametric pvalue')
                d_sig <- data.frame(
                    start = levels(g)[1],
                    end = levels(g)[2],
                    y = 1.1 * max(d[, var.i]) - 0.1 * min(d[, var.i]),
                    label = formatC(dh[var.i, p_col], digits = 2)
                )
                
                suppressWarnings(
                    p <- p + ggsignif::geom_signif(
                        data = d_sig,
                        aes(xmin = start, xmax = end,
                            annotations = label, y_position = y),
                        manual = T,
                        tip_length = min(
                            0.01 * (max(d[, var.i]) - min(d[, var.i])),
                            0.01),
                        family = font_family
                    ) + scale_y_continuous(expand = c(0.15, 0))
                )
            }
            
        }
        
        ## y-axis label
        ylab <- ifelse(unit == ' ',
                       'Concentration',
                       parse(text = paste0('Concentration (', unit, ')'))
        )
        
        p <- p + labs(x = '', title = var.i, 
                      y = ylab) +
            ggtheme +
            theme(axis.text.x = element_text(angle = x_angle, vjust = v_just, hjust = h_just))
        
        p
    })
    
    return(p_list)
}

##--------------------------------------------------
## barplot 
##--------------------------------------------------

batch_barplot <- function(d, g, dh, add_sig, x_angle, h_just, v_just, 
                          posthoc, unit, conf.int = 0.95, 
                          ggtheme = theme_bw(), font_family = 'Arial') {
    # d: data.frame
    # g: group, length is the same as row number of d
    # dh: hypothesis test resutls
    # conf.int: confidence interval
    d.colnames <- colnames(d)
    d <- cbind(d, group = g)
    mult <- qnorm((1 + conf.int)/2)
    if (abs(conf.int - 0.68) < 0.01)
        mult = 1
    
    ## check dh is numeric matrix
    if (is.data.frame(dh)) {
        rn <- rownames(dh)
        x <- as.matrix(dh) %>% 
            apply(2, as.numeric) %>%
            `rownames<-`(rn)
    } else if (is.matrix(dh)) {
        if (!is.numeric(dh)) {
            rn <- rownames(dh)
            dh <- apply(dh, 2, as.numeric) %>%
                `rownames<-`(rn)
        }
    } else {
        stop('dh is not data.frame not matrix.')
    }
    
    p_list <- lapply(d.colnames, function(var.i) {
        
        p <- ggplot(
            d, 
            aes(x = group, 
                y = eval(parse(text = paste0('`', var.i, '`'))))) +                                stat_summary(fun.y = mean, geom = 'bar', 
                                                                                                                fill = 'gray80', color = 'black') +
            stat_summary(fun.data = mean_se, fun.args = list(mult = mult),
                         geom = 'errorbar', width = 0.2)
        
        
        if (add_sig) {
            if (nlevels(g) == 1) {
                
            } else if(nlevels(g) > 2) {
                p_matrix <- matrix(
                    1, nrow = nlevels(g), ncol = nlevels(g),
                    dimnames = list(levels(g), levels(g))    
                )
                ## multcompLetters work on lower.tri
                p_matrix[lower.tri(p_matrix)] <- dh[var.i, stringr::str_detect(colnames(dh), posthoc)]
                
                sig_letters <- multcompView::multcompLetters(
                    p_matrix, threshold = 0.05
                )$Letters
                ## compute y limits
                ymax <- vaggregate(d[, var.i], g, mean_se) %>%
                    unlist %>% max
                d_sig <- data.frame(group = names(sig_letters), 
                                    # y = ymax * 1.05,
                                    # y = 1.05 * max(d[, var.i]) - 0.05 * min(d[, var.i]),
                                    y = 1.1 * max(vaggregate(d[, var.i], g, mean_cl_normal)[3, ] %>% unlist),
                                    label = sig_letters)                    
                
                suppressWarnings(
                    p <- p + geom_text(
                        data = d_sig,
                        aes(x = group, y = y, label = label),
                        family = font_family
                    ) + scale_y_continuous(
                        expand = c(0.1, 0)
                    )
                )
            } else {
                ## 2 groups
                p_col <- ifelse(posthoc == 'Parametric',
                                'parametric pvalue', 'non-parametric pvalue')
                d_sig <- data.frame(
                    start = levels(g)[1],
                    end = levels(g)[2],
                    y = 1.1 * max(vaggregate(d[, var.i], g, mean_cl_normal)[3, ] %>% unlist),
                    label = formatC(dh[var.i, p_col], digits = 2)
                )
                
                suppressWarnings(
                    p <- p + ggsignif::geom_signif(
                        data = d_sig,
                        aes(xmin = start, xmax = end, y_position = y,
                            annotations = label),
                        manual = T,
                        tip_length = min(
                            0.01 * max(vaggregate(d[, var.i], g, mean_cl_normal)[3, ] %>% unlist),
                            0.01),
                        family = font_family
                    ) + scale_y_continuous(expand = c(0.15, 0))
                )
            }
            
        }
        
        ## y-axis label
        ylab <- ifelse(unit == ' ',
                       'Concentration',
                       parse(text = paste0('Concentration (', unit, ')'))
        )
        
        p <- p + labs(x = '', title = var.i, 
                      y = ylab) +
            ggtheme +
            theme(axis.text.x = element_text(angle = x_angle, vjust = v_just, hjust = h_just))
        
        p
    })
    
    return(p_list)
}

##--------------------------------------------------
## grid_fam
## change font family of grid plot 
## 20190530
##--------------------------------------------------
grid_fam <- function(p, font_family="Arial") 
{
    if (class(p)[1] == 'pheatmap') {
        ## return gtable if p is of clas pheatmap
        ## colname and rowname
        # print('1')
        g <- p$gtable
        # print('2')
        px <- which(g$layout$name %in% c('col_names', 'row_names'))
        # print('3')
        for (i in px) {
            g$grobs[[i]]$gp$fontfamily <- font_family
            g$grobs[[i]]$gp$fontface <- 1L
            g$grobs[[i]]$gp$font <- 1L
        }
        # print('4')
        ## legend
        px <- which(g$layout$name == 'legend')
        # print('5')
        id <- grep('text', names(g$grobs[[px]]$children))
        # print('6')
        for (i in id) {
            g$grobs[[px]]$children[[i]]$gp$fontfamily <- font_family
            g$grobs[[px]]$children[[i]]$gp$fontface <- 'plain'
            g$grobs[[px]]$children[[i]]$gp$font <- structure(1L, names = 'plain')
        }
        # print('7')
        return(invisible(g))
    }
    
    ## https://stackoverflow.com/questions/17012518/why-does-this-r-ggplot2-code-bring-up-a-blank-display-device
    ## pdf(NULL) prevents ggplotGrob from generating Rplots.pdf
    pdf(NULL)
    suppressWarnings(g <- ggplotGrob(p))
    
    # dev.off()
    
    px <- which(g$layout$name=="panel")
    
    id <- grep("text", names(g$grobs[[px]]$children))
    
    for(i in id)  {
        n <- length(g$grobs[[px]]$children[[i]]$gp$fontfamily)
        g$grobs[[px]]$children[[i]]$gp$fontfamily <- rep(font_family, n)
        g$grobs[[px]]$children[[i]]$gp$fontface <- rep('plain', n)
        g$grobs[[px]]$children[[i]]$gp$font <- structure(rep(1L, n), 
                                                         names = rep('plain', n))
    }
    # grid::grid.newpage()
    # grid::grid.draw(g)
    
    invisible(g)
}


##--------------------------------------------------
## OPLS-DA splot
## 20190606
## Calculate p1 and pcorr1 for splot
##--------------------------------------------------
splotCal <- function(opls.obj, x) {
    ## opls.obj: ropls::opls object
    ## x: normalized data passed to opls
    ## normalized data
    # s <- as.matrix(mSetObj$dataSet$norm)
    s <- as.matrix(x)
    ## p1
    # T <- as.matrix(mSetObj$analSet$oplsda$scoreMN)
    T <- as.matrix(opls.obj@scoreMN)
    ## p1 is pearson correlation between p1 and variable
    p1 <- c()
    for (i in 1:ncol(s)) {
        scov <- cov(s[, i], T)
        p1 <- matrix(c(p1, scov), ncol = 1)
    }
    ## pcorr1 is pearson correlation divided by standard deviation of each variable
    pcorr1 <- c()
    for (i in 1:nrow(p1)) {
        den <- apply(T, 2, sd) * sd(s[, i])
        corr1 <- p1[i, ]/den
        pcorr1 <- matrix(c(pcorr1, corr1), ncol = 1)
    }
    data.frame(p1, pcorr1, vip = attr(opls.obj, 'vipVn'))
}
