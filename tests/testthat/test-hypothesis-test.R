context("Testing hypothesis test")

library(testthat)

test_that('Hypothesis test', {
    
    d.data <- read.data('inst/extdata/multiple group/data LX.csv', type = 'data')
    d.sample <- read.data('inst/extdata/multiple group/sample LX.csv', type = 'sample')
    d.var <- read.data('inst/extdata/multiple group/var LX.csv', type = 'var')
    
    res.ht <- hypothesisTest(d.data, d.sample$Group)
    
    res.ht.2 <- hypothesisTest2(d.data, d.sample$Group)
    
    testthat::expect_equal(res.ht, res.ht.2)
    
    group <- d.sample$Group
    g_pair <- combn(levels(group), 2)
    
    
    ## Pairwise T
    res.ht.3 <- hypothesisTest2(d.data, d.sample$Group, test = 'Pairwise T')
    
    test_results <- sapply(d.data, function(x.i) {
        res.t <- apply(g_pair, 2, function(pair.i) {
            t.test(x.i[group == pair.i[1]], x.i[group == pair.i[2]])$p.value
        })
        
        names(res.t) <- apply(g_pair, 2, function(pair.i) {
            paste0('T: ', paste0(rev(pair.i), collapse = '/'))
        })
        
        res.t
    }) %>% t
    
    expect_equal(res.ht.3, test_results)
    
    ## Pairwise Wilcox
    res.ht.3 <- hypothesisTest2(d.data, d.sample$Group, test = 'Pairwise Wilcox')
    
    test_results <- sapply(d.data, function(x.i) {
        res.t <- apply(g_pair, 2, function(pair.i) {
            wilcox.test(x.i[group == pair.i[1]], x.i[group == pair.i[2]])$p.value
        })
        
        names(res.t) <- apply(g_pair, 2, function(pair.i) {
            paste0('Mann-Whitney-U: ', paste0(rev(pair.i), collapse = '/'))
        })
        
        res.t
    }) %>% t
    
    expect_equal(res.ht.3, test_results)
})
