context("Testing read data")

library(testthat)

test_that("Read csv", {
    
    d.data <- read.data('inst/extdata/multiple group/data LX.csv', type = 'data')
    d.sample <- read.data('inst/extdata/multiple group/sample LX.csv', type = 'sample')
    d.var <- read.data('inst/extdata/multiple group/var LX.csv', type = 'var')
    
    expect_equal(dim(d.data), c(24, 163))
    expect_equal(dim(d.sample), c(24, 2))
    expect_equal(dim(d.var), c(163, 1))
})