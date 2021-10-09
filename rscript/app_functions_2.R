#' Read csv and excel data for data, sample and var
#'
#' @param filename file name
#' @param type type of data, either data, sample or var
#' @param ... only sheet is actually used
#'
#' @return data.frame
#' @export
#'
#' @examples
#' read.data('data.xlsx', type = 'data', sheet = 'data')
#' read.data('data.xlsx', type = 'sample', sheet = 'sample')
#' read.data('data.xlsx', type = 'var', sheet = 'var')
#' read.data('data.csv', type = 'data')
#' read.data('sample.csv', type = 'sample')
#' read.data('var.csv', type = 'var')
read.data <- function(filename, type = c('data', 'sample', 'var'), ...) {
    if(!file.exists(filename)) stop(filename, ' does not exist!')
    if(file.access(filename, mode = 4) == -1) stop(filename, ' is not readable!')
    
    filename.ext <- fs::path_ext(filename)
    dots <- list(...)
    type <- match.arg(type)
    
    ## first time read without rownames and colnames
    if (filename.ext == 'csv') {
        d0 <- read.csv(filename, header = F, strip.white = T, check.names = F)
    } else if (filename.ext == 'xlsx' | filename.ext == 'xls') {
        sheet <- dots$sheet
        if(is.null(sheet)) stop('sheet is required for excel.')
        
        d0 <- XLConnect::readWorksheetFromFile(
            filename,
            sheet = sheet,
            header = F,
            check.names = F
        )
    } else {
        stop('Only csv, xlsx and xls file formats are supported!')
    }
    
    
    ## blank rownames or is.na
    rownames.blank <- stringr::str_trim(d0[-1, 1]) %>% 
        purrr::map_lgl(~ .x == '' | is.na(.x))
    if (any(rownames.blank)) {
        stop('Blank rownames at row ', paste0(which(rownames.blank), collapse = ','))
    }
    
    ## duplicated rownames
    rownames.duplicated <- duplicated(d0[-1, 1])
    if (any(rownames.duplicated)) {
        stop(length(which(rownames.duplicated)), ' duplicated rownames: ',
             paste(d0[-1, ][rownames.duplicated, 1], collapse = ','), '.')
    }
    
    ## valid rownames
    tryCatch(
        {
            make.names(d0[-1, 1])
        },
        error = function(e) {
            stop('Unable to make syntactically valid rownames.')
        }
    )
    
    ## blank colnames or is.na
    colnames.blank <- unlist(d0[1, -1]) %>% 
        purrr::map_lgl(~ .x == '' | is.na(.x))
    if (any(colnames.blank)) {
        stop('Blank colnames at column ', paste0(which(colnames.blank), collapse = ','))
    }    
        
    ## duplicated colnames
    colnames.duplicated <- duplicated(unlist(d0[1, -1]))
    if (any(colnames.duplicated)) {
        stop(length(which(colnames.duplicated)), ' duplicated colnames: ',
             paste(d0[, -1][1, colnames.duplicated], collapse = ','), '.')
    }
    
    ## valid colnames
    tryCatch(
        {
            make.names(d0[1, -1])
        },
        error = function(e) {
            stop('Unable to make syntactically valid colnames.')
        }
    )
    
    rm(d0)
    
    ## second time read with rownames and colnames
    if (filename.ext == 'csv') {
        d1 <- read.csv(filename, header = T, row.names = 1,
                       strip.white = T, check.names = F)
    } else if (filename.ext == 'xlsx' | filename.ext == 'xls') {
        sheet <- dots$sheet
        if(is.null(sheet)) stop('sheet is required for excel.')
        
        d1 <- XLConnect::readWorksheetFromFile(
            filename,
            sheet = sheet,
            header = T,
            rownames = 1,
            check.names = F
        )
    } else {
        stop('Only csv, xlsx and xls file formats are supported!')
    }
    
    if (type == 'data') {
        ## error message for NA
        if (any(is.na(d1))) {
            stop("Missing value exists, ",
                    paste(colnames(d1)[apply(d1, 2, function(x) {any(is.na(x))})],
                          collapse = ',')
            )
        }
        
        ## check column type
        if (!all(sapply(d1, is.numeric))) {
            ind_cols <- which(!sapply(d1, is.numeric))
            stop(
                "Non-numeric colum exists, ",
                paste(paste(colnames(d1)[ind_cols], collapse = ','),
                      'is not numeric.')
            )
        }
        
        if (any(apply(d1, 2, sd) == 0)) {
            stop(
                "Column with 0 variation exists, ",
                paste(colnames(d1)[apply(d1, 2, sd) == 0],
                      collapse = ',')
            )
        }
    }
    
    if (type == 'sample') {
        testthat::expect_true('Group' %in% colnames(d1), info = 'There exists no column named Group.')
        if ('Order' %in% colnames(d1)) {
            d1$Group <- factor(
                d1$Group,
                levels = dplyr::arrange(d1[!duplicated(d1$Group), ], Order)[, 'Group']
            )
        } else {
            d1$Group <- factor(d1$Group)
        }
    }
    
    d1
}

