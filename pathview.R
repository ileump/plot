rm(list=ls())
suppressMessages(library(RCurl))
suppressMessages(library(rvest))
library(downloader)
library(stringr)
library(magrittr)
library(RJSONIO)
library(readr)
library(lubridate)
library(dplyr)
library(openxlsx)
# setwd("C:\\Users\\liqis\\Desktop\\1\\1")

get_map_info <- function (object, short.name = TRUE) {
    cobj = class(object)[1]
    if (cobj == "character") {
        object <- pathview:::parseKGML2(object)
        ndata = graph:::nodes(object)
    }else if (cobj == "KEGGPathway") {
        ndata = nodes(object)
    }else if (cobj == "graphNEL") {
        ndata = getKEGGnodeData(object)
    } else {
        stop("object should be either a filename, KEGGPathway, or graphNEL!")
    }
    #info_url<- purrr::map(ndata, 'link')
    nodeNames = sapply(ndata, KEGGgraph:::getName)
    nodeType = sapply(ndata, KEGGgraph:::getType)
    nodeComp = sapply(ndata, KEGGgraph:::getComponent)
    node.size = sapply(nodeComp, length)
    # ninfo = purrr::map(nodeComp, 'link')
    ninfo = purrr::map(ndata, 'link')
    grs1 = sapply(ndata, function(x) {
        grs = x@graphics
        url_link = x@link
        title = str_c(x@name,collapse = ",")
        #url_link = x[[1]]@link[1]
        c(title = title ,labels = grs@name, shape = grs@type,link_out = url_link )
    })
    coords = sapply(ndata, function(x) {
        grs = x@graphics
        coords = paste0(grs@x,",",grs@y,",",grs@width,",",grs@height)
        c(coords=coords)
    })
    grs1 = t(grs1)
    #grs2 = t(grs2)
    graphic.data = cbind(data.frame(grs1, stringsAsFactors = F), 
                         data.frame(coords))
    nd.list = cbind(kegg.names = nodeNames, type = nodeType, 
                    component = nodeComp, size = node.size)
    
    nd.list  = as.data.frame(nd.list )
    
    nd.list = cbind(nd.list, graphic.data)
    if (short.name) {
        gnames = sapply(strsplit(nd.list$labels, ", "), "[[", 
                        1)
        map.idx = nd.list$type == "map"
        gnames[map.idx] = nd.list$labels[map.idx]
        gnames[is.na(gnames)] = ""
        gnames = gsub("[.][.][.]", "", gnames)
        nd.list$labels = gnames
        nd.list$kegg.names = lapply(nd.list$kegg.names, function(x) gsub("^.*:", 
                                                                         "", x))
    }
    return(nd.list)
}

kegg.file <- "output/110/KEGG/mmu00010.xml"
info_file <- get_map_info(kegg.file)

#diff_info <-  read.table("C:\\Users\\liqis\\Desktop\\1\\1\\test.csv",header=T,sep=",",fill=TRUE,quote = '"',encoding='UTF-8',
                         # check.names=F,stringsAsFactors = FALSE)
#colnames(diff_info) <- c("Kegg_ID","Describe","FC","Pvalue")


get_html <- function(info_file,diff_info){
    result<- list()
    for (i in 1:nrow(info_file)) {
        shape <- info_file$shape[[i]]
        coords <- info_file$coords[[i]]
        title <- info_file$title[[i]]
        href <- info_file$link_out[[i]]
        if (info_file$labels[i]  %in% diff_info$Kegg_ID){
            pro_ids <- diff_info[diff_info["Kegg_ID"]==as.character(info_file$labels[i]),]
            info <- list()
            for(j in 1:nrow(pro_ids)){
                keggid <- pro_ids[j,1]
                name <- pro_ids[j,2]
                metaid2_Des = paste0("log2_FC = ",round( as.numeric(pro_ids[j,3]),2),"，Pvalue = ",pro_ids[j,4])
                map_col <- ifelse(pro_ids[j,3]>0,"up","down")
                info_a <- list(c(list(title = keggid), list(subTitle = name),
                                 list(text = metaid2_Des),list(color_subTitle = map_col)))
                info = c(info,info_a)
            }
            result =  dplyr::combine(result, list(c(shape = shape, coords = coords,
                                                    href = href, title = title,info = info)))
        }else{
            result =  dplyr::combine(result, list(c(shape = shape, coords = coords, href = href, title = title)))
        }
    }
    return(result)
}

result <- get_html(info_file,diff_info)
dd <- toJSON(result)
cat(dd, file = 'json.txt', fill = FALSE, labels = NULL, append = FALSE)



##----------------------------------------
## Bowen
##----------------------------------------
## 测试数据 data_cpd
# id   log2_fc     pvalue
# C00033 -0.4940827 0.35038712
# C00031 -0.5135232 0.07047529
# C00103 -0.7641529 0.62683875
# C00631  0.4366954 0.23831009
# C00267 -0.4766890 0.94459220
# C00221 -0.3237748 0.31589037

## 测试数据 data_gene
# id    log2_fc    pvalue
# 11674 -0.6923123 0.6850653
# 11676 -0.4056056 0.2629593
# 230163 -1.1822312 0.3310975
# 353204  0.1766734 0.4957342
# 79459 -0.1298072 0.3858522
# 110695 -1.4311203 0.5417462

kegg.file <- "output/110/KEGG/mmu00010.xml"

node_info <- node.info(kegg.file)

node_info$href <- local({
    object <- pathview:::parseKGML2(kegg.file)
    ndata <- graph:::nodes(object)
    purrr::map(ndata, 'link')
})

df_node <- lapply(names(node_info$kegg.names), function(id) {
    node_data <- purrr::map(node_info, id)
    
    coords <- NULL
    if (node_data$type == 'gene') {
        coords <- paste(node_data$x, node_data$y, node_data$x + node_data$width, node_data$y + node_data$height, sep = ',')
        
        info <- lapply(node_data$kegg.names, function(gene.id) {
            ind <- which(data_gene$id == gene.id)
            if (length(ind) > 0) {
                list(title = gene.id,
                     text = sprintf('log2_FC = %s, P-value = %s',
                                    round(data_gene[ind, 'log2_fc'], 2),
                                    round(data_gene[ind, 'pvalue'], 3)),
                     color_subTitle = ifelse(data_gene[ind, 'log2_fc'] > 0, 'up', 'down')
                )
            } else {
                NULL
            }
        }) %>% Filter(Negate(is.null), .)
        
        if (length(info) == 0) {info <- NULL}
        
    } else if (node_data$type == 'compound') {
        coords <- paste(node_data$x, node_data$y, node_data$width, sep = ',')
        
        info <- lapply(c(node_data$kegg.names, 'fasdfasf'), function(cpd.id) {
            ind <- which(data_cpd$id == cpd.id)
            if (length(ind) > 0) {
                list(title = cpd.id,
                     text = sprintf('log2_FC = %s, P-value = %s',
                                    round(data_cpd[ind, 'log2_fc'], 2),
                                    round(data_cpd[ind, 'pvalue'], 3)),
                     color_subTitle = ifelse(data_cpd[ind, 'log2_fc'] > 0, 'up', 'down')
                )
            } else {
                NULL
            }
        }) %>% Filter(Negate(is.null), .) 
        
        if (length(info) == 0) {info <- NULL}
    }
    
    
    c(id = id, node_data, coords = coords, info = list(info)) %>%
        Filter(Negate(is.null), .)
})

dd <- toJSON(df_node[1:5])
cat(dd, file = 'json.txt', fill = FALSE, labels = NULL, append = FALSE)
