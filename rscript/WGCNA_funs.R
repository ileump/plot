##-----------------------------------------------------------------------
## MEGENA network
##-----------------------------------------------------------------------
require(igraph)
megenaGraph <- function(x, v, MEGENA.output, main, mark.module = F, vlabel = F,
                        plotTop = 0) {
    # Plot network of modified MEGENA output
    #
    # Args:
    #   x: A data frame with columns: row, col, zScoreDiff, Class, color. row, col,
    #      zScoreDiff are output from MEGENA::calculate.PFN. Class is modified from
    #      Classes, append '-' or '+' to '-/-' and '+/+' to indicate decrease or
    #      increase of correlation, determined based on zScoreDiff.
    #   v: A data frame of vertice information
    #   MEGENA.output: Ouptut of do.MEGENA
    #   mark.module: mark modules in the network
    #   vlabel: label all the nodes
    #
    # Returns:
    #   Plot network
    
    ## construct network
    g <- igraph::graph.data.frame(x, directed = FALSE)
    
    V(g)$color <- v[attr(V(g), 'names'), 'color']
    
    if ('shape' %in% colnames(v))
        V(g)$shape <- v[attr(V(g), 'names'), 'shape']
    
    
    if (vlabel) {
        V(g)$cex <- 0.6
        
        # if (!is.na(MEGENA.output$node.summary))
        # V(g)$size <- as.numeric(MEGENA.output$node.summary[attr(V(g), 'names'), 'node.degree']) ** 0.3 * 3
        V(g)$size <- degree(g) ** 0.3 * 3
    } else {
        ## Vertice label size is 0 if degree < 4
        V(g)$cex <- 0.6
        
        # vertice size is proportional to sqrt of degree
        # if (!is.na(MEGENA.output$node.summary))
        # V(g)$size <- as.numeric(MEGENA.output$node.summary[attr(V(g), 'names'), 'node.degree']) ** 0.3 * 1.5
        V(g)$size <- degree(g) ** 0.3 * 1.5
    }
    ## Vertice size is 1, 1/2, 1/3 for hub at scale 1, 2, 3
    if (!is.na(MEGENA.output$hub.output)) {
        for(i in 1:length(MEGENA.output$hub.output$hub.list)) {
            V(g)$cex <- ifelse(
                attr(V(g), 'names') %in% MEGENA.output$hub.output$hub.list[[i]],
                1 / i**0.2, V(g)$cex) * 1.2
        }
    }
    
    
    
    # In MEGENA.output$module.output$modules, the first one is parent module which
    # includes all the nodes
    groupList <- MEGENA.output$module.output$modules[-1]
    groupColours <- c(rainbow(length(groupList), alpha = 0.1))
    
    if (length(groupList) == 0)
        mark.module <- F
    # apply layout
    set.seed(5)
    l <- layout_with_dh(g)
    l <- norm_coords(l, ymin = -1, ymax = 1, xmin = -1, xmax = 1)
    
    # replace name by label
    if ('label' %in% colnames(v))
        V(g)$name <- v[attr(V(g), 'names'), 'label']
    
    # par(mar = c(5.1, 4.1, 4.1, 10.1))
    vertex.label <- ''
    if (vlabel)
        vertex.label <- V(g)$name
    
    if (mark.module) {
        plot(g,
             mark.groups=groupList, # Mark the groups
             mark.col= groupColours,
             mark.border = NA,
             edge.width = 0.6,
             vertex.size = V(g)$size,
             #vertex.color = vertex.color,
             vertex.label = vertex.label,
             vertex.label.color = 'black',
             vertex.label.cex = V(g)$cex,
             vertex.label.dist = 0.3,
             rescale = TRUE, layout = l, axes = FALSE,
             main = main)
        
    } else {
        plot(g,
             # mark.groups=groupList, # Mark the groups
             # mark.col= groupColours,
             # mark.border = NA,
             edge.width = 0.6,
             vertex.size = V(g)$size,
             #vertex.color = vertex.color,
             vertex.label = vertex.label,
             vertex.label.color = 'black',
             vertex.label.cex = V(g)$cex,
             vertex.label.dist = 0.3,
             rescale = TRUE, layout = l, axes = FALSE,
             main = main)
    }
    
    ## vertice legend label
    v.legend <- data.frame(
        label = unique(v[, 'class']),
        color = unique(v[, 'color'])
    ) %>% dplyr::filter(
        label %in% v[c(pfn_res$row, pfn_res$col), 'class']
    )
    v.legend <- v.legend[order(v.legend$label), ]
    
    ## Vertices
    legend(x=1.2, y=1.2,
           y.intersp = 1,
           v.legend$label, pch=21,
           pt.bg = v.legend$color,
           col="#777777", pt.cex=2, cex=1, bty="n", ncol=1)
    
    if (plotTop > 0) {
        degreeStat <- apply(g[], 2, function(x) {sum(x > 0)})
        degreeStat <- data.frame(
            name = names(degreeStat),
            degree = degreeStat
        )
        degreeStat <- degreeStat[order(degreeStat$degree, decreasing = T), ]
        
        for (node.i in degreeStat[1:min(plotTop, nrow(degreeStat)), 'name']) {
            testthat::expect_equal(
                colnames(g[]), V(g)$name
            )
            ## diagnal entry is 0
            g.sub <- induced_subgraph(
                g, g[node.i] > 0 | colnames(g[]) == node.i
            )
            ## subgraph layout
            l.sub <- l[g[node.i] > 0 | colnames(g[]) == node.i, ]
            plot(g.sub,
                 # mark.groups=groupList, # Mark the groups
                 # mark.col= groupColours,
                 # mark.border = NA,
                 edge.width = 0.6,
                 # vertex.size = V(g)$size,
                 #vertex.color = vertex.color,
                 # vertex.label = vertex.label,
                 vertex.label.color = 'black',
                 vertex.label.cex = V(g.sub)$cex,
                 vertex.label.dist = 0.3,
                 rescale = TRUE, layout = l.sub, axes = FALSE,
                 main = paste(main, node.i))
            
            ## vertice legend label
            v.legend <- data.frame(
                label = unique(v[, 'class']),
                color = unique(v[, 'color'])
            ) %>% dplyr::filter(
                label %in% v[V(g.sub)$name, 'class']
            )
            v.legend <- v.legend[order(v.legend$label), ]
            
            ## Vertices
            legend(x=1.2, y=1.2,
                   y.intersp = 1,
                   v.legend$label, pch=21,
                   pt.bg = v.legend$color,
                   col="#777777", pt.cex=2, cex=1, bty="n", ncol=1)
        }
    }
}

moduleStats <- function(g, MEGENA.output, v, module) {
    ## g is igraph
    ## MEGENA.output is output from do.MEGENA
    ## v is output from getAttrs
    ## module: integer, module ID
    
    ## calculate degree of each node
    node.degree <- apply(g[], 2, function(x) {sum(x > 0)})
    ## sort by degree
    node.degree <- node.degree[order(node.degree, decreasing = T)]
    ## edge list to data.frame
    ## format: node1 node2
    edges <- stringr::str_split(attr(E(g), 'vnames'), '\\|') %>%
        do.call('rbind', .) %>% `colnames<-`(c('node1', 'node2')) %>%
        as.data.frame
    
    sapply(names(node.degree), function(node.i) {
        
        ## all the connected node to node.i
        connected.nodes <- c(edges$node1[edges$node2 == node.i],
                     edges$node2[edges$node1 == node.i])
        
        ## assert that the number of connected nodes is equal to
        ## the degree
        # testthat::expect_equal(
        #     length(nodes.i), node.degree[node.degree$gene == node.i, 'degree']
        # )
        
        stats.i <- lipidStats(connected.nodes, v)
        stats.i[, 'node'] <- node.i
        out <- tidyr::spread(
            stats.i, key = 'class', value = 'V1'
        )
        out <- cbind(module = module, 
                     degree = node.degree[node.i],
                     out)
    }) %>% t
    
}

lipidStats <- function(x, v) {
    ## output: average number of carbon chain length and double bonds
    if (is.data.frame(v) & (!is.data.table(v)))
        v <- data.table(v)
    
    ## we want all the classes so that we get a table in the end
    out0 <- data.frame(
        class = v[, unique(class[!is.na(class)])],
        V1 = ''
    )
    
    out <- v[name %in% x, sprintf(
        '%d, %.1f (%.1f), %.1f (%.1f)',
        .N,
        round(mean(carbon, na.rm = T), 1), round(sd(carbon, na.rm = T), 1),
        round(mean(dbond, na.rm = T), 1), round(sd(dbond, na.rm = T), 1)
    ),
    by = class
    ]
    
    ## rbind with out0
    out <- rbind(out, subset(out0, !(class %in% out$class)))
    
    ## sort by alphabetical order
    out <- out[order(out$class), ]
    
    out
}
