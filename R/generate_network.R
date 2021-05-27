# generate cell populations for simulation
#' @export
create_cellpops <- function(npops=1, edge_prob=0.2, min_mod=10, max_mod=20) {
    n_cellpops <- 1:npops
    cellpop_names <- paste0('CellPop', n_cellpops, '_')
    
    cellpop_meta <- construct_populations(n_cellpops, cellpop_names, edge_prob, min_mod, max_mod)
    branch_meta <- create_branches(n_cellpops, cellpop_names, cellpop_meta$module_info)
    
    start_branch <- bblego_start('B', type = "simple", num_modules = 2)
    start_branch$expression_patterns$time = 25
    backbone <- bblego(start_branch, branch_meta, cellpop_meta)
    
    return (list(backbone=backbone))
}

construct_populations <- function(n_cellpops, cellpop_names, edge_prob, min_mod, max_mod) {
    # For every cell population specified generate a random module network unique to that cell pop.
    cellpop_meta <- c()
    modules_range <- min_mod:max_mod
    
    for (i in n_cellpops) {
        # create cell population name
        cellpop_name <- cellpop_names[i]

        # select number of modules in cell pop and generate graph
        num_modules <- sample(modules_range, 1)
        generated_net <- erdos.renyi.game(num_modules, edge_prob, directed=TRUE, loops=TRUE)

        # create cell population metadata
        module_info <- create_module_info(generated_net, cellpop_name)
        expression_patterns <- create_expression_pattern(module_info, cellpop_name)
        module_network <- create_module_network(generated_net, module_info, cellpop_name)

        # adding current cell pop meta to cellpop_meta vector
        meta <- list(module_info=module_info, module_network=module_network, expression_patterns=expression_patterns)
        cellpop_meta <- c(cellpop_meta, list(meta))
    }
    
    cellpop_meta <- list(module_info = map_df(cellpop_meta, "module_info"), 
                        module_network = map_df(cellpop_meta, "module_network"), 
                        expression_patterns = map_df(cellpop_meta, "expression_patterns"))
    
    return (cellpop_meta)
}

create_module_info <- function(cellpop_network, cellpop_name) {
    # generate the module info 
    module_info <- list()
    module_info$module_id <- as.numeric(V(cellpop_network))
    module_info$in_degree <- degree(cellpop_network, mode='in')
    module_info$out_degree <- degree(cellpop_network, mode='out')
    
    module_info <- module_info %>% 
                    as_tibble %>%
                    mutate(module_id = paste0(cellpop_name, module_id),
                           basal = ifelse(in_degree == 0, 1, 0),
                           burn = FALSE,
                           independence = 1)
    
    return (module_info)
}

create_expression_pattern <- function(module_info, cellpop_name) {
    # creating an expression patterns
    sphase <- paste0("s", cellpop_name)
    ephase <- paste0("e", cellpop_name)
    modules_in_state <- paste("+", module_info$module_id, sep='', collapse=',')

    expression_patterns <- tribble(
      ~from, ~to, ~module_progression, ~start, ~burn, ~time,
      sphase, ephase, modules_in_state, FALSE, FALSE, 100,
    )
    
    return (expression_patterns)
}

create_module_network <- function(cellpop_net, module_info, cellpop_name) {
    # creating the module networks metadata
    module_network <- igraph::as_edgelist(cellpop_net) %>%
                                as.data.frame %>%
                                rename(from = V1, to = V2) %>%
                                mutate(to = paste0(cellpop_name, to),
                                       from = paste0(cellpop_name, from), 
                                       effect = 1L, 
                                       strength = sample(1:10, gsize(cellpop_net), replace = TRUE),
                                       hill = 2) %>%
                                do(randomize_module_effect(.data, module_info)) %>%
                                as_tibble
    
    return (module_network)
}

randomize_module_effect <- function(module_expr, modules_info) {
    # randomly assign whether in degree positively or negatively regulate
    # this module. Unless there is 1 or 0 in degree. In this case level the effect
    # as positive regulation.
    
    nodes_in_degree <- modules_info$in_degree
    names(nodes_in_degree) <- modules_info$module_id
    
    for (i in 1:nrow(module_expr)) {
        module <- module_expr[i, ]$to
        in_degree <- nodes_in_degree[[module]]
        
        if (in_degree == 1){
            next()
        }

        nodes_in_degree[[module]] <- in_degree - 1 
        module_expr[i, ]$effect <- sample(c(1L, -1L), size=1)   
    }
    
    return (module_expr)
}


create_branches <- function(n_cellpops, cellpop_names, module_info) {
    network <- data.frame()
    info <- data.frame(module_id=c("B1"), basal=c(0), burn=c(TRUE), independence=c(1))
    expr_pat <- data.frame(from=c("sB","sBmid"),to=c("sBmid","sBend"),
                          module_progression=c('+B1', paste0('+B1,',paste('+B', n_cellpops + 1, sep='', collapse= ','))),
                          start=c(FALSE),burn=c(TRUE),time=c(25))

    for (i in n_cellpops) {
        cellpop_name <- cellpop_names[i]
        branch <- paste0('B', i+1)
        start <- paste0('s', cellpop_name)
        cellpop_module <- module_info %>% filter(grepl(cellpop_name, module_id))  %>% pull(module_id)
        
        info <- create_branches_info(info, branch)
        expr_pat <- create_branches_expr_pat(expr_pat, branch, start)
        network <- create_branches_network(network, branch, n_cellpops, cellpop_module)
    }

    branch_meta <- list(module_info = info, module_network = network, expression_patterns = expr_pat)
    return (branch_meta)
}

find_largest_out_degree <- function(cellpop_name, module_info) {
    cell_info <- module_info %>% filter(grepl(cellpop_name, module_id)) %>% pull(module_id)
    cell_out_degree <- cell_info %>% pull(out_degree) %>% max()
    max_modules <- cell_info %>% filter(out_degree == cell_out_degree) %>% select(module_id)
    
    if (nrow(max_modules) == 1) {
        cellpop_module <- max_modules$module_id
    } else {
        cellpop_module <- max_modules$module_id[1]
    }
    
    return (cellpop_module)
}

create_branches_info <- function(info, branch) {
    info_row <- data.frame(module_id=c(branch), basal=c(0), burn=c(TRUE), independence=c(1))
    info <- rbind(info, info_row)
    
    return (info)
}

create_branches_expr_pat <- function(expr_pat, branch, start) {
    
    expr_pat_row <- data.frame(from=c("sBend"), to=c(start), 
                              module_progression=c(paste0('+',branch)),
                              start=c(FALSE), burn=c(TRUE), time=c(25))
    
    expr_pat <- rbind(expr_pat, expr_pat_row)
    
    return (expr_pat)
}

create_branches_network <- function(network, branch, n_cellpops, cellpop_module) {
    # adding branch metadata to module network
    burn_row <- data.frame(from=c("B1"), to=c(branch), effect=c(1L), strength=c(1), hill=c(2))
    cellpop_row <- data.frame(from=c(branch), to=c(cellpop_module), effect=c(1L), strength=c(1), hill=c(2))
    
    network <- rbind(network, burn_row)
    network <- rbind(network, cellpop_row)
    
    # adding inhibitory relations between different branches in the module network
    for (i in n_cellpops) {
        other_branch <- paste0('B', i + 1)
        
        # we don't want a given branch to inhibit itself.
        if (branch == other_branch) {
            next()
        }
        
        network_row <- data.frame(from=c(branch),to=c(other_branch),effect=c(-1L),strength=c(1),hill=c(2))
        network <- rbind(network, network_row)
    }
    
    return (network)
}