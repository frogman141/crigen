#' @export
create_network <- function(num_tfs, num_egenes, num_hks, cache_dir="~/.cache/crigen", bio_net_name=sample(fantom5$dataset_id, size=1)) {
    
    # if no biological network provided select one
    bio_net <- get_network(bio_net_name, cache_dir)
       
    # sample grn and housekeeping network from a biological network
    message("Creating GRN...")
    grn <- create_grn(bio_net, num_tfs, num_egenes) %>% igraph::as_data_frame()
    message("GRN Created...")
    
    message("Creating HK Network...")
    hk_net <- create_hk_net(bio_net, num_hks) %>% igraph::as_data_frame()
    message("HK Network Created...")
    
    # merging sampled GRN and housekeeping network and extracting metadata
    feature_network <- bind_rows(grn, hk_net) %>%
                                mutate(valid_edge=ifelse(grepl('EGene', from) & grepl('EGene', to), TRUE,
                                                  ifelse(grepl('TF', from) | grepl('HK', from), TRUE, FALSE))) %>%
                                filter(valid_edge) %>% 
                                select(-weight)
    
    
    feature_info <- data.frame(feature_id=unique(c(feature_network$from, feature_network$to))) %>%
                        mutate(is_tf=ifelse(grepl('TF', feature_id), TRUE, FALSE),
                               is_hk=ifelse(grepl('HK', feature_id), TRUE, FALSE),
                               is_egene=ifelse(grepl('EGene', feature_id), TRUE, FALSE),
                               burn=TRUE)
    
    simulation_system <- lst(cache_dir)
    return (lst(feature_network, feature_info, simulation_system))
}

get_network <- function(bio_net_name, cache_dir) {
    url <- fantom5 %>% filter(dataset_id == bio_net_name) %>% pull(links.download)
    adj_mat <- download_cacheable_file(url, cache_dir, TRUE)
    
    bio_net <- adj_mat %>%
                    Matrix::summary() %>% 
                    as.data.frame() %>% 
                    transmute(
                      i = rownames(adj_mat)[.data$i], 
                      j = colnames(adj_mat)[.data$j],
                      weight = .data$x) %>%
                    igraph::graph_from_data_frame(vertices = colnames(adj_mat))
    
    return (bio_net)
}

create_grn <- function(bio_net, num_tfs, num_egenes) {
    grn <- NA
    pass_qc <- FALSE
    
    # keep sampling TFs and E-Genes until we have created a network with
    # the desired number of Transcription Factors and E-Genes
    while (isFALSE(pass_qc)) {
        # sample biological network for TFs and E-Genes
        tfs <- sample_tfs(bio_net, num_tfs)
        egenes <- sample_egenes(bio_net, tfs, num_egenes)
        grn <- get_grn(tfs, egenes, bio_net)
        
        # rename genes so there anonymized
        genes_in_grn <- c(tfs, egenes)
        tf_mapper <- rename_genes(tfs, gene_type='TFs')
        egene_mapper <- rename_genes(egenes, gene_type='EGenes')
        
        gene_mapper <- append(tf_mapper, egene_mapper)
        V(grn)$name <- ifelse(V(grn)$name %in% genes_in_grn, gene_mapper[V(grn)$name], V(grn)$name)
        
        # check grn meets user parameters
        pass_qc <- grn_qc(grn, num_tfs, num_egenes)
    }
    
    return (grn)
}


create_hk_net <- function(bio_net, num_hks) {
    hk_net <- extract_housekeeping_net(bio_net, num_hks)
    hk_mapper <- rename_genes(V(hk_net)$name, gene_type='HK')
    
    # assigning new names to housekeeping genes
    V(hk_net)$name <- ifelse(V(hk_net)$name %in% names(hk_mapper), hk_mapper[V(hk_net)$name], V(hk_net)$name)
    return (hk_net)
}

get_grn <- function(tfs, egenes, bio_net) {
    
    grn <- bio_net %>% igraph::induced_subgraph(c(tfs, egenes))
    
    # remove any edges to my current tfs
    grn_df <- grn %>% igraph::as_data_frame() %>% select(from, to)
    
    # identify egenes that have no in degree edges
    missing_targets <- setdiff(egenes, grn_df$to)
    
    if (length(missing_targets) > 0) {
        # if there are egenes with missing edges reassign using degree distribution
        deg <- igraph::degree(grn)
        new_regs <- sample(names(deg), length(missing_targets), prob = deg, replace = TRUE)
        grn_df <- bind_rows(grn_df, tibble(from = new_regs, to = missing_targets))
    }
    
    grn <- grn_df %>% igraph::graph_from_data_frame(directed = TRUE)
    return (grn)
}

grn_qc <- function(grn, num_tfs, num_egenes){
    # check grn meets user parameters
    pass_qc <- FALSE
    genes_in_grn <- V(grn)$name
    
    tf_count <- length(genes_in_grn[grepl('TF', genes_in_grn)])
    egene_count <- length(genes_in_grn[grepl('EGene', genes_in_grn)])
    
    if (tf_count == num_tfs & egene_count == num_egenes){
        pass_qc <- TRUE
    }
    
    return (pass_qc)
}

rename_genes <- function(genes, gene_type="HK") {
    # rename all of the nodes sampled from bio_net to generic labels
    num_of_genes <- length(genes)
    
    if (gene_type == 'TFs') {
        gene_mapper <- set_names(paste0('TF', 1:num_of_genes), genes)
    } else if (gene_type == 'EGenes') {
        gene_mapper <- set_names(paste0('EGene', 1:num_of_genes), genes)
    } else {
        gene_mapper <- set_names(paste0('HK', 1:num_of_genes), genes)
    }
    
    return (gene_mapper)
}