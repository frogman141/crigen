#' @export
sample_tfs <- function(bio_net, num_tfs) {
    # Extract Transcription Factor Interaction Network (TFI Net) from Biological Network. 
    # Then Run Modular Sampling across the TFI net to create an in silico TFI Net. 
    tf_net <- get_tf_net_from_bio_net(bio_net)
    
    # need to add test to check the num of tfs requested is lower than number of tfs in bio net provided
    assert_that(
        num_tfs <= gorder(tf_net),
        msg = paste0("Number of Transcription Factors in Biological Network (", nrow(realnet), ") is not large enough (>= ", nrow(model$feature_info), ")")
    )
    
    
    tfs_nodes <- modular_sample_net(tf_net, num_tfs)
    return (tfs_nodes)
}

get_tf_net_from_bio_net <- function(bio_net) {
    # identify nodes with greater than 1 out degree.
    nodes <- bio_net %>% 
                igraph::degree(mode='out') %>%
                as.data.frame %>%
                filter(. > 0) %>%
                rownames
    
    # extract the TF network and convert to undirected graph
    # this makes it easier to run modular sampling
    tf_net <- bio_net %>%
                igraph::induced_subgraph(nodes) %>%
                igraph::as.undirected()
    
    return (tf_net)
}