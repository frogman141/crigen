#' @export
sample_egenes <- function(bio_net, tfs, num_egenes) {
    downstream_egenes <- get_nieghboring_nodes(bio_net, tfs, mode='out')
    page_rank <- run_page_rank(bio_net, c(tfs, downstream_egenes))
    egenes <- sample_page_rank_results(page_rank, tfs, num_egenes)
    
    return (egenes)
}

run_page_rank <- function(bio_net, gene_universe) {
    
    page_rank <- bio_net %>%
                    igraph::page_rank(directed = TRUE,
                                      weights = igraph::E(bio_net)$weight,
                                      damping=0.1)
    
    return (page_rank)
}

sample_page_rank_results <- function(page_rank, tfs, num_egenes) {
    
    sampled_egenes <- enframe(page_rank$vector, "feature_id", "score") %>% 
                          mutate(score = .data$score + runif(n(), 0, 1e-15)) %>% 
                          filter(!.data$feature_id %in% tfs) %>% 
                          sample_n(num_egenes, weight = .data$score) %>% 
                          pull(.data$feature_id) %>% unique()
    
    return (sampled_egenes)
}