#' @export
extract_housekeeping_net <- function(bio_net, num_hks) {
    # extract housekeeping network using breadth first sampling
    seed_node <- sample(names(V(bio_net)), size=1)

    hk_names <- bio_net %>% 
                    igraph::bfs(seed_node, neimode = "all") %>% 
                    `[[`("order") %>% 
                    `[`(seq_len(num_hks)) %>% 
                    names()

    hk_net <- bio_net %>%
              igraph::induced_subgraph(hk_names) %>%
              igraph::as_data_frame() %>%
              filter(.data$from != .data$to) %>% # removing loops
              igraph::graph_from_data_frame(directed=TRUE)
    
    return (hk_net)
}