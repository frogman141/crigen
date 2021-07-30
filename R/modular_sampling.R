#' @export
modular_sample_net <- function(net, num_nodes) {  
    # adding nodes to the extracted sub network until i read users specified size
    for (i in 1:num_nodes) {
        
        if (i == 1) {
            # first node is the seed node with is randomly sampled forom TF network
            sampled_nodes <- c(sample(V(net)$name, size=1))
        } else {
            node_to_add <- find_node_to_add(net, sampled_nodes)
            sampled_nodes <- c(sampled_nodes, node_to_add)
        }
    
    }
    
    return (sampled_nodes)
}

#' @export
get_nieghboring_nodes <- function(net, sampled_nodes, mode='out') {
    neighboring_nodes <- c()
    neighbors <- igraph::adjacent_vertices(net, sampled_nodes, mode='out')

    for (node in sampled_nodes) {
        node_neighbors <- neighbors[[node]]
        num_of_neighbors <- length(node_neighbors)

        for (i in 1:num_of_neighbors) {
            neighbor <- node_neighbors[i]$name
            neighboring_nodes <- c(neighboring_nodes, neighbor)
        }
    }
    
    return (neighboring_nodes)
}

find_node_to_add <- function(net, sampled_nodes) {
    neighboring_nodes <- get_nieghboring_nodes(net, sampled_nodes)
    modularity_df <- calculate_modularity(net, sampled_nodes, neighboring_nodes)
    node_to_add <- get_node_with_max_modularity(modularity_df)
    
    return (node_to_add)
}

calculate_modularity <- function(net, sampled_nodes, neighboring_nodes) {
    q <- c()
    nodes <- V(net)$name
    
    for (neigbor in neighboring_nodes) {
        proposed_nodes <- c(sampled_nodes, neigbor)
        membership_ids <- get_membership(proposed_nodes, nodes)

        proposed_q <- igraph::modularity(net, membership_ids)
        q <- c(q, proposed_q)
    }
    
    modularity_df <- data.frame(nodes=neighboring_nodes, modularity=as.double(q), stringsAsFactors=FALSE)
    return (modularity_df)
}

get_node_with_max_modularity <- function(modularity_df) {
    node_to_add <- modularity_df %>% filter(modularity == max(modularity)) %>% pull(nodes)
    
    # if there are multiple nodes with the same modularity randomly select one.
    if (length(node_to_add) > 1){
        node_to_add <- sample(node_to_add, size=1)
    }
    
    return (node_to_add)
}

get_membership <- function(proposed_nodes, net_nodes) {
    membership_id <- c()
    
    for (node in net_nodes){
        if (node %in% proposed_nodes){
            membership_id <- c(membership_id, 2)
        } else {
            membership_id <- c(membership_id, 1)
        }
    }
    
    return (membership_id)
}