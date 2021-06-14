library(igraph)
library(nem)
library(data.table)

##TODO:

#1. Add names DONE
#2. Add replicates  DONE
#3. Add p-value density matrices (first without dynamic errors) DONE
#4. Add sampling from KEGG to randomPhi DONE
#5. Add dynamic errors (alpha/beta) with increasing path length DONE
#6. Add dynamic errors as above, but by path length in original KEGG graph  DONE
#7. Add deterministic E-gene attachment DONE
#8. Add leaf-biased E-gene attachment DONE

### Simulate data for NEMs ###

## BUM helper functions used in nem package, but not exported:

# cdf
pbum <- function(x,a,lambda){					
  lambda[1]*x + lambda[2]*pbeta(x,a[1],1) + lambda[3]*pbeta(x,1,a[2])
}

# density function
dbum <- function(x,a,lambda){		
  lambda[1] +  lambda[2]*dbeta(x,a[1],1) + lambda[3]*dbeta(x,1,a[2])
}

# cdf of alternative distribution
bum.palt <- function(x,a,lambda){
  pihat = dbum(1,a,lambda) # uniform component
  (pbum(x,a,lambda) - pihat*x)/(1-pihat)
}

# quantile function
bum.qalt <- function(p,a,lambda,nbisect=20){				
  n <- length(p)
  mid <- rep(0,n)
  top <- rep(1,n)
  bot <- rep(0,n)
  gohigher <- rep(FALSE,n)
  for (j in 1:nbisect){
    mid <- (top+bot)/2
    gohigher <- (bum.palt(mid,a,lambda)<p)		
    bot[gohigher] <- mid[gohigher]
    top[!gohigher] <- mid[!gohigher]
  }
  return(mid)
}

# draw n random samples from alternative distribution
bum.ralt <- function(n,a,lambda){
  u <- runif(n)
  return(bum.qalt(u,a,lambda))
}

#quick plotting wrapper

plotAdj <- function(adj_matrix, no_self_loops = TRUE, colour_signs = TRUE, ...){
  
  if(no_self_loops == TRUE){
    diag(adj_matrix) <- 0
  }
  
  g <- graph_from_adjacency_matrix(adj_matrix, weighted = TRUE)
  if(colour_signs == TRUE){
    weights <- E(g)$weight
    E(g)$colours <- ifelse(weights == 1, "darkgray", "red")
    plot(g, edge.color = E(g)$colours, ...)
  }else{
    plot(g, ...)
 }
}

#preset here is the cancer networks KEGG xml

#Functions for randomPhi below:

simple_input_sample  <- function(n_nodes, input_graph,
                                 sample_transitive_closed_input = TRUE,
                                 return_input_distances = FALSE){
  
  stopifnot(!is.null(input_graph))
  subgraph_nodes <- sample(V(input_graph), n_nodes)
  
  if(sample_transitive_closed_input == TRUE){
    
    input_distances <- distances(input_graph, v = subgraph_nodes, 
                                 to = subgraph_nodes, mode = 'out')
    
    trans_subgraph_mat <- (!is.infinite(input_distances))*1
    phi <- suppressWarnings(transitive.reduction(trans_subgraph_mat))
    
    if(return_input_distances == TRUE){
      return(list(phi, input_distances))
    }
    
  }else{
    
    out <- as_adj(induced_subgraph(input_graph, subgraph_nodes), sparse = FALSE)
    
    if(return_input_distances == TRUE){
      input_distances <- distances(input_graph, v = subgraph_nodes, 
                                   to = subgraph_nodes, mode = 'out')
      out <- list(out, input_distances)
    }
    
    return(out)
  }
  
}

seed_neighbours <- function(n_nodes, input_graph, 
                            sample_transitive_closed_input = TRUE,
                            return_input_distances = FALSE){
  
  stopifnot(!is.null(input_graph))
  max_dist <- n_nodes - 1
  seed_node <- sample(V(input_graph), 1)
  
  to_v <- V(input_graph)
  to_v <- to_v[to_v != seed_node]
  
  seed_dists <- distances(input_graph, seed_node, to = to_v)
  seed_dists[seed_dists > max_dist] <- max_dist + 1
  seed_weights <- 1/(seed_dists)
  
  subgraph_nodes <- c(seed_node, 
                      sample(to_v, size = n_nodes - 1, prob = seed_weights))
  
  if(sample_transitive_closed_input == TRUE){
    
    input_distances <- distances(input_graph, v = subgraph_nodes, 
                                 to = subgraph_nodes, mode = 'out')
    
    trans_subgraph_mat <- (!is.infinite(input_distances))*1
    out <- suppressWarnings(transitive.reduction(trans_subgraph_mat))
    
    if(return_input_distances == TRUE){
      out <- list(out, input_distances)
    }
    
  }else{
    
    out <- as_adj(induced_subgraph(input_graph, subgraph_nodes), sparse = FALSE)
    
    if(return_input_distances == TRUE){
      input_distances <- distances(input_graph, v = subgraph_nodes, 
                                   to = subgraph_nodes, mode = 'out')
      out <- list(out, input_distances)
    }
  }
  
  return(out)
}

walk_sample <- function(n_nodes, input_graph,
                        sample_transitive_closed_input = TRUE,
                        reset_probability = 0.1,
                        return_input_distances = FALSE){
  
  stopifnot(!is.null(input_graph))
  nodes_sampled <- 1
  subgraph_nodes <- names(sample(V(input_graph), 1))
  current_node <- subgraph_nodes
  
  while(nodes_sampled != n_nodes){
    
    new_node <- ifelse(runif(1) > reset_probability, 
                       tail(names(igraph::random_walk(graph = input_graph, 
                                                      start = current_node, 
                                                      steps = 2, mode = 'all')), 1),
                       sample(names(V(input_graph)), 1))
    
    subgraph_nodes <- unique(c(new_node, subgraph_nodes))
    current_node <- sample(subgraph_nodes, 1)
    nodes_sampled <- length(subgraph_nodes)
  }
  
  if(sample_transitive_closed_input == TRUE){
    
    input_distances <- distances(input_graph, v = subgraph_nodes, 
                                 to = subgraph_nodes, mode = 'out')
    trans_subgraph_mat <- (!is.infinite(input_distances))*1
    out <- suppressWarnings(transitive.reduction(trans_subgraph_mat))
    
    if(return_input_distances == TRUE){
      out <- list(out, input_distances)
    }
    
  }else{
    out <- as_adj(induced_subgraph(input_graph, subgraph_nodes), sparse = FALSE)
    
    if(return_input_distances == TRUE){
      input_distances <- distances(input_graph, v = subgraph_nodes, 
                                   to = subgraph_nodes, mode = 'out')
      out <- list(out, input_distances)
      
    }
    
  }
  
  return(out)
  
}

transitive_closure <- function(g, self_loops = TRUE){
    
    ##This is just a faster version of the transitive.closure function from the nem package.
    ##It's faster because, under-the-hood, igraph is written in C, so even though this code
    ##strictly isn't that smart a way of doing it, the breadth-first-search in igraph is very
    ##fast. 
    
    if(class(g) != 'igraph'){
        ig <- graph_from_adjacency_matrix(g)
    }else{
        ig <- g
    }
    
    reset_names <- FALSE
    if(is.null(V(ig)$name)){
        reset_names <- TRUE
        V(ig)$name <- paste0("S", 1:length(V(ig)))
    }
    
    if(is.directed(ig)){
        reachability <- lapply(V(ig), function(v) subcomponent(ig, v, "out")$name)
    }else{
        reachability <- lapply(V(ig), function(v) subcomponent(ig, v, "all")$name)
    }
    
    edges_to_add <- c(t(stack(reachability)[,c(2,1)]))
    ig <- add_edges(ig, edges_to_add)
    if(self_loops == TRUE){
        ig <- simplify(ig, remove.loops = FALSE)
    }else{
        ig <- simplify(ig)
    }
    
                               
    if(reset_names == TRUE){
        ig <- delete_vertex_attr(ig, "name")
    }
    return(ig)
                               
}

signed_transitive_closure <- function(adj_mat){
    
    weighted_graph <- graph_from_adjacency_matrix(adj_mat, weighted = TRUE)
    abs_trans_close <- as_adj(transitive_closure(abs(adj_mat)), sparse = FALSE)
    pairs_to_check <- get.edgelist(graph_from_adjacency_matrix(abs(adj_mat) != abs_trans_close))
    
    if(nrow(pairs_to_check) == 0){
        return(adj_mat)
    }
    
    shortcuts <- unlist(lapply(1:nrow(pairs_to_check), 
                        function(i){
                          shortest_paths(weighted_graph, from = pairs_to_check[i,1], 
                                         to = pairs_to_check[i,2],
                                         mode = "out", weights = NA, output = "epath")$epath
                          }), recursive = FALSE)
    shortcut_signs <- unlist(lapply(shortcuts, function(x) prod(E(weighted_graph)$weight[x])))
    edges_to_add <- c(t(pairs_to_check))
    weighted_graph <- add_edges(graph = weighted_graph, edges = edges_to_add, weight = shortcut_signs)
    out <- as_adj(weighted_graph, sparse = FALSE, attr = "weight")
    return(out)
}

#Functions for randomTheta

collapse_scc <- function(g){
  
  comps <- components(g, mode = 'strong')
  dists <- distances(g, v = names(comps$membership), mode = 'out')
  colnames(dists) <- comps$membership
  rownames(dists) <- comps$membership
  dists[is.infinite(dists)] <- 0
  collapsed_g <- graph_from_adjacency_matrix(dists)
  #Got this next line from https://cutt.ly/gkM49x5 which removes duplicated nodes
  collapsed_g <- simplify(contract.vertices(collapsed_g, factor(V(collapsed_g)$name), 
                                            vertex.attr.comb = list(name='first')))
  return(list(collapsed_g, comps$membership))
  
}

leaf_cut <- function(phi){
  
  phi_g <- graph_from_adjacency_matrix(phi)
  comp_g <- collapse_scc(phi_g)
  
  g_to_cut <- comp_g[[1]]
  iter <- 1
  hierarchies <- list()
  
  while(length(V(g_to_cut)) != 0){
    
    degrees <- degree(g_to_cut, mode = 'out')
    remove_nodes <- names(degrees[degrees == 0])
    orig_nodes <- data.table::as.data.table(names(comp_g[[2]][comp_g[[2]] %in% remove_nodes]))
    orig_nodes[,position:=iter]
    hierarchies[[iter]] <- orig_nodes
    g_to_cut <- delete.vertices(g_to_cut, remove_nodes)
    iter <- iter + 1
  }
  
  hierarchies <- data.table::rbindlist(hierarchies)
  data.table::setnames(hierarchies, 'V1', 'node')
  
  return(hierarchies)
}

##Main

randomPhi <- function(n_nodes = 5, n_edges = 2*n_nodes, phi_type = 'barabasi', 
                      custom_phi = NULL, input_graph = NULL, 
                      sample_transitive_closed_input = TRUE, 
                      reset_probability = 0.1, 
                      return_input_distances = FALSE, ...){
  
  if(!is.null(custom_phi) & phi_type != 'custom'){
    warning('Set phi_type = "custom" to use a custom phi')
  }
  
  if(!is.null(input_graph) & !(phi_type %in% c('input_simple', 'input_near_seed', 'input_walk'))){
    warning('Set phi_type = "input_simple", "input_near_seed" or "input_walk" to sample an input graph')
  }
  
  stopifnot(phi_type %in% c('erdos_renyi', 'barabasi', 'preset', 'input_simple',
                            'input_near_seed', 'input_walk'))
  if(phi_type == 'erdos_renyi'){
    phi <- as_adj(sample_gnm(n_nodes, n_edges, directed = TRUE, ...),
                  sparse = FALSE)
    
  }else if(phi_type == 'custom'){
    stopifnot(!is.null(custom_phi))
    phi <- custom_phi
    stopifnot(all(diag(phi) == 1))
    
  }else if(phi_type == 'barabasi'){
    phi <- as_edgelist(sample_pa(n_nodes, directed = TRUE, ...))
    phi <- as_adj(graph_from_edgelist(phi[,c(2,1)]), sparse=FALSE)
    
  }else if(phi_type == 'input_simple'){
    phi <- simple_input_sample(n_nodes, input_graph, sample_transitive_closed_input,
                               return_input_distances)
    if(return_input_distances == TRUE){
      input_distances <- phi[[2]]
      phi <- phi[[1]]
    }
    
  }else if(phi_type == 'input_near_seed'){
    phi <- seed_neighbours(n_nodes, input_graph, sample_transitive_closed_input,
                           return_input_distances)
    if(return_input_distances == TRUE){
      input_distances <- phi[[2]]
      phi <- phi[[1]]
    }
    
  }else if(phi_type == 'input_walk'){
    phi <- walk_sample(n_nodes, input_graph, sample_transitive_closed_input,
                       reset_probability, return_input_distances)
    if(return_input_distances == TRUE){
      input_distances <- phi[[2]]
      phi <- phi[[1]]
    }
    
  }
  diag(phi) <- 1
  if(phi_type != 'custom'){
    rownames(phi) <- paste0('S', 1:n_nodes)
    colnames(phi) <- rownames(phi)
  }
  
  if(return_input_distances == TRUE){
    return(list(phi, input_distances))
  } else{
    return(phi)
  }
  
}

randomTheta <- function(phi, no_egenes = 3*dim(phi)[1], theta_type = 'uniform',
                        custom_theta = NULL){
  if(!is.null(custom_theta) & theta_type != 'custom'){
    warning('Set theta_type = "custom" to make use of a custom theta')
  }
  
  no_sgenes <- unique(dim(phi))
  
  stopifnot(theta_type %in% c('uniform', 'custom', 'deterministic', 'leaf_biased'))
  if(theta_type == 'uniform'){
    
    #Creating the adj matrix like this is faster than more obvious ways (say with table)
    attachments <- sample(no_sgenes, size = no_egenes, replace = TRUE)
    edj_list <- matrix(c(attachments, 1:no_egenes), ncol = 2)
    theta <- as_adj(graph_from_edgelist(edj_list), sparse = FALSE)[1:no_sgenes, 1:no_egenes]
    
  }else if(theta_type == 'custom'){
    stopifnot(!is.null(custom_theta))
    theta <- custom_theta
    
  }else if(theta_type == 'deterministic'){
    stopifnot('no_egenes must be a multiple no_sgenes' = !(no_egenes %% no_sgenes))
    
    egenes_per_sgene <- no_egenes/no_sgenes
    attachments <- rep(1:no_sgenes, egenes_per_sgene)
    edj_list <- matrix(c(attachments, 1:no_egenes), ncol = 2)
    theta <- as_adj(graph_from_edgelist(edj_list), sparse = FALSE)[1:no_sgenes, 1:no_egenes]
    
  }else if(theta_type == 'leaf_biased'){
    
    empty_net <- matrix(0, nrow = no_sgenes + no_egenes,
                        ncol = no_sgenes + no_egenes)
    rownames(empty_net) <- c(rownames(phi), 1:no_egenes)
    colnames(empty_net) <- rownames(empty_net)
    empty_net <- graph_from_adjacency_matrix(empty_net)
    
    hierarchies <- leaf_cut(phi)
    attachments <- sample(hierarchies$node, size = no_egenes, replace = TRUE, prob = 1/hierarchies$position)
    edj_list <- matrix(c(attachments, 1:no_egenes), ncol = 2)
    edges_to_add <- c(t(edj_list))
    
    theta <- add_edges(empty_net, edges_to_add)
    theta <- as_adj(theta, sparse = FALSE)[1:no_sgenes, no_sgenes + 1:no_egenes]
    colnames(theta) <- paste0('E', 1:no_egenes)
  }
  
  if(!(theta_type %in% c('custom', 'leaf_biased'))){
    rownames(theta) <- paste0('S', 1:no_sgenes)
    colnames(theta) <- paste0('E', 1:no_egenes)
  }
  
  stopifnot(all(colSums(theta) == 1))
  return(theta)
}

randomSignedPhi <- function(n_nodes = 5, n_edges = 2*n_nodes, phi_type = "barabasi", 
                          inhibit_prob = 0.3, input_graph = NULL, reset_probability = 0.1){
    
    phi <- randomPhi(n_nodes = n_nodes, phi_type = phi_type, input_graph = input_graph, reset_probability)
    diag(phi) <- 0
    ones <- which(phi == 1)
    signed_edges <- sample(c(-1,1), length(ones), replace = TRUE, prob = c(inhibit_prob, 1-inhibit_prob))
    phi[ones] <- signed_edges
    diag(phi) <- 1
    return(phi)
    
}

randomSignedTheta <- function(phi, no_egenes, theta_type = "uniform", inhibit_prob = 0.3){
    
    theta <- randomTheta(phi = phi, no_egenes = no_egenes, theta_type = theta_type)
    ones <- which(theta == 1)
    signed_edges <- sample(c(-1,1), length(ones), replace = TRUE, prob = c(inhibit_prob, 1-inhibit_prob))
    theta[ones] <- signed_edges
    return(theta)
}      

attachSignedTheta <- function(phi, negenes, nparents = 1, type = "stochastic", inhibit_prob = 0.3,
                             allow_saturate = FALSE){
    
    if(class(phi) != 'igraph'){
        phi_g <- graph_from_adjacency_matrix(phi, weighted = TRUE)
    }else{
        phi_g <- phi
    }
    
    nsgenes <- length(V(phi_g))
    
    if(nparents == 1){
        
        attachments <- sample(V(phi_g)$name, size = negenes, replace = TRUE)
        edge_list <- data.frame(attachments, paste0("E", 1:negenes))
        signed_edges <- sample(c(-1,1), nrow(edge_list), replace = TRUE, prob = c(inhibit_prob, 1-inhibit_prob))
        
    }else if(nparents > 1){
        
        if(type == "stochastic"){
            
            if(allow_saturate == FALSE){
                fail_n <- ppois(length(V(phi_g)), lambda = nparents)
                fail_n <- exp(negenes*log(fail_n))
                fail_n <- round(1/(1-fail_n))
            
                if(fail_n <= 100000){
                    warning(paste0("This will fail in roughly 1 in every ", fail_n, " runs.",
                                  " Consider either increasing the number of S-genes, decreasing nparents,",
                                   " decreasing number of E-genes,",
                                   " or switching to deterministic attachment."))
                }
            }else{
                warning("allow_saturate = TRUE means that the average number of parents could be pretty different to nparents")
            }
            
            nparents_vec <- 1 + rpois(negenes, nparents - 1)
            
        }else if(type == "deterministic"){
            stopifnot("nparents must be an integer for deterministic attachment" = nparents%%1 == 0)
            nparents_vec <- rep(nparents, times = negenes)
        }
        
        edge_list <- lapply(nparents_vec, function(e) sample(V(phi_g)$name, size = e, replace = allow_saturate))
        edge_list <- reshape2::melt(edge_list)
        edge_list$L1 <- paste0("E", edge_list$L1)
        signed_edges <- sample(c(-1,1), nrow(edge_list), replace = TRUE, prob = c(inhibit_prob, 1-inhibit_prob))
        
    }else{
        stop("Invalid number of parents: must be > 1")
    }
                            
    phi_g <- add_vertices(phi_g, negenes, attr = list("name" = paste0("E", 1:negenes)))
    phi_g <- add_edges(phi_g, c(t(edge_list)), attr = list("weight" = signed_edges))
                            
    theta <- as_adj(phi_g, sparse = FALSE, attr = "weight")[1:nsgenes, nsgenes + 1:negenes]
    return(list("full_net" = phi_g, "theta" = theta))
    
}    
    
distanceDependentErrors <- function(phi_orig, theta, mode = 'phi',
                                    input_distances = NULL, max_prop_dist = nrow(phi_orig)){
  
  stopifnot("mode must be either 'phi' or 'input'" = mode %in% c('phi', 'input'))
  
  #Generates the Fmat but with entries = to prior probability of observing
  #an effect
  
  if(mode == 'input'){
    stopifnot('Please provide input graph distances' = !is.null(input_distances))
    input_distances[is.infinite(input_distances)] <- 0 #turn it into weighted adjacency matrix
    diag(input_distances) <- 1
    phi_orig <- input_distances
  }
  
  full_net <- cbind(phi_orig, theta)
  the_rest <- matrix(0, nrow = ncol(theta), ncol = ncol(full_net))
  rownames(the_rest) <- colnames(theta)
  full_net <- rbind(full_net, the_rest)
  
  full_net_g <- graph_from_adjacency_matrix(full_net, weighted = TRUE) #so for input we take weights for distance
  
  path_lengths <- distances(full_net_g, v = rownames(phi_orig),
                            to = colnames(theta),
                            mode = 'out')
  
  path_lengths <- path_lengths - 1 #we consider an E-gene directly attached as having distance 0
  prior_effects <- 1 - path_lengths/max_prop_dist
  prior_effects[prior_effects < 0 | is.na(prior_effects)] <- 0 #any Infs are now negative so also removes them
  noisy_out <- (prior_effects > runif(length(prior_effects)))*1
  return(noisy_out)
}

## Simulate binary is the fundamental simulation function as with alpha = beta = 0,
## the function returns the F = ???? matrix.

simulate_binary <- function(alpha = 0.15, beta = 0.05, no_sgenes = 5, no_egenes = 3*no_sgenes, 
                            no_phi_edges = 2*no_sgenes, no_replicates = 1, phi_type = 'barabasi', 
                            theta_type = 'uniform', custom_phi = NULL, custom_theta = NULL, 
                            transitive_closure = TRUE, input_graph = NULL, 
                            distance_depn_error = FALSE, dist_error_mode = 'phi', 
                            max_prop_dist = no_sgenes, sample_transitive_closed_input = TRUE, 
                            reset_probability = 0.1){
  
  if(phi_type == 'custom'){
    stopifnot(!is.null(custom_phi))
    no_sgenes <- dim(custom_phi)[1]
  }
  
  if(theta_type == 'custom'){
    stopifnot(!is.null(custom_theta))
    no_egenes <- ncol(custom_theta)[1]
  }
  
  if(distance_depn_error == TRUE && dist_error_mode == 'input'){
    stopifnot("Set phi_type = 'input_*' if you want input distant dependent errors" = phi_type %in% c('input_simple',
                                                                                                      'input_walk',
                                                                                                      'input_simple'))
    return_input_distances <- TRUE
    
    out <- randomPhi(no_sgenes, no_phi_edges, phi_type, custom_phi, input_graph,
                     sample_transitive_closed_input, reset_probability,
                     return_input_distances)
    phi <- out[[1]]
    input_distances <- out[[2]]
  }else{
    
    phi <- randomPhi(no_sgenes, no_phi_edges, phi_type, custom_phi, input_graph,
                     sample_transitive_closed_input, reset_probability,
                     return_input_distances = FALSE)
    input_distances <- NULL
  }
  
  
  phi_orig <- phi
  phi_transitive_closure <- NULL
  
  theta <- randomTheta(phi, no_egenes, theta_type, custom_theta)
  
  if(transitive_closure == TRUE){
    phi <- suppressWarnings(transitive.closure(phi))
    phi <- as(phi, 'matrix')
    phi_transitive_closure <- phi
  }
  
  if(distance_depn_error == TRUE){
    stopifnot("transitive_closure must be TRUE for distance dependent errors" = transitive_closure)
    Fmat <- distanceDependentErrors(phi_orig, theta, dist_error_mode, input_distances,
                                    max_prop_dist)
  } else{
    Fmat <- phi %*% theta
  }
  
  Fmat <- Fmat[rep(1:no_sgenes, each=no_replicates),]
  rownames(Fmat) <- paste0(rownames(Fmat), rep(paste0("_", 1:no_replicates), no_sgenes))
  
  ones <- which(Fmat == 1)
  zeros <- which(Fmat == 0)
  
  typeI <- (runif(length(ones)) > alpha)*1
  typeII <- (runif(length(zeros)) < beta)*1
  
  Fmat[ones] <- typeI
  Fmat[zeros] <- typeII
  out <- list('simulated_data' = Fmat, 'phi' = phi_orig,
              'theta' = theta, 'phi_transitive_closure' = phi_transitive_closure,
              parameters = list('alpha' = alpha, 'beta' = beta,
                                'no_sgenes' = no_sgenes, 'no_egenes' = no_egenes,
                                'no_phi_edges' = no_phi_edges,
                                'phi_type' = phi_type,
                                'theta_type' = theta_type,
                                'distance_depn_error' = distance_depn_error,
                                'gave_input' = !is.null(input_graph),
                                'dist_error_mode' = dist_error_mode, 
                                'max_prop_dist' = max_prop_dist, 
                                'sample_transitive_closed_input' = sample_transitive_closed_input, 
                                'reset_probability' = reset_probability))
  return(out)
}

simulate_logratio_2008 <- function(alpha = 0.25, no_sgenes = 5, no_egenes = 3*no_sgenes, 
                                   no_phi_edges = 2*no_sgenes, no_replicates = 1, 
                                   phi_type = 'barabasi', theta_type = 'uniform', 
                                   custom_phi = NULL, custom_theta = NULL, 
                                   transitive_closure = TRUE, distance_depn_error = FALSE,
                                   dist_error_mode = 'phi', max_prop_dist = no_sgenes,
                                   input_graph = NULL, sample_transitive_closed_input = TRUE, 
                                   reset_probability = 0.1){
  
  out <- simulate_binary(alpha = 0, beta = 0,
                         no_sgenes = no_sgenes,
                         no_egenes = no_egenes,
                         no_phi_edges = no_phi_edges,
                         no_replicates = no_replicates,
                         phi_type = phi_type,
                         theta_type = theta_type,
                         custom_phi = custom_phi,
                         custom_theta = custom_theta,
                         transitive_closure = transitive_closure, 
                         distance_depn_error = distance_depn_error,
                         input_graph = input_graph,
                         dist_error_mode = dist_error_mode, 
                         max_prop_dist = max_prop_dist,
                         sample_transitive_closed_input = sample_transitive_closed_input,
                         reset_probability = reset_probability)
  
  sim_data <- out$simulated_data - 0.5
  sim_data <- sim_data + rnorm(length(sim_data), 0, alpha)
  
  out$simulated_data <- sim_data
  out$parameters$alpha <- alpha
  out$parameters$beta <- NULL
  return(out)
}

.return_pvaldensity_row <- function(row, lambda_sd, alpha_sd, beta_sd){
  
  #This is a strange function. It uniformly samples parameter values for a beta-
  #uniform mixture model. It then saves these parameters, then samples from the
  #beta-uniform mixture using the parameters with normal noise added. The resulting
  #'p-values' are then converted to densities using the mixture model parameterised
  #without the normal noise. In this way, it attempts to simulate a discrepancy
  #between the fitted distribution and the true distribution.
  
  lambda <- numeric(3)
  #Sampling pik2 and pik3 from (0,0.5) is what Frolich did. Not sure why...
  #I suppose to ensure their sum is never > 1, but there other ways of doing this
  #which don't constrain both values from being <=0.5 which seems arbitrary to me
  lambda[2:3] <- runif(2, min = 0, max = 0.5) 
  lambda[1] <- 1 - lambda[2] - lambda[3]
  
  a <- numeric(2)
  a[1] <- runif(1)
  a[2] <- runif(1, min = 5, max = 50)
  
  lambda_blur <- numeric(3)
  lambda_blur[2:3] <- abs(lambda[2:3] + rnorm(2, 0, lambda_sd))
  lambda_blur[1] <- 1 - lambda_blur[2] - lambda_blur[3]
  
  a_blur <- numeric(2)
  a_blur_error <- rnorm(1, 0, alpha_sd)
  #alpha can't be < 0 or > 1 otherwise will give p-value = 0 and Inf density
  a_blur[1] <- ifelse(a[1] + a_blur_error < 0, a[1] - a_blur_error,
                      ifelse(a[1] + a_blur_error > 1, a[1] - a_blur_error,
                             a[1] + a_blur_error))
  a_blur[2] <- a[2] + rnorm(1, 0, a[2]*beta_sd)
  
  ones <- which(row == 1)
  zeros <- which(row == 0)
  
  alt_samples <- bum.ralt(sum(row), a_blur, lambda_blur)
  null_samples <- runif(sum(!row))
  
  row[ones] <- alt_samples
  row[zeros] <- null_samples
  
  row <- dbum(row, a, lambda)
  
  return(row)
}

simulate_pvaldensity_2008 <- function(no_sgenes = 5, no_egenes = 3*no_sgenes, 
                                      no_phi_edges = 2*no_sgenes, no_replicates = 1,
                                      lambda_sd = 0.05, alpha_sd = 0.05, beta_sd = 0.1,
                                      phi_type = 'barabasi', theta_type = 'uniform', 
                                      custom_phi = NULL, custom_theta = NULL, 
                                      transitive_closure = T, distance_depn_error = FALSE,
                                      dist_error_mode = 'phi', max_prop_dist = no_sgenes,
                                      input_graph = NULL, sample_transitive_closed_input = TRUE, 
                                      reset_probability = 0.1){
  
  out <- simulate_binary(alpha = 0, beta = 0,
                         no_sgenes = no_sgenes,
                         no_egenes = no_egenes,
                         no_phi_edges = no_phi_edges,
                         no_replicates = no_replicates,
                         phi_type = phi_type,
                         theta_type = theta_type,
                         custom_phi = custom_phi,
                         custom_theta = custom_theta,
                         transitive_closure = transitive_closure,
                         distance_depn_error = distance_depn_error,
                         input_graph = input_graph,
                         dist_error_mode = dist_error_mode, 
                         max_prop_dist = max_prop_dist,
                         sample_transitive_closed_input = sample_transitive_closed_input,
                         reset_probability = reset_probability)
  
  sim_data <- out$simulated_data
  out$simulated_data <- t(apply(sim_data, 1, .return_pvaldensity_row, #key function here
                                lambda_sd, alpha_sd, beta_sd)) 
  out$simulated_data <- log(out$simulated_data)
  out$parameters$alpha <- NULL
  out$parameters$beta <- NULL
  
  out$parameters$alpha_sd <- alpha_sd
  out$parameters$beta_sd <- beta_sd
  out$parameters$lambda_sd <- lambda_sd
  
  return(out)
  
}
    
simulate_count_matrix <- function(meanDispPairs, ngenes = 15000, nsgenes = 15, ndeggenes = 1000, nreps = 3,
                                 foldchange_range = c(1.5,6), prop_informative = 0.95, 
                                 inhibit_prob = 0.3, nparents = 1, nedges = 2*nsgenes, 
                                 phi_type = "barabasi", theta_type = "uniform",
                                 eparent_attachment_type = "stochastic",
                                 input_graph = NULL, reset_probability = 0.1,
                                 mean_col = "norm_mean", disp_col = "disp",
                                 transitive_closure = TRUE, allow_saturated = FALSE){
    
    negenes <- round(prop_informative*ndeggenes)
    
    #Generate random graph with signed edges:
    phi <- randomSignedPhi(n_nodes = nsgenes, phi_type = phi_type, input_graph = input_graph, inhibit_prob = inhibit_prob)
    
    phi_orig <- phi
    transitive_closure_phi <- NULL
    
    if(nparents == 1){
        
        if(transitive_closure == TRUE){
            phi <- signed_transitive_closure(phi)
            transitive_closure_phi <- phi
        }
        
        theta <- randomSignedTheta(phi, no_egenes = negenes, theta_type = theta_type, inhibit_prob = inhibit_prob)
        t_Fmat <- t(-1*(phi %*% theta)) #The informative DEGs
        
    }else{
        
        full_net_list <- attachSignedTheta(phi, negenes = negenes, type = eparent_attachment_type, nparents = nparents,
                                          allow_saturate = allow_saturate)
        theta <- full_net_list$theta
        
        extended_Fmat <- signed_transitive_closure(as_adj(full_net_list$full_net, attr = "weight", sparse = FALSE))
        transitive_closure_phi <- extended_Fmat[1:nsgenes, 1:nsgenes]
        
        t_Fmat <- t(-1*extended_Fmat[1:nsgenes, nsgenes + 1:negenes])
    }
    
    if(prop_informative < 1){ #Adds in some uniformative DEGs
    notinform <- t(randomSignedTheta(phi, no_egenes = ndeggenes - negenes, theta_type = "uniform"))
    t_Fmat <- rbind(t_Fmat, notinform)
    }
                                     
    degs <- which(t_Fmat != 0)
    degFolds <- runif(n = length(degs), min = foldchange_range[1], max = foldchange_range[2])
    t_Fmat[degs] <- t_Fmat[degs]*log2(degFolds)
    
    #Add the non-differential genes
    notdegs <- matrix(rep(0, nsgenes*(ngenes - nrow(t_Fmat))), ncol = nsgenes)
    beta <- rbind(t_Fmat, notdegs)
    rownames(beta) <- paste0("gene", 1:nrow(beta))
    
    #Sample from mean-dispersion pairs as basis for new dataset:
    idx <- sample(nrow(meanDispPairs), ngenes, replace=TRUE)
    mu0 <- meanDispPairs[idx,mean_col]
    disp <- meanDispPairs[idx,disp_col]
    
    #Add baseline "null" expression column
    beta <- cbind(log2(mu0), beta)
    
    #Now we need to create a design matrix to apply the fold changes
    cond_names <- c("CTRL", paste0("S", 1:nsgenes))
    condition <- factor(rep(cond_names, each = nreps), levels = cond_names)
    mod_mat <- model.matrix(~ condition)
    mu <- 2^(beta %*% t(mod_mat))
    
    #Finally, sample from Negative Binomial to get data:
    sim_data <- matrix(rnbinom(prod(dim(mu)), mu = mu, size = 1/disp), nrow = ngenes)
    return(list("sim_data" = sim_data, "col_data" = condition, "mu0" = mu0, "disp" = disp,
               "phi" = phi_orig, "transitive_closure_phi" = transitive_closure_phi,
               "theta" = theta))
}
    
###Maybe not their permanent home, but here are some method functions:
    
his <- function(matData, intMin = 1.5, intMax = 5, 
                intSteps = 200, strTail = "both",
                returnMatHISPerAssay = FALSE){
  
  if(any(is.na(matData))){
    warning("Imputing NAs for zeros!")
    matData[is.na(matData)] <- 0
  }
  min_bool <- abs(matData) >= intMin
  keepRows <- rowSums(min_bool) > 0
  keepCols <- colSums(min_bool) > 0
  matData <- matData[keepRows, keepCols]
  
  if(is.null(dim(matData))){
    out <- matrix(0, nrow = length(keepRows), ncol = length(keepCols))
    out <- Matrix::Matrix(out, sparse = TRUE)
    return(out)
  }
  
  if(strTail == 'single'){
    matTs <- seq(from = intMin, to = intMax,
                 length.out = intSteps/2)
  }else{
    matTs <- c(seq(to = -intMin, from = -intMax,
                 length.out = intSteps/2),
               seq(from = intMin, to = intMax,
                   length.out = intSteps/2))
  }
  
  data_dims <- dim(matData)
  iMaxVal <- length(matTs)
  iNodes <- data_dims[1]
  iAssays <- data_dims[2]
  
  if(returnMatHISPerAssay == TRUE){
    matSumMinHitCounts <- array(0, c(iNodes, iNodes, iAssays))
  }else{
    matSumMinHitCounts <- matrix(0, nrow = iNodes, ncol = iNodes)
  }
  
  for(t in 1:iMaxVal){
    
    if(matTs[t] <= 0){
      matHits <- matData <= matTs[t]
    }else{
      matHits <- matData >= matTs[t]
    }
    
    matPatterns <- unique(matHits)
    
    #First divergence from MATLAB code here - making a list with each pattern
    #as a column so I can associate an index with each row in matData. This
    #bizarre code loops through every row of matHits, then for each row I
    #iterate through each row of matPatterns and check which row of matPatterns
    #matches the given row in matHits and return the index.
    
    #Unfortunately this is a much slower line than the corresponding one in MATLAB
    
    matGenePatternMapping <- unlist(lapply(1:nrow(matHits), 
           function(i) which(unlist(lapply(1:nrow(matPatterns), 
           function(j, i) all(matPatterns[j,] == matHits[i,]), i = i)))))
    
    matPatternHitCount <- rowSums(matPatterns)
    
    boolChildren <- matrix(FALSE, nrow = nrow(matPatterns), ncol = nrow(matPatterns))
    boolParents <- boolChildren
    
    for(i in 1:nrow(matPatterns)){
      
      x2 <- matPatterns[i,]
      
      #Reason for t() is to ensure output is a matrix so colSums never fails
      matChildCount <- colSums(t(matPatterns[,x2]))
      matNotChildCount <- matPatternHitCount - matChildCount
      
      matCandidateChildren <- matChildCount>0 & matNotChildCount==0
      matCandidateParents <- matChildCount==sum(x2) & matNotChildCount>0
      
      boolChildren[i, matCandidateChildren] <- TRUE
      boolParents[i, matCandidateParents] <- TRUE
      
    }
    
    for(i in 1:nrow(matPatterns)){
      c <-  which(boolChildren[i,])
      for(c2 in c){
        p <- which(boolParents[c2,])
        p <- p[!(p == i)] #This is what p(p==i) = []; does in MATLAB
        
        if(!any(p %in% c)){
          
          #This next line is needed in R so that if c2 is empty
          #`xj <- matGenePatternMapping == c2` returns a list of FALSEs
          #as occurs in MATLAB
          c2 <- ifelse(length(c2) == 0, 0, c2)
          xi <- matGenePatternMapping == i
          xj <- matGenePatternMapping == c2
          if(returnMatHISPerAssay == TRUE){
            matSumMinHitCounts[xi,xj,matPatterns[c2,]] <- matSumMinHitCounts[xi,xj,matPatterns[c2,]] + 1
          }else{
            matSumMinHitCounts[xi,xj] <- matSumMinHitCounts[xi,xj] + matPatternHitCount[c2]
          }
        }
      }
    }
  }
  
  if(returnMatHISPerAssay == TRUE){
    
    #There is a bug in the MATLAB code here which I haven't fixed yet
    #because it's difficult to know what the authors wanted the code to do
    matHISPerAssay <- array(0, c(iNodes,iNodes,iAssays))
    matHISPerAssay[keepRows, keepRows, keepCols] = matSumMinHitCounts
    
    matSumMinHitCounts <- rowSums(matSumMinHitCounts, 2)
  }
  
  diag(matSumMinHitCounts) <- 0
  idxs <- which(matSumMinHitCounts != 0, arr.ind = TRUE)
  vals <- matSumMinHitCounts[matSumMinHitCounts != 0]
  
  d <- which(keepRows)
  matHIS <- Matrix::sparseMatrix(i = d[idxs[,1]], j = d[idxs[,2]],
                                 x = vals/iMaxVal,
                                 dims = c(length(keepRows), 
                                          length(keepRows)))
  
  if(returnMatHISPerAssay == TRUE){
    return(list("matHIS" = matHIS, "matHISPerAssay" = matHISPerAssay)) 
  }else{
    return(matHIS)
  }
}