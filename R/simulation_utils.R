#' @export
precompile_reactions <- function(model) {
    # fetch paraneters and settings
    sim_system <- model$simulation_system
    reactions <- sim_system$reactions 
    state_ids <- names(sim_system$initial_state)
    buffer_ids <- unique(unlist(map(reactions, "buffer_ids")))
    parameters <- extract_parameters(model$feature_info, model$feature_network) %>% deframe()
    
    # compile prop funs
    comp_funs <- GillespieSSA2::compile_reactions(reactions = reactions,
                                                  buffer_ids = buffer_ids,
                                                  state_ids = state_ids,
                                                  params = parameters,
                                                  hardcode_params = FALSE,
                                                  fun_by = 1000L)
    
    model$simulation_system$compiled_reactions <- comp_funs
    return (model)
}

#' @export
burn_phase <- function(model, sim, params) {
    sim_system <- model$simulation_system
    
    # save results of burn phase and set sim_time to NULL
    model$simulation <- list()
    burn_time <- sim_system$burn_time
    initial_burn_state <- sim_system$initial_state
    
    burn_out <- run_gillispie(sim, 'burn', params, burn_time, initial_burn_state, sim_system)
    sim_system$initial_perturb_state <- burn_out$counts[nrow(burn_out$counts),]
    
    if (sim_system$keep_burn) {
        burn_out$meta$crispr_type <- model$crispr_param$crispr_type
        model$simulation$burn <- burn_out
    }
    
    model$simulation_system <- sim_system
    
    return (model)
}
     
#' @export
perturb_phase <- function(model, sim, params) {
    sim_system <- model$simulation_system
    
    # save results of burn phase and set sim_time to NULL
    perturb_time <- sim_system$perturb_time
    initial_perturb_state <- sim_system$initial_perturb_state
    
    perturb_out <- run_gillispie(sim, 'perturb', params, perturb_time, initial_perturb_state, sim_system)
    sim_system$initial_sampling_state <- perturb_out$counts[nrow(perturb_out$counts),]
    
    if (sim_system$keep_perturb) {
        perturb_out$meta$crispr_type <- model$crispr_param$crispr_type
        model$simulation$perturb <- perturb_out
    }
    
    model$simulation_system <- sim_system
    
    return (model)
}
 
#' @export
sampling_phase <- function(model, sim, params) {
    sim_system <- model$simulation_system
    
    # save results of burn phase and set sim_time to NULL
    sampling_time <- sim_system$sampling_time
    initial_sampling_state <- sim_system$initial_sampling_state
    
    sampling_out <- run_gillispie(sim, 'sampling', params, sampling_time, initial_sampling_state, sim_system)
    sampling_out$meta$crispr_type <- model$crispr_param$crispr_type
    model$simulation$sampling <- sampling_out
    
    return (model)
}


run_gillispie <- function(sim, phase, params, sim_time, initial_state, sim_system) {
    
    ssa_algo <- ssa_etl(tau=sim_system$tau)
    reactions <- sim_system$compiled_reactions
    census_interval <- sim_system$census_interval
    
    # simulate cell
    out <- GillespieSSA2::ssa(sim_name = sim,
                              verbose = FALSE,
                              params = params,
                              method = ssa_algo,
                              reactions = reactions,
                              final_time = sim_time,
                              stop_on_neg_state = FALSE,
                              initial_state = initial_state,
                              census_interval = census_interval)
    
    out <- process_ssa_output(sim, phase, out, sim_time)
    colnames(out$counts) <- names(initial_state)
    
    return (out)
}
                        
process_ssa_output <- function(sim, phase, out, sim_time) {
    meta <- tibble(time = c(head(out$time, -1), sim_time)) %>%
                mutate(phase=phase, sim_name=sim, time=ceiling(time),
                       target_gene=str_split(sim_name, '-')[[1]][1]) %>%
                rowwise %>%
                mutate(grna=ifelse(target_gene == 'CTRL', str_split(sim_name, '-')[[1]][1],
                                   paste(str_split(sim_name, '-')[[1]][1:2], collapse='-')))
    
    counts <- out$state %>% Matrix::Matrix(sparse = TRUE)
    
    return(lst(meta, counts))
}
     
