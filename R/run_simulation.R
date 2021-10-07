#' @export
run_simulation <- function(model,
                           noise_mean=1,
                           noise_sd=0.005,
                           nsims_per_grna=100,
                           burn_time=200,
                           perturb_time=400,
                           sampling_time=400,
                           keep_burn=FALSE,
                           keep_perturb=FALSE,
                           tau=30/3600,
                           census_interval=2) {
    
    # setting up local or cluster parallelisms 
    model <- initialize_simulation_params(model, noise_mean, noise_sd,
                                          nsims_per_grna, burn_time,
                                          perturb_time, sampling_time,
                                          keep_burn, keep_perturb, tau,
                                          census_interval)
    
    message ("Precompile Reactions...")
    model <- precompile_reactions(model)
    sim_names <- model$simulation_system$sim_names
    
    message ("Running Simulation...")
    simulation <- pbapply::pblapply(sim_names, simulate_cell, model=model, compile_reactions=FALSE)
    
    model <- merge_simulations(model, simulation)
    return (model)
}

#' @export
initialize_simulation_params <- function(model, noise_mean, noise_sd, nsims_per_grna,
                                         burn_time, perturb_time, sampling_time, keep_burn,
                                         keep_perturb, tau, census_interval) {
    
    sim_system <- model$simulation_system
    sim_names <- lapply(rownames(sim_system$multiplier),
                        function(x) paste0(x, '-sim.', 1:nsims_per_grna)) %>%
                        unlist()
    
    sim_system$tau <- tau
    sim_system$noise_sd <- noise_sd
    sim_system$sim_names <- sim_names
    sim_system$burn_time <- burn_time
    sim_system$keep_burn <- keep_burn
    sim_system$noise_mean <- noise_mean
    sim_system$keep_perturb <- keep_perturb
    sim_system$perturb_time <- perturb_time
    sim_system$sampling_time <- sampling_time
    sim_system$nsims_per_grna <- nsims_per_grna
    sim_system$census_interval <- census_interval
    
    model$simulation_system <- sim_system
    return (model)
}

#' @export
simulate_cell <- function(sim, model, compile_reactions) {
    # get reusable parameters
    if (compile_reactions) {
        model <- precompile_reactions(model)
    }
    
    cell_params <- inject_sim_noise(sim, model) %>% extract_row_to_vector(sim)
    perturb_params <- calculate_perturb_parameters(sim, model, cell_params)
       
    model <- model %>% burn_phase(sim, cell_params)
    model <- model %>% perturb_phase(sim, perturb_params) 
    model <- model %>% sampling_phase(sim, perturb_params)
    
    return (model$simulation)
}

#' @export
merge_simulations <- function(model, simulation) {
    sim <- list()
    sim_system <- model$simulation_system
    
    if (sim_system$keep_burn) {
        sim$burn <- extract_simulation_results(simulation, 'burn')
    }
    
    if (sim_system$keep_perturb) {
        sim$perturb <- extract_simulation_results(simulation, 'perturb')
    }
    
    sim$sampling <- extract_simulation_results(simulation, 'sampling')
    sim$sampling$meta <- sim$sampling$meta %>% as.data.frame %>% mutate(step_ix=row_number())
    
    kinetics_parameters <- extract_simulation_results(simulation, 'kinetics', FALSE)
    sim_system$kinetics_parameters <- kinetics_parameters %>% bind_rows()
    
    model$simulation <- sim
    model$simulation_system <- sim_system
    return (model)
} 
                        
calculate_perturb_parameters <- function(sim, model, cell_params) {
    
    if (grepl('grna', sim)) {
        
        if (grepl('PRT|NT', sim)) {
            grna <- paste(str_split(sim, '-')[[1]][1:3], collapse='-')
        } else {
            grna <- paste(str_split(sim, '-')[[1]][1:2], collapse='-')
        }
        
        cell_multi <- model$simulation_system$multiplier %>% extract_row_to_vector(grna)
        kd_wprs <- names(cell_multi)
        perturb_params <- cell_params
        perturb_params[kd_wprs] <- perturb_params[kd_wprs] * cell_multi
    } else {
        perturb_params <- cell_params 
    }
    
    return (perturb_params)
}

extract_row_to_vector <- function(df, row_name) {
    row_vector <- df[row_name, ] %>%
                    t %>%
                    as.data.frame %>%
                    rename(param_value=colnames(.)[1]) %>% 
                    mutate(label=rownames(.)) %>% 
                    select(label, param_value) %>%
                    deframe()
    
    return (row_vector)
}
                        
extract_simulation_results <- function (sim, extractor, map_results=TRUE) {
    results <- list()
    
    for (sim_i in 1:length(sim)) {
        results[[sim_i]] <- sim[[sim_i]][[extractor]]
    }
    
    if (isTRUE(map_results)) {
        results <- lst(meta=map_df(results, "meta"), counts=do.call(rbind, map(results, "counts")))
        results$meta <- results$meta %>% as.data.frame %>% mutate(step_ix=row_number())
    }
    
    return (results)
} 