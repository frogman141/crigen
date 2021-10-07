#' @export
run_simulation <- function(model,
                           noise_mean=1,
                           noise_sd=0.005,
                           nsims_per_grna=100,
                           burn_time=100,
                           sampling_time=500,
                           keep_burn=FALSE,
                           tau=30/3600,
                           census_interval=1,
                           verbose=FALSE) {
    
    # setting up local or cluster parallelisms 
    model <- initialize_simulation_params(model, noise_mean, noise_sd,
                                          nsims_per_grna, burn_time, sampling_time,
                                          keep_burn, tau, census_interval)
    
    message ("Precompile Reactions...")
    model <- precompile_reactions(model)
    sims <- get_sims(model)
    
    message ("Simulating Control Population...")
    simulations <- pbapply::pblapply(X = sims$ctrl_sims, FUN = simulate_cell, model=model, compile_reactions=FALSE, ctrl_sim=TRUE, verbose=verbose)
    model$ctrl_sim <- extract_simulation_results(simulations, 'sampling')
    
    message ("Simulating Perturbations...")
    simulation <- pbapply::pblapply(X = sims$grna_sims, FUN = simulate_cell, model=model, compile_reactions=FALSE, ctrl_sim=FALSE, verbose=verbose)
    model$perturb_sim <- extract_simulation_results(simulation, 'sampling')
    
    model <- merge_simulations(model, model$ctrl_sim, model$perturb_sim)
    return (model)
}

#' @export
initialize_simulation_params <- function(model, noise_mean, noise_sd, nsims_per_grna, burn_time,
                                         sampling_time, keep_burn, tau, census_interval) {
    
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
    sim_system$sampling_time <- sampling_time
    sim_system$nsims_per_grna <- nsims_per_grna
    sim_system$census_interval <- census_interval
    
    model$simulation <- list()
    model$simulation_system <- sim_system
    return (model)
}

#' @export
simulate_cell <- function(sim, model, compile_reactions, ctrl_sim, verbose) {
    # get reusable parameters
    if (compile_reactions) {
        model <- precompile_reactions(model)
    }
    
    perturb_params <- get_parameters(sim, model)
       
    model <- model %>% burn_phase(sim, perturb_params, verbose)
    model <- model %>% sampling_phase(sim, perturb_params, ctrl_sim, verbose)
    
    return (model$simulation)
}

#' @export
merge_simulations <- function(model, ctrl_sim, perturb_sim) {
    sim <- list()
    
    sim$meta <- bind_rows(ctrl_sim$meta, perturb_sim$meta) %>% as.data.frame %>% mutate(step_ix=row_number())
    sim$counts <- rbind(ctrl_sim$counts, perturb_sim$counts)
    
    model$simulation <- sim
    model <- model[!names(model) %in% c('ctrl_sim', 'perturb_sim')]
    return (model)
} 

#' @export
get_sims <- function(model) {
    sim_names <- model$simulation_system$sim_names
    ctrl_sims <- sim_names[grepl("CTRL", sim_names)]
    grna_sims <- sim_names[!grepl("CTRL", sim_names)]
    
    return (lst(ctrl_sims, grna_sims))
}

get_parameters <- function(sim, model) {
    cell_params <- inject_sim_noise(sim, model) %>% extract_row_to_vector(sim)
    
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
       
#' @export
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