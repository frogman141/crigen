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
                           census_interval=2,
                           compute_rna_velocity=FALSE,
                           ncpus=1,
                           mem='4G',
                           cluster='local',
                           time='24:00:00',
                           template="../templates/slurm.tmpl",
                           rscript=system("which Rscript", intern=TRUE)) {
    
    # setting up local or cluster parallelisms 
    model <- initialize_simulation_params(model, noise_mean, noise_sd,
                                           nsims_per_grna, burn_time,
                                           perturb_time, sampling_time,
                                           keep_burn, keep_perturb, tau,
                                           census_interval, compute_rna_velocity,
                                           cluster, ncpus, mem, time, template, rscript)
    
    sim_names <- model$simulation_system$sim_names
    
    if (cluster == 'local') {
        model <- precompile_reactions(model)
        cl <- model$simulation_system$resources$cl
        simulation <- pblapply(sim_names, simulate_cell, model=model, cl=cl)
    } else {
        njobs <- model$simulation_system$resources$njobs
        simulation <- run_jobs(sim_names, simulate_cell, njobs=njobs,
                               model=model, inject_noise=TRUE, compile_reactions=TRUE)
    }
    
    model <- merge_simulations(model, simulation)
    return (model)
}

initialize_simulation_params <- function(model, noise_mean, noise_sd, nsims_per_grna,
                                 burn_time, perturb_time, sampling_time, keep_burn,
                                 keep_perturb, tau, census_interval, compute_rna_velocity,
                                 cluster, ncpus, mem, time, template, rscript) {
    
     # todo list: need to add torque, lsf, sge
    assert_that(cluster %in% c('local', 'slurm'),
                msg=paste0("Parallelism Option Provided is not valid. Options are: local, Slurm"))
    
    sim_system <- model$simulation_system
    sim_names <- lapply(rownames(sim_system$multiplier),
                        function(x) paste0(x, '-sim.', 1:nsims_per_grna)) %>%
                        unlist()
    
    default_resources <- list(ncpus=ncpus, njobs=length(sim_names), memory=mem, walltime=time, rscript=rscript)
    
    setup_cluster(cluster, template, default_resources)
    
    sim_system$tau <- tau
    sim_system$cluster <- cluster
    sim_system$noise_sd <- noise_sd
    sim_system$sim_names <- sim_names
    sim_system$burn_time <- burn_time
    sim_system$keep_burn <- keep_burn
    sim_system$noise_mean <- noise_mean
    sim_system$keep_perturb <- keep_perturb
    sim_system$perturb_time <- perturb_time
    sim_system$resources <- default_resources
    sim_system$sampling_time <- sampling_time
    sim_system$nsims_per_grna <- nsims_per_grna
    sim_system$census_interval <- census_interval
    sim_system$compute_rna_velocity <- compute_rna_velocity
    
    model$simulation_system <- sim_system
    return (model)
}
                     
simulate_cell <- function(sim, model, compile_reactions, inject_noise) {
    # get reusable parameters
    sim_system <- model$simulation_system
    cell_params <- inject_sim_noise(sim, model) %>% extract_row_to_vector(sim)
    
    if (compile_reactions) {
        model <- precompile_reactions(model)
    }
    
    if (grepl('grna', sim)) {
        grna <- paste(str_split(sim, '-')[[1]][1:2], collapse='-')
        cell_multi <- sim_system$multiplier %>% extract_row_to_vector(grna)
        perturb_params <- calculate_perturb_parameters(cell_multi, cell_params)
    } else {
        perturb_params <- cell_params 
    }
    
    model <- model %>% 
                burn_phase(sim, cell_params) %>%
                perturb_phase(sim, perturb_params) %>%  
                sampling_phase(sim, perturb_params)
    
    model$simulation$kinetics <- cell_params
    return (model$simulation)
}

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
                        
calculate_perturb_parameters <- function(cell_multi, cell_params){
    kd_wprs <- names(cell_multi)
    perturb_params <- cell_params
    perturb_params[kd_wprs] <- perturb_params[kd_wprs] * cell_multi
    
    return (perturb_params)
}
          
extract_simulation_results <- function (sim, extractor, map_results=TRUE) {
    results <- list()
    
    for (sim_name in names(sim)) {
        results[[sim_name]] <- sim[[sim_name]][[extractor]]
    }
    
    if (isTRUE(map_results)) {
        results <- lst(meta=map_df(results, "meta"), counts=do.call(rbind, map(results, "counts")))
        results$meta <- results$meta %>% as.data.frame %>% mutate(step_ix=row_number())
    }
    
    return (results)
} 