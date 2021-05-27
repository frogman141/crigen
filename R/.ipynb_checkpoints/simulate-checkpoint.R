#' @export
run_simulation <- function(models, slurm=FALSE, ncores=1, max_njobs=200, time='1-00:00:00', mem="6G") {
    # run generate cells for all models in the experiment. Depending on user input we can run the simulations 
    # on a slurm cluster or locally.
    if (slurm == TRUE) {
        models_fp <- run_slurm(models, max_njobs, time, mem)
    } else {
        models_fp <- run_locally(models, ncores)
    }
    
    return (models_fp)
}

run_locally = function(models, ncores) {
    ndetected = detectCores()
    
    if(ncores > ndetected){
        message("ncores specified greater than the number of cores detected...")
        message("Using all cores avaliable instead")
        ncores <- ndetected
    }
    
    models = mclapply(models, generate_cells, mc.cores=ncores)
    return (models)
}


run_slurm = function(models, max_njobs, time, mem) {
    # submitting a job array to slurm to parallel process
    # setting dyngen job array job count, sbatch options and submitting array to slurm
    nmodels = length(models)
    job_opts = list(time=time, mem=mem)
    
    # setting the number of jobs to be submitted to slurm
    if (nmodels > max_njobs) {
        njobs = max_njobs
    } else {
        njobs = nmodels
    }
    
    models = Slurm_lapply(models, generate_cells, mc.cores = 1, njobs = njobs, sbatch_opt = job_opts, plan = "collect")
    return  (models)
}

generate_cells = function(model) {
    model_fp = model$cache_fp
    cache_dir = model$cache_dir
    
    # generate cells and experiments
    message("Using Dyngen to Generate Experiments...")
    model = model %>% 
                generate_cells() %>%
                generate_experiment()
    
    # update cell metadata from dyngen 
    message("Update Cell Metadata and Convert to SCE...")
    model = update_cell_metadata(model)
    sce = convert_dyngen_to_sce(model)
    
    # save dyngen and sce objects
    message("Cache Simulation Results...")
    save_output(model, model$cache_fp, model$sce_fp)
    save_output(model, model$cache_fp, model$dyngen_fp)
    return (list(sce_fp=model$sce_fp, dyngen_fp=model$dyngen_fp))
}

update_cell_metadata = function(model) {
    # update cell metadata with the gRNA/Perturbation Data
    cell_meta = model$experiment$cell_info
    kd_multiplier = model$simulations$kd_multiplier
    
    # extracting which cell population a cell comes from along with what it's gRNAs target gene is.
    cell_meta = cell_meta %>%
                    left_join(kd_multiplier, by='simulation_i') %>%
                    mutate(grna_name = model,
                           cell_population = str_extract(model, ".+?(?=_)"), 
                           target_gene = gsub('*|_gRNA_[1-9]', '', model),
                           multiplier = ifelse(is.na(multiplier), 1, multiplier))

    model$experiment$cell_info = cell_meta
    return (model)
}

convert_dyngen_to_sce = function(model){
    # take the provided dyngen model (must have all perturbations merged) and convert it to a SCE object
    sce = as_sce(model)
    sce = update_sce_coldata(sce, model)
    
    return (sce)
}

save_output = function(obj, cache_dir, cache_fp) {
    
    if (!dir.exists(cache_dir)) {
        dir.create(cache_dir, recursive=TRUE)
    }
    
    saveRDS(obj, cache_fp)
}

update_sce_coldata = function(sce, model, cols_of_interest=c('cell_id', 'target_gene', 'grna_name', 'cell_population')) {
    # mergining metadata from dyngen into sce object. That were previously left out such as target_gene, grna_name, and cells pop.
    
    cell_info = model$experiment$cell_info
    colData(sce)$cell_id = rownames(colData(sce))
    
    temp = colData(sce) %>%
            as.data.frame %>%    
            left_join(cell_info[, cols_of_interest], by='cell_id')

    colData(sce) = cbind(colData(sce), temp)
    return(sce)
}