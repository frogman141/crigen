#' @export
simulate_experiment <- function(crigen_obj, slurm=FALSE, ncores=1, max_njobs=100, time='1-00:00:00', mem="8G") {
    # run generate cells for all models in the experiment. Depending on user input we can run the simulations 
    # on a slurm cluster or locally.
    models <- crigen_obj$simulator_models$model
    check_cache(crigen_obj$experiment_meta$dyngen_meta)
    
    
    if (isTRUE(slurm)) {
        crigen_obj <- use_slurm(crigen_obj, max_njobs, time, mem)
    } else {
        crigen_obj <- use_local(crigen_obj, ncores)
    }
    
    return (crigen_obj)
}

check_cache <- function(meta) {
    
    for (row_i in 1:nrow(meta)) {
        row <- meta[row_i, ]
        
        if (!dir.exists(row$model_dir)) {
            dir.create(row$model_dir, recursive=TRUE)
        }
    }
}

use_local = function(crigen_obj, ncores) {
    ndetected <- detectCores()
    
    if(ncores > ndetected){
        message("ncores specified greater than the number of cores detected...")
        message("Using all cores avaliable instead...")
        ncores <- ndetected
    }
    
    sces <- mclapply(crigen_obj$simulator_models$models, run_dyngen, mc.cores=ncores)
    crigen_obj$sce <- merge_sces(sces)
    return (crigen_obj)
}

use_slurm <- function(crigen_obj, max_njobs, time, mem) {
    # submitting a job array to slurm to parallel process
    # setting dyngen job array job count, sbatch options and submitting array to slurm
    models <- crigen_obj$simulator_models$models
    nmodels <- length(models)
    job_opts <- list(time=time, mem=mem)
    
    # setting the number of jobs to be submitted to slurm
    if (nmodels > max_njobs) {
        njobs <- max_njobs
    } else {
        njobs <- nmodels
    }
    
    sces <- Slurm_lapply(models, run_dyngen, mc.cores=1, njobs=njobs, sbatch_opt=job_opts, plan="collect")
    crigen_obj$sce <- merge_sces(sces)
    return  (crigen_obj)
}

run_dyngen <- function(dyngen_model) {
    
    update_cell_metadata <- function(sim_model) {
        # update cell metadata with the gRNA/Perturbation Data
        grna <- sim_model$model_meta$grna
        model_name <- sim_model$model_meta$model_name
        target_gene <- sim_model$model_meta$genes[[1]][1]
        
        sim_model$experiment <- rename_cells_index(sim_model$experiment, model_name)

        # rename cell index and extracting which cell population a cell belongs to
        cell_meta <- sim_model$experiment$cell_info %>%
                                        mutate(grna_name = grna, 
                                               model = model_name,
                                               target_gene = target_gene,
                                               cell_population = gsub('^.|.$', '', from))
        print (cell_meta)
        sim_model$experiment$cell_info <- cell_meta
        return (sim_model)
    }

    rename_cells_index <- function(experiment, model_name) {
        old_cell_id <- experiment$cell_info$cell_id
        new_cell_ids <- paste0(old_cell_id , "_", model_name)

        experiment$cell_info$cell_id <- new_cell_ids
        rownames(experiment$counts_mrna) <- new_cell_ids
        rownames(experiment$counts_premrna) <- new_cell_ids
        rownames(experiment$counts_protein) <- new_cell_ids

        return (experiment)
    }
    
    # generate cells and experiments
    print("Using Dyngen to Generate Experiments...")
    dyngen_model = dyngen_model %>% generate_cells() %>% generate_experiment()

    # update cell metadata from dyngen 
    print("Update Cell Metadata...")
    dyngen_model <- update_cell_metadata(dyngen_model)
    
    print ("Convert to SCE...")
    sce <- as_sce(dyngen_model)
    
    # save dyngen
    print("Cache Simulation Results...")
    saveRDS(sce, dyngen_model$model_meta$sce_fp)
    saveRDS(dyngen_model, dyngen_model$model_meta$sim_fp)
    
    print ("Downsampling Results...")
    sampled_cells <- downsample_cells(colnames(sce), dyngen_model$model_meta$sample_percentage)
    sce <- sce[, sampled_cells]
    
    return (sce)
}