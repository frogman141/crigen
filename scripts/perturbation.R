############################# Wrapper Function to Generate Perturbations of Dyngens GRN #############################

run_perturbation_experiment = function(common_model, target_genes, off_target_genes=c(),
                                      num_of_grnas=1, off_target_vec=c(0), random_kd_multiplier=FALSE,
                                      kd_multiplier=0, num_of_sims=100, timepoint=0, ctrl_label="CTRL", cluster=FALSE){
    
    target_genes = c(ctrl_label, target_genes)
    
    # creating perturbation models, generate cells for each model, generate experiment, and update cell metadata
    perturb_models = create_perturbation_models(common_model, target_genes,  off_target_genes, num_of_grnas, off_target_vec,
                                                random_kd_multiplier, kd_multiplier, num_of_sims, timepoint, ctrl_label)
    
    perturb_models = generate_models_cells(perturb_models, cluster)
    model_comb = combine_models(perturb_models) %>% generate_experiment()
    model_comb = update_cell_metadata(model_comb)

    return (model_comb)
}

############################# Wrapper Functions to Convert Dyngen Object to SCE #############################

convert_dyngen_to_sce = function(model){
    # take the provided dyngen model (must have all perturbations merged) and convert it to a SCE object
    sce = as_sce(model)
    sce = update_sce_coldata(sce, model)
    
    return (sce)
}

update_sce_coldata = function(sce, model) {
    # mergining metadata from dyngen into sce object. That were previously left out such as target_gene, grna_name, and cells pop.
    
    cell_info = model$experiment$cell_info %>% select(cell_id, target_gene, grna_name, cell_population)
    colData(sce)$cell_id = rownames(colData(sce))
    
    temp = colData(sce) %>%
            as.data.frame %>%    
            left_join(cell_info, by='cell_id') %>% 
            select(cell_id, target_gene, grna_name, cell_population)

    colData(sce) = cbind(colData(sce), temp)
    return(sce)
}

############################# Utility Function for Creating gRNA Level perturbations of Dyngens GRN #############################

create_perturbation_models = function(common_model, target_genes,  off_target_genes,
                                      num_of_grnas, off_target_vec, random_kd_multiplier,
                                      kd_multiplier, num_of_sims, timepoint, ctrl_label) {
    # creating individual instances of both CTRL and target gene gRNA KO models
    models = list()

    for (target_gene in target_genes) {
        if (target_gene == 'CTRL'){
            models[[target_gene]] = common_model

        } else {
            # create KO models for each Target Gene gRNAs
            genes_grnas = create_gene_grnas(target_gene, off_target_genes, num_of_grnas,
                                            off_target_vec, random_kd_multiplier, kd_multiplier)
            grna_models = create_grna_models(common_model, genes_grnas, num_of_sims, timepoint)
            models = c(models, grna_models)
        }
    }
    
    return (models)
}

generate_models_cells = function(perturb_models, cluster) {
    # run generate cells for all models in the experiment. Depending on user input we can run the simulations 
    # on a slurm cluster or locally.
    if (cluster == TRUE) {
        perturb_models = run_slurm(perturb_models)
    } else {
        perturb_models = run_locally(perturbed_models)
    }
    
    return (perturb_models)
}

create_gene_grnas = function(target_gene, off_target_genes, num_of_grnas, off_target_vec, random_kd_multiplier, kd_multiplier) {
    
    target_grnas = data.frame()
    not_target_genes = off_target_genes[off_target_genes != target_gene]
    
    for (grna_i in 1:num_of_grnas){
        # for each grna get the off and on target activity
        off_target = off_target_effect(off_target_vec, not_target_genes)
        
        grna = data.frame(grna_name = paste0(target_gene, "_grna_", grna_i),
                          perturbed_genes = paste(c(target_gene, off_target$off_target_genes), collapse=','),
                          num_genes = off_target$num_of_genes,
                          num_off_target =  off_target$num_off_targets)
        
        target_grnas = rbind(target_grnas, grna)
    }
    
    target_grnas = target_grnas %>%
                        rowwise %>%
                        mutate(target_gene = target_gene, 
                               perturbed_gene_list = strsplit(perturbed_genes, ','),
                               kd_multiplier = on_target_effect(kd_multiplier, random_kd_multiplier, num_genes)) %>%
                        select(grna_name, target_gene, num_genes, num_off_target, kd_multiplier, perturbed_gene_list)
    
    return (target_grnas)
}

create_grna_models = function(common_model, gene_grnas, num_of_sims, timepoint) {
    # generating individual ko models for a given target genes grna
    grna_models = list()
    grna_names = gene_grnas$grna_name

    for (grna_i in 1:nrow(gene_grnas)) {
        model_grna = common_model
        grna = gene_grnas[grna_i, ]
        grna_name = grna_names[grna_i]
        grna = duplicate_grnas(grna, num_of_sims)

        expr_params = simulation_type_knockdown(num_simulations = num_of_sims,
                                                timepoint = timepoint, 
                                                num_genes = grna$num_genes,
                                                genes = grna$perturbed_gene_list,
                                                multiplier = grna$kd_multiplier)

        model_grna$simulation_params$experiment_params = expr_params
        grna_models[[grna_name]] = model_grna
    }
    
    return (grna_models)
}

update_cell_metadata = function(comb_model) {
    # update cell metadata with the gRNA/Perturbation Data
    cell_meta = comb_model$experiment$cell_info
    kd_multiplier = comb_model$simulations$kd_multiplier

    cell_meta = cell_meta %>%
                    left_join(kd_multiplier, by='simulation_i') %>%
                    mutate(grna_name = model,
                           cell_population=str_extract(grna_name, ".+?(?=_)"), 
                           target_gene=gsub('*|_gRNA_[1-9]', '', model),
                           multiplier=ifelse(is.na(multiplier), 1, multiplier))

    comb_model$experiment$cell_info = cell_meta
    return (comb_model)
}

duplicate_grnas = function(dat, sim_count_per_grnas){
    grnas = data.frame()
    
    for (i in 1:sim_count_per_grnas){
        grnas = rbind(grnas, dat)
    }
    
    return (grnas)
}

off_target_effect = function(off_target_vec, genes){
    #select the number off target genes through random sampling
    num_off_targets = sample(off_target_vec, 1, replace=TRUE)
    num_of_genes = num_off_targets + 1
    
    # randomly select the genes that are going to be off targets
    off_target_genes = sample(genes, num_off_targets)
    return (list(num_off_targets=num_off_targets, off_target_genes=off_target_genes, num_of_genes=num_of_genes))
}

on_target_effect = function(kd_multiplier, random_kd_multiplier, num_of_targets){
    # randomly samples the on target activity of a gRNA for both it's target genes and off-target genes
        
    if (random_kd_multiplier == TRUE){
        on_target = runif(num_of_targets)
    } else {
        on_target = rep(kd_multiplier, num_of_targets / length(kd_multiplier))
    } 
    
    return (list(on_target))
}

run_slurm = function(perturb_models){
    # submitting a job array to slurm to parallel process
    # setting dyngen job array job count, sbatch options and submitting array to slurm
    num_of_models = length(perturb_models)
    dyngen_job_opts = list(time='4-00:00:00', mem="4G")
    perturb_models = Slurm_lapply(perturb_models, run_generate_cells, mc.cores=1,
                                  njobs=num_of_models, sbatch_opt=dyngen_job_opts, plan="collect")
    
    return  (perturb_models)
}

run_locally = function(perturb_models) {
    # code to generate cells of dyngen models locally
    for (model_name in names(perturb_models)) {
        message(paste0("Generating Cells for ", model_name, "..."))
        
        # select model and generate cells for said model
        model = perturb_models[[model_name]]
        model = run_generate_cell(model)

        # reassign model with the generated cells
        perturb_models[[model_name]] = model
    }
    
    return (perturb_models)
}

run_generate_cells = function(model) {
    model = model %>% generate_cells()
    return (model)
}