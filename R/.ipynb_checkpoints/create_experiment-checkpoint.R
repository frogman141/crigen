#' @export
create_experiment <- function(crigen_obj){
    # Configuring Dyngen
    crigen_obj$experiment_meta = list()
    
    crigen_obj <- crigen_obj %>%
                        create_cellpops() %>%
                        create_common_model() %>%
                        extract_genes_meta() %>%
                        create_grnas() %>%
                        create_models()
#     crigen_obj$experiment_meta$dyngen_meta <- crigen_obj$experiment_meta$dyngen_meta %>%
#                                                     mutate(sample_percentage=ifelse(grna == 'CTRL', 1, sample_percentage *))
    return (crigen_obj)
}

# configure_dyngen is responsible for configuring the common model that'll be used to simulate perturbations. It takes a population backbone that
# is created by create_cellpops along with some parameters such as num of cells, egenes, and housing keeping (hk) genes, along with other technical parameters.
create_common_model <- function(crigen_obj) {
    # setting up variables that will be used later
    crigen_obj$simulator_models <- list()
    params <- crigen_obj$simulator_params
    cache_dir <- crigen_obj$experiment_params$cache_dir
    ctrl_label <- crigen_obj$experiment_params$ctrl_label
    model_dir <- file.path(cache_dir, ctrl_label)
    ctrl_meta <- create_ctrl_meta(ctrl_label, model_dir)
     
    # create gold standard and simulation params
    gold_standard <- gold_standard_default(census_interval = 1, tau = params$tau)
    sim_params <- simulation_default(census_interval = params$census_interval,
                                     ssa_algorithm = ssa_etl(tau = params$tau),
                                     experiment_params = simulation_type_wild_type(
                                               num_simulations=params$num_simulations
                                     ))
    
    # configure and initialize common dyngen model
    config <- initialise_model(backbone = params$backbone,
                               num_tfs = params$num_tfs,
                               num_hks = params$num_hks,
                               num_cells = params$num_cells,
                               num_targets = params$num_targets,
                               simulation_params = sim_params,
                               gold_standard_params = gold_standard)
    
    # simulate gold standard for common model
    common_model = config %>%
                      generate_tf_network() %>%
                      generate_feature_network() %>% 
                      generate_kinetics() %>%
                      generate_gold_standard()
    
    common_model$model_meta <- ctrl_meta
    crigen_obj$experiment_meta$dyngen_meta <- ctrl_meta
    crigen_obj$simulator_models$common_model <- common_model
    crigen_obj$simulator_models$models[[ctrl_label]] <- common_model
    return (crigen_obj)
}

extract_genes_meta <- function(crigen_obj) {
    # extract and randomly select target genes. This function heavily relays upon dyngen creating the 
    # gene regulatory network.
    expr_params <- crigen_obj$experiment_params
    common_model <- crigen_obj$simulator_models$common_model
    
    burn_tfs <- common_model$feature_info %>% filter(burn == TRUE & is_tf == TRUE & is_hk == FALSE) %>% pull(feature_id) %>% as.character
    grn_network <- common_model$feature_info %>% filter(!feature_id %in% burn_tfs)
    genes <- grn_network %>% pull(feature_id) %>% as.character
    
    if (expr_params$target_tf_only == TRUE) {
        target_genes <- grn_network %>% filter(is_tf == TRUE) %>% pull(feature_id) %>% as.character
    } else {
        target_genes <- genes
    }
    
    if (expr_params$ntargets != 0 & expr_params$ntargets <= length(target_genes)) {
        target_genes <- sample(target_genes, expr_params$ntargets)
    }
    
    crigen_obj$experiment_meta$gene_meta <- list(target_genes=target_genes, genes=genes)
    return (crigen_obj)
}

create_grnas <- function(crigen_obj) {
    # generating grna metadata based on genes metadata
    grna_meta <- data.frame()
    params <- crigen_obj$experiment_params
    gene_meta <- crigen_obj$experiment_meta$gene_meta
    
    for (target_gene in gene_meta$target_genes){
        grna = create_grna(target_gene, gene_meta$genes,
                           params$ngrnas_per_target,
                           params$grna_library,
                           params$experiment_type,
                           params$grna_to_ctrl_ratio)
        
        grna_meta = rbind(grna_meta, grna)
    }
    
    crigen_obj$experiment_meta$grna_meta <- grna_meta
    return (crigen_obj)
}


create_grna <- function(target_gene, genes, num_of_grnas, library, experiment_type, grna_to_ctrl_ratio){
    target_gene_grna <- data.frame()
    potiental_off_target <- genes[target_gene != genes]
    
    for (grna_i in 1:num_of_grnas) {
        off_target_meta <- off_target(potiental_off_target, library)
        grna_name <- paste0(target_gene, '_grna_', grna_i)
        perturb_genes <- c(target_gene, off_target_meta$off_target_genes)
        on_target_activity <- on_target(perturb_genes, off_target_meta$num_of_genes, library, experiment_type)
        
        grna <- data.frame(grna = grna_name,
                          gene = perturb_genes,
                          is_target = FALSE,
                          off_target_count = off_target_meta$num_off_targets,
                          on_target = on_target_activity, 
                          experiment_type = experiment_type,
                          library = library,
                          grna_to_ctrl_ratio = grna_to_ctrl_ratio)
        
        grna <- grna %>% mutate(is_target=ifelse(gene == target_gene, TRUE, is_target))
        target_gene_grna <- rbind(target_gene_grna, grna)
    }
    
    return (target_gene_grna)
}

create_ctrl_meta <- function(ctrl_label, model_dir){
    ctrl_meta <- data.frame(model_name = ctrl_label, sample_percentage = 1, genes = NA,
                        num_genes = 0, kd_multiplier = NA, model_dir = model_dir,
                        sim_fp = file.path(model_dir, paste0(ctrl_label, '.Rds')),
                        sce_fp = file.path(model_dir,  paste0(ctrl_label, '_sce.Rds')),
                        grna = ctrl_label) %>%
                    rowwise %>% mutate(genes=list(c(ctrl_label)), kd_multiplier=list(c(1)))
    
    return (ctrl_meta)
}