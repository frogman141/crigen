# configure_dyngen is responsible for configuring the common model that'll be used to simulate perturbations. It takes a population backbone that
# is created by create_cellpops along with some parameters such as num of cells, egenes, and housing keeping (hk) genes, along with other technical parameters.
#' @export
configure_dyngen <- function(crigen_obj, num_of_cells=1000, num_of_tfs=nrow(backbone$module_info), num_of_egenes=500,
                             num_of_hk=9000, num_of_sims=100, census_interval=10, tau=100/3600) {
    
    # create gold standard and simulation params
    gold_standard <- gold_standard_default(census_interval = 1, tau = tau)
    sim_params <- simulation_default(census_interval = census_interval, ssa_algorithm = ssa_etl(tau = tau),
                                     experiment_params = simulation_type_wild_type(num_simulations = num_of_sims))
    
    # configure and initialize common dyngen model
    config = initialise_model(backbone = backbone,
                              num_tfs = num_of_tfs,
                              num_cells = num_of_cells,
                              num_targets = num_of_egenes,
                              num_hks = num_of_hk,
                              simulation_params = sim_params,
                              gold_standard_params = gold_standard)
    
    # simulate gold standard for common model
    crigen_obj$common_model = config %>%
                                  generate_tf_network() %>%
                                  generate_feature_network() %>% 
                                  generate_kinetics() %>%
                                  generate_gold_standard()

    return (crigen_obj)
}

#' @export
design_experiment <- function(crigen_obj, num_of_targets=0, target_tf_only=TRUE, num_of_grnas=3, library='default', experiment_type='KO'){
    crigen_obj$grna_meta <- data.frame()
    genes_meta <- extract_genes_meta(crigen_obj$common_model, num_of_targets, target_tf_only)
    
    for (target_gene in target_genes){
        grna = create_grna(genes_meta$target_genes, genes_meta$genes, num_of_grnas, library, experiment_type)
        crigen_obj$grna_meta = rbind(crigen_obj$grna_meta, grna)
    }
    
    return (crigen_obj)
}


extract_genes_meta <- function(common_model, num_of_targets, target_tf_only) {
    # extract and randomly select target genes. This function heavily relays upon dyngen creating the 
    # gene regulatory network.
    burn_tfs <- common_model$feature_info %>% filter(burn == TRUE & is_tf == TRUE & is_hk == FALSE) %>% pull(feature_id) %>% as.character
    grn_network <- common_model$feature_info %>% filter(!feature_id %in% burn_tfs)
    genes <- grn_network %>% pull(feature_id) %>% as.character
    
    if (target_tf_only == TRUE) {
        target_genes <- grn_network %>% filter(is_tf == TRUE) %>% pull(feature_id) %>% as.character
    } else {
        target_genes <- genes
    }
    
    if (num_of_targets != 0 & num_of_targets <= length(target_genes)) {
        target_genes <- sample(target_genes, num_of_targets)
    } 
    
    return (list(target_genes=target_genes, genes=genes))
}

create_grna <- function(target_gene, genes, num_of_grnas, library, experiment_type){
    target_gene_grna = data.frame()
    
    for (grna_i in 1:num_of_grnas) {
        
        off_target_meta <- off_target(genes, library)
        grna_name <- paste0(target_gene, '_grna_', grna_i)
        perturb_genes <- c(target_gene, off_target_meta$off_target_genes)
        on_target_activity <- on_target(perturb_genes, off_target_meta$num_of_genes, library, experiment_type)
        
        
        grna <- data.frame(grna = grna_name,
                          gene = perturb_genes,
                          is_target = FALSE,
                          off_target_count = off_target_meta$num_off_targets,
                          on_target = on_target_activity, 
                          experiment_type = experiment_type, 
                          library = library)
        
        grna <- grna %>% mutate(is_target=ifelse(gene == target_gene, TRUE, is_target))
        target_gene_grna <- rbind(target_gene_grna, grna)
    }
    
    return (target_gene_grna)
}