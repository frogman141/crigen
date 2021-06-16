#' @export
create_models <- function(crigen_obj) {
    # generating individual ko models for a given target genes grna
    crigen_obj <- crigen_obj %>% fetch_metadata() %>% create_crispr_models()
    return (crigen_obj)
}

fetch_metadata <- function(crigen_obj) {
    
    grna_meta <- crigen_obj$experiment_meta$grna_meta
    cache_dir <- crigen_obj$experiment_params$cache_dir
    dyngen_meta <- grna_meta %>%
                        group_by(grna) %>%
                        do(bind_rows(create_dyngen_meta(.data, cache_dir)))
    
    crigen_obj$experiment_meta$dyngen_meta <- bind_rows(crigen_obj$experiment_meta$dyngen_meta, dyngen_meta)
    return(crigen_obj)
}

create_crispr_models <- function (crigen_obj) {
    nsims <- crigen_obj$simulator_params$num_simulations
    dyngen_meta <- crigen_obj$experiment_meta$dyngen_meta
    
    for(row_i in 1:nrow(dyngen_meta)) {
        dyngen_row <- dyngen_meta[row_i, ]
        sim_df <- duplicate_grnas(dyngen_row, nsims)
        
        expr_params <- simulation_type_knockdown(genes = sim_df$genes,
                                                 num_genes = sim_df$num_genes,
                                                 num_simulations = nsims,
                                                 multiplier = sim_df$kd_multiplier)
        
        # adding caching filepaths to the dyngen model
        dyngen_model <- crigen_obj$simulator_models$common_model
        dyngen_model$model_meta <- dyngen_row
        
        # adding experiment parameters to simulate knockdown
        dyngen_model$simulation_params$experiment_params <- expr_params
        crigen_obj$simulator_models$models[[dyngen_row$model_name]] <- dyngen_model
    }
    
    return (crigen_obj)
}

create_dyngen_meta <- function(grna_meta, cache_dir) {
    # creating metadata for dyngen knock down models. This function works on grna level. As such
    # the metadata provided must grouped by unique grna.
    experiment_type <- unique(grna_meta$experiment_type)
    grna_name <- unique(grna_meta$grna)
    
    if (toupper(experiment_type) == 'KO'){
        dyngen_meta <- create_knockout_meta(grna_meta)
    } else {
        dyngen_meta <- create_other_expr_meta(grna_meta)
    }
    
    dyngen_meta <- add_cache_info(dyngen_meta, grna_name, cache_dir)
    return (dyngen_meta)
}

create_knockout_meta <- function(grna_meta) {
    # create the metadata for simulated crispr knockout experiments
     # generate KD Multiplier matrix. If KO Matrix will be Multiplier rows if interference matrix will have 1 row.
    grna_name <- unique(grna_meta$grna)
    off_target_count <- unique(grna_meta$off_target_count)
    grna_to_ctrl_ratio <- unique(grna_meta$grna_to_ctrl_ratio)
    
    if (off_target_count == 0) {
        meta <- data.frame(model_name = grna_name, sample_percentage = grna_to_ctrl_ratio) %>%
                    mutate(genes = list(grna_meta$gene),
                           num_genes = length(grna_meta$gene),
                           kd_multiplier = list(0))
    } else {
        ko_combos <- knockout_combinations(grna_meta$gene)
        model_names <- name_ko_models(grna_name, ko_combos)
        sample_percentages <- calc_sample_percentages(ko_combos, grna_meta$on_target)
        rownames(ko_combos) <- model_names
        
        meta <- data.frame(model_name = model_names, sample_percentage = (sample_percentages * grna_to_ctrl_ratio)) %>%
                    rowwise %>%
                    mutate(genes = list(grna_meta$gene),
                           num_genes = length(grna_meta$gene),
                           kd_multiplier = kd_multiplier_to_list(model_name, ko_combos))
    }
    
    return (meta)
}

create_other_expr_meta <- function(grna_meta) {
    # create the metadata for simulated crispr interference and activation experiments
    grna_name <- unique(grna_meta$grna)
    grna_to_ctrl_ratio <- unique(grna_meta$grna_to_ctrl_ratio)
    
    meta <- data.frame(model_name = grna_name, sample_percentage = grna_to_ctrl_ratio) %>%
                mutate(genes = list(grna_meta$gene),
                       num_genes = length(grna_meta$gene),
                       kd_multiplier = list(grna_meta$on_target))
    
    return (meta)
}

add_cache_info <- function(meta, grna_name, cache_dir) {
    # adding cache directory and file paths to simulator metadata
    meta <- meta %>% 
                rowwise %>%
                mutate(model_dir = file.path(cache_dir, grna_name),
                       sim_fp = file.path(model_dir, paste0(model_name, '.Rds')),
                       sce_fp = file.path(model_dir,  paste0(model_name, '_sce.Rds')))

    return (meta)
}

knockout_combinations <- function(perturbed_genes) {
    # generating all of the potiential combinations of KO effects for a given gRNA
    # THERE IS AN ISSUE HERE BECAUSE WHEN YOU ONLY HAVE 1 GENE YOU get a weird bug
    combinations <- list()

    for (gene in perturbed_genes){
        combinations[[gene]] = 0:1
    }
    
    ko_combos <- combinations %>% as.data.frame %>% expand.grid()
    return (ko_combos)
}

calc_sample_percentages <- function(ko_combos, on_target){
    # calculate the percentage of cells that should be sampled from 
    # a given model.
    probs_mat <- data.frame()
    
    for (i_row in 1:nrow(ko_combos)) {
        row <- ko_combos[i_row, ]
        
        for (i_col in 1:ncol(ko_combos)) { 
            gene_on_target <- on_target[i_col]
            value <- ko_combos[i_row, i_col]

            if (value == 1) {
                # probability of grna edit
                row[i_col] <- gene_on_target
            } else {
                # probability of no edit
                row[i_col] <- 1 - gene_on_target
            }
        }

        probs_mat <- rbind(probs_mat, row)
    }

    sample_percentages <- rowProds(as.matrix(probs_mat))
    return (sample_percentages)
}

name_ko_models <- function(grna, ko_combos) {
    model_names <- c()

    for (i in 1:nrow(ko_combos)) {
        model_name <- grna

        for (col in colnames(ko_combos)) {
            value <- ko_combos[i, col]
            
            if (value == 1) {
                model_name <- paste0(model_name, "_", col, "PRT")
            } else {
                model_name <- paste0(model_name, "_", col, "NT")
            }
        }
        
        model_names <- c(model_names, model_name)
    }
    
    return (model_names)
}

kd_multiplier_to_list <- function(row_name, df){
    row_values <- c()
    row <- df[row_name, ]
    
    for (col in colnames(row)){
        row_values <- c(row_values, row[, col])
    }
    
    return (list(row_values))
}

duplicate_grnas <- function(dat, sim_count_per_grnas){
    # duplicates the metadata provided by the number of simuluations we want
    # dyngen to run for the given gRNA model
    
    grnas <- data.frame()
    
    for (i in 1:sim_count_per_grnas){
        grnas <- rbind(grnas, dat)
    }
    
    return (grnas)
}