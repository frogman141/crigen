#' @export
create_models <- function(crigen_obj, ctrl_label='CTRL', cache_dir='../data') {
    # generating individual ko models for a given target genes grna
    crigen_obj <- setup(crigen_obj, ctrl_label, cache_dir)
    crigen_obj <- fetch_metadata(crigen_obj)
    crigen_obj <- create_crispr_models(crigen_obj)
    
    return (crigen_obj)
}

setup <- function(crigen_obj, ctrl_label, cache_dir) {
    crigen_obj$models <- list()
    crigen_obj$cache_dir <- cache_dir
    crigen_obj$ctrl_label <- ctrl_label
    crigen_obj$models[[ctrl_label]] <- crigen_obj$common_model
    
    return (crigen_obj)
}

fetch_metadata <- function(crigen_obj) {
    crigen_obj$dyngen_meta <- crigen_obj$grna_meta %>%
                            group_by(grna) %>%
                            do(bind_rows(create_dyngen_meta(.data, crigen_obj$cache_dir)))
    
    return(crigen_obj)
}

create_crispr_models <- function (crigen_obj) {
    nsims <- nrow(crigen_obj$common_model$simulation_params$experiment_params)
    
    for(row_i in 1:nrow(crigen_obj$dyngen_meta)) {
        dyngen_model <- crigen_obj$common_model
        dyngen_crispr_meta <- crigen_obj$dyngen_meta[row_i, ]
        sim_df <- duplicate_grnas(dyngen_crispr_meta, nsims)
        
        expr_params <- simulation_type_knockdown(genes = sim_df$genes,
                                                 num_genes = sim_df$num_genes,
                                                 num_simulations = nsims,
                                                 multiplier = sim_df$kd_multiplier)

        dyngen_model$simulation_params$experiment_params <- expr_params
        crigen_obj$models[[dyngen_kd_meta$model_name]] <- dyngen_model
    }
    
    return (crigen_obj)
}

create_dyngen_meta <- function(meta, cache_dir) {
    # creating metadata for dyngen knock down models. This function works on grna level. As such
    # the metadata provided must grouped by unique grna.
    experiment_type <- unique(meta$experiment_type)
    
    if (experiment_type == 'KO'){
        dyngen_meta <- create_knockout_meta(meta)
    } else {
        dyngen_meta <- create_interfernce_or_activation_meta(meta)
    }
    
    dyngen_meta <- add_cache_info(dyngen_meta)
    return (dyngen_meta)
}

create_knockout_meta <- function(meta) {
    # create the metadata for simulated crispr knockout experiments
     # generate KD Multiplier matrix. If KO Matrix will be Multiplier rows if interference matrix will have 1 row.
    grna_name <- unique(meta$grna)
    ko_combos <- knockout_combinations(meta$gene)
    model_names <- name_ko_models(grna_name, ko_combos)
    sample_percentages <- calc_sample_percentages(ko_combos, meta$on_target)
    
    meta <- data.frame(model_name = model_names,
                       sample_percentage = sample_percentages) %>%
                rowwise %>%
                mutate(genes = list(meta$gene),
                       num_genes = length(meta$gene),
                       kd_multiplier = list(ko_combos[model_name, ]))
    
    return (meta)
}

create_interfernce_or_activation_meta <- function(meta) {
    # create the metadata for simulated crispr interference and activation experiments
    grna_name <- unique(meta$grna)
    
    meta <- data.frame(model_name = grna_name,
                       sample_percentage = 1,
                       genes = NA) %>%
                mutate(genes = list(meta$gene),
                       num_genes = length(meta$gene),
                       kd_multiplier = list(meta$on_target))
    
    return (meta)
}

add_cache_info <- function(meta) {
    # adding cache directory and file paths to simulator metadata
    meta <- meta %>% 
                rowwise %>%
                mutate(model_dir = file.path(cache_dir, model_name),
                       sim_fp = file.path(cache_dir, paste0(model_name, '.Rds')),
                       sce_fp = file.path(cahce_dir,  paste0(model_name, '_sce.Rds')))

    return (meta)
}

knockout_combinations <- function(perturbed_genes) {
    # generating all of the potiential combinations of KO effects for a given gRNA
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

duplicate_grnas <- function(dat, sim_count_per_grnas){
    # duplicates the metadata provided by the number of simuluations we want
    # dyngen to run for the given gRNA model
    
    grnas <- data.frame()
    
    for (i in 1:sim_count_per_grnas){
        grnas <- rbind(grnas, dat)
    }
    
    return (grnas)
}