#' @export
# add ctrl grnas to the datasets
create_grnas <- function(model,
                         ntargets=NA,
                         on_target=NA,
                         off_target=NA,
                         ctrl_label='CTRL',
                         ngrna_per_target=3,
                         target_tf_only=TRUE,
                         crispr_type='Interference',
                         library=sample(names(grna_libraries_meta), size=1)) {
    
    assert_that(crispr_type %in% c('Interference', 'Activation', 'Knockout'),
                msg=paste0("Crispr Type provided not: Interference, Activation, or Knockout"))
    
    model$simulation_system$crispr_params <- list(crispr_type=crispr_type, on_target=on_target,
                                                  off_target=off_target, ctrl_label=ctrl_label)
    
    # get the number of targets and select the targets  
    message("Selecting Target Genes...")
    ntargets <- get_ntargets(model, ntargets, target_tf_only)
    model$target_genes <- select_target_genes(model, ntargets, target_tf_only)
    
    message("Checking gRNA Parameters...")
    pdf <- grna_libraries_meta[[library]]
    grna_params <- get_grna_params(pdf)
    on_target <- get_on_target(on_target, grna_params)
    off_target <- get_off_target(off_target, grna_params)
    
    # generate grna, multiplier_df, and model metadata
    message("Creating gRNA Characteristics...")
    model <- model %>%
                generate_grnas(ngrna_per_target, on_target, off_target, pdf) %>%
                generate_multiplier(crispr_type) %>%
                generate_simulation_meta()
    
    return (model)
}

get_grna_params <- function(pdf) {
    row_i <- sample(1:nrow(pdf), size=1, prob=pdf$probability)
    grna_params <- pdf[row_i, ]
    
    return (grna_params)
}

get_on_target <- function(on_target, grna_params) {
    if (is.na(on_target)) { 
        on_target <- (1 - (grna_params$on_target / 100))
    }
    
    return (on_target)
}

get_off_target <- function(off_target, grna_params) {
    if (is.na(off_target)) { 
        off_target <- grna_params$off_target
    }
    
    return (off_target)
}

get_ntargets <- function(model, ntargets, target_tf_only) {
    # number of potential targets
    feature_info <- model$feature_info
    
    if (isTRUE(target_tf_only)) {
       npotential_targets <- feature_info %>% filter(is_tf == TRUE & is_hk == FALSE) %>% nrow
    } else {
       npotential_targets <- feature_info %>% nrow
    }
    
    # get the number of target genes 
    if (is.na(ntargets)) {
        if (isTRUE(target_tf_only)) {
            ntargets <- feature_info %>% filter(is_tf == TRUE & is_hk == FALSE) %>% nrow
        } else {
            ntargets <- feature_info %>% nrow
        }
    }
    
    assert_that(
        ntargets <= npotential_targets,
        msg = paste0("Number of Target Genes Provided (", ntargets, ") is too large enough (>= number of avaliable genes ", npotential_targets, ")")
    )
    
    return (ntargets)
}

select_target_genes <- function(model, ntargets, target_tf_only) {
    # extract and randomly select target genes. This function heavily relays upon dyngen creating the 
    # gene regulatory network.
    feature_info <- model$feature_info
    
    if (target_tf_only == TRUE) {
        target_genes <- feature_info %>% filter(is_tf == TRUE) %>% pull(feature_id) %>% as.character
    } else {
        target_genes <- feature_info %>% pull(feature_id) %>% as.character
    }
    
    if (ntargets < length(target_genes)) {
        target_genes <- sample(target_genes, size=ntargets)
    }
    
    return (target_genes)
}

generate_grnas <- function(model, ngrna_per_target, on_target, off_target, pdf) {
    grna_seq <- 1:ngrna_per_target
    genes <- model$feature_info$feature_id
    grna_names <- lapply(model$target_genes, function (x) paste0(x, "-grna.", grna_seq)) %>% unlist()
    model$grnas_meta <- data.frame(grna=model$simulation_system$crispr_params$ctrl_label, gene=NA, is_target=NA, on_target_activity=NA)
                         
    for (grna in grna_names) {
        target_gene <- str_split(grna, '-')[[1]][1]
        perturbed_genes <- c(target_gene, sample(genes, size=off_target))
        genes_on_target <- sample_on_targets(perturbed_genes, target_gene, on_target, pdf)
        
        meta <- data.frame(grna=rep(grna, length(perturbed_genes)),
                           gene=perturbed_genes,
                           is_target=ifelse(perturbed_genes == target_gene, TRUE, FALSE),
                           on_target_activity=genes_on_target)
        
        model$grnas_meta <- rbind(model$grnas_meta, meta)
    }
    
    return (model)
}

generate_multiplier <- function(model, crispr_type) { 
    multiplier <- NULL
    genes <- model$feature_info$feature_id
    base_multi <- set_names(rep(1, length(genes)), genes)
    grnas_names <- model$grnas_meta %>% pull(grna) %>% unique
      
    for (grna_name in grnas_names) {
        
        if (grna_name == model$simulation_system$crispr_params$ctrl_label) {
             grna_multi <- base_multi %>% t %>% as.data.frame
             rownames(grna_multi) <- c(model$simulation_system$crispr_params$ctrl_label)
        } else {
            grna_meta <- model$grnas_meta %>% filter(grna == grna_name)
            prtb_multi <- get_grna_multiplier(grna_meta, crispr_type)
            grna_multi <- update_grna_multiplier(prtb_multi, base_multi)
        }
        
        multiplier <- rbind(multiplier, grna_multi)
    }
    
    colnames(multiplier) <- paste0("transcription_rate_", colnames(multiplier))
    model$simulation_system$multiplier <- multiplier
    return (model)
}
                         
generate_simulation_meta <- function(model) {
    
    crispr_type <- model$simulation_system$crispr_params$crispr_type
    
    if (crispr_type != 'Knockout') {
        sample_percentages <- rep(1, nrow(model$simulation_system$multiplier))
    } else {
        sample_percentages <- calculate_sample_percentages(model)
    }
    
    model$sim_meta <- data.frame(sim_name=rownames(model$simulation_system$multiplier),
                                 sample_percentage=sample_percentages) %>%
                            rowwise %>%
                            mutate(grna=str_split(sim_name, "_")[[1]][1],
                                   grna=ifelse('' == grna, sim_name, grna)) %>%
                            filter(sample_percentage != 0)
    
    return (model)
}
                               
sample_on_targets <- function(perturbed_genes, target_gene, on_target, pdf) {
    genes_on_target <- c()
    
    for (gene in perturbed_genes) {
        if (gene == target_gene & !is.na(on_target)) {
            sampled_activity <- on_target
        } else {
            # if user provided on target activity is na sample it
            sampled_activity <- sample(pdf$on_target, size = 1, prob = pdf$probability)
            sampled_activity  <- (1 - (sampled_activity / 100))
        }
        
        genes_on_target <- c(genes_on_target, sampled_activity)
    }
    
    return (genes_on_target)
}

get_grna_multiplier <- function(grna, crispr_type) { 
    ngenes <- grna %>% nrow
    grna_name <- grna %>% pull(grna) %>% unique
    
    if (crispr_type == 'Knockout' & ngenes > 1) {
        grna_multipliers <- all_ko_combinations(grna$gene)
        model_names <- name_ko_models(grna_name, grna_multipliers)

    } else if (crispr_type == 'Knockout' & ngenes == 1) {
        grna_multipliers <- data.frame(gene=c(0, 1))
        colnames(grna_multipliers) <- grna$gene
        
        model_names <- c(paste0(grna_name, '-', grna$gene, '_PRT'),
                         paste0(grna_name, '-', grna$gene, '_NT'))
    } else {
        grna_multipliers <- grna %>% select(on_target_activity) %>% t()
        colnames(grna_multipliers) <- grna %>% pull(gene)
        model_names <- grna_name
    }

    rownames(grna_multipliers) <- model_names
    return (grna_multipliers)
}
                         
update_grna_multiplier <- function (prtb_multi, base_multi) {
    multi <- NULL
    
    for (row in rownames(prtb_multi)) {
        sim_multi <- base_multi
        
        for (col in colnames(prtb_multi)) {
            sim_multi[[col]] <- prtb_multi[row, col]
        }
        
        sim_multi <- sim_multi %>% t %>% as.data.frame
        rownames(sim_multi) <- row
        
        multi <- rbind(multi, sim_multi)
    }
    
    return (multi)
}
                         
all_ko_combinations <- function(perturbed_genes) {
    # generating all of the potiential combinations of KO effects for a given gRNA
    combinations = list()

    for (gene in perturbed_genes) {
        combinations[[gene]] = 0:1
    }
    
    ko_combinations = combinations %>% as.data.frame %>% expand.grid()
    return (ko_combinations)
}

name_ko_models <- function(grna, ko_combos) {
    model_names <- c()

    for (i in 1:nrow(ko_combos)) {
        model_name <- grna

        for (col in colnames(ko_combos)) {
            value <- ko_combos[i, col]
            spacer <- ifelse(model_name == grna, '-', '_')
            model_name <- paste0(model_name, spacer, col)
            
            if (value == 0) {
                model_name <- paste0(model_name, "_PRT")
            } else {
                model_name <- paste0(model_name, "_NT")
            }
        }
        
        model_names <- c(model_names, model_name)
    }
    
    return (model_names)
}
                         
calculate_sample_percentages <- function(model){
    sample_percentages <- c()
    multiplier <- model$simulation_system$multiplier
    
    for (sim_name in rownames(multiplier)) {
        
        if (sim_name == model$simulation_system$crispr_params$ctrl_label) {
            sample_percentages <- c(sample_percentages, 1)
            next()
        }
        
        grna_name <- paste(str_split(sim_name, "-")[[1]][1:2], collapse='-')
        grna <- model$grnas_meta %>% filter(grna == grna_name)
        ngenes <- grna %>% nrow
        
        if (ngenes == 1) {
            on_target <- grna$on_target_activity
            sim_perc <- ifelse(grepl('PRT', sim_name), 1 - on_target, on_target)
        } else {
            sim_perc <- multiple_genes_percentage(grna, multiplier, sim_name)
        }
        
        sample_percentages <- c(sample_percentages, sim_perc)
    }
    
    return (sample_percentages)
}
                         
multiple_genes_percentage <- function (grna_meta, multiplier, sim_name) {
    genes <- paste0("transcription_rate_", grna_meta$gene)
    grna_meta$gene <- genes
    ko_combos <- multiplier[sim_name, genes]

    for (col_name in colnames(ko_combos)) { 
        gene_on_target <- grna_meta %>%
                            filter(col_name == gene) %>%
                            pull(on_target_activity)

        value <- ko_combos[sim_name, col_name]

        if (value == 1) {
            # probability of grna edit
            ko_combos[, col_name] <- gene_on_target
        } else {
            # probability of no edit
            ko_combos[, col_name] <- 1 - gene_on_target
        }
    }

    sim_perc <- matrixStats::rowProds(as.matrix(ko_combos))
    return (sim_perc)
}