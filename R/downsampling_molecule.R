#' @export
downsampling_molecule <- function(model, ncells=NA, ncells_per_grna=100, map_reference_ls=TRUE, map_reference_cpm=TRUE) {
    realcount <- model %>% get_realcounts()
    check_experiment_parameters(model, ncells, ncells_per_grna)
    
    # separate out mRNA and pre-mRNA
    spliced_mrna <- sim_counts[, model$feature_info$mol_mrna, drop = FALSE]
    unspliced_mrna <- sim_counts[, model$feature_info$mol_premrna, drop = FALSE]
    
    # rename columns by genes and rownames by cell
    dimnames(spliced_mrna) <- list(cell_info$cell_id, model$feature_info$feature_id)
    dimnames(unspliced_mrna) <- list(cell_info$cell_id, model$feature_info$feature_id)
    
    # add together unspliced and spliced mrna
    sim_counts <- unspliced_mrna + spliced_mrna
    sim_counts <- simulate_counts(realcount, sim_counts, map_reference_ls, map_reference_cpm)
    
    cell_info <- sample_cells(model, ncells, ncells_per_grna, model$simulation$sampling$meta)
    model$experiment <- list(meta=cell_info, counts=sim_counts)
    
    return (model)
}

#' @export
check_experiment_parameters <- function(model, ncells, ncells_per_grna) {
    
    meta <- model$simulation$sampling$meta
    
    if (is.na(ncells) & is.na(ncells_per_grna)) {
        stop("Niether parameter ncells and ncells_per_grna have been defined. \n",
             "Chose one of the parameters and provide an integer.")
    } else if (is.integer(ncells) & is.integer(ncells_per_grna)) {
        stop("Both parameter ncells and ncells_per_grna have been defined. \n",
             "Chose one of the parameters and provide an integer.")
    } else if (is.numeric(ncells)) {
        num_cells <- ncells
    } else if (is.numeric(ncells_per_grna)){
        num_cells <- ncells_per_grna * length(unique(meta$grna))
    }
    
    if (num_cells > nrow(meta)) {
        stop(
          "Simulations only generated ", nrow(meta), " different datapoints, whereas ", num_cells, " cells were requested.\n",
          "Increase the number of simulations or decreasing the census interval (see `?run_simulation`)."
        )
    }
    
}

#' @export
get_realcounts <- function(model) {
    dataset <- sample(realcounts$dataset_id, size=1)
    url <- realcounts %>% filter(dataset_id == dataset) %>% pull(links.download)
    realcount <- download_cacheable_file(url, model$simulation_system$cache_dir, TRUE)
    
    assert_that('dgCMatrix' %in% class(realcount),
                msg='Provided Realcount is not a matrix. Please provide a single cell expression matrix.')
    
    return (realcount)
}

#' @export
simulate_counts <- function(realcount, tsim_counts, map_reference_ls, map_reference_cpm) {
  
  lib_size <- calculate_library_size(realcount, tsim_counts, map_reference_ls)
  tsim_counts_cpm <- calculate_cpm(realcount, tsim_counts, map_reference_cpm)
  
  # simulate sampling of molecules
  tsim_counts_t <- sample_molecules(tsim_counts_cpm, lib_size)
  sim_counts <- Matrix::drop0(Matrix::t(tsim_counts_t))
  
  return (sim_counts)
}

calculate_library_size <- function(realcount, tsim_counts, map_reference_ls) {
  realcount_ls <- Matrix::rowSums(realcount)
  tsim_counts_ls <- Matrix::rowSums(tsim_counts)
  
  # simulate library size variation from real data
  if (map_reference_ls) {
    lib_size <- round(stats::quantile(realcount_ls, sort(stats::runif(length(tsim_counts_ls)))[order(order(tsim_counts_ls))])) %>% unname
  } else {
    lib_size <- tsim_counts_ls
  }
  
  return (lib_size)
}

calculate_cpm <- function(realcount, tsim_counts, map_reference_cpm) {
    realcount_cpm <- realcount
    tsim_counts_cpm <- tsim_counts
    realcount_ls <- Matrix::rowSums(realcount)
    tsim_counts_ls <- Matrix::rowSums(tsim_counts)
    realcount_cpm@x <- realcount@x / realcount_ls[realcount@i+1] * 1e6
    tsim_counts_cpm@x <- tsim_counts@x / tsim_counts_ls[tsim_counts@i+1] * 1e6

    # map real density on tsim counts
    tsim_counts_cpm_new <- tsim_counts_cpm

    if (map_reference_cpm) {
        tsim_counts_cpm_new@x <- stats::quantile(realcount_cpm@x, sort(stats::runif(length(tsim_counts_cpm@x)))[order(order(tsim_counts_cpm@x))]) %>% unname
    }
  
  return (tsim_counts_cpm_new)
}

sample_molecules <- function(tsim_counts_cpm, lib_size_new) {
  tsim_counts_t <- Matrix::t(tsim_counts_cpm)
  
  new_vals <- unlist(map(seq_along(lib_size_new), function(cell_i) {
    pi <- tsim_counts_t@p[[cell_i]]
    pj <- tsim_counts_t@p[[cell_i + 1]]
    
    if (pi != pj) {
      pix <- seq(pi + 1, pj, 1)
      gene_is <- tsim_counts_t@i[pix] + 1
      gene_vals <- tsim_counts_t@x[pix]
      lib_size <- lib_size_new[[cell_i]]
      
      # sample 'lib_size' molecules for each of the genes, weighted by 'gene_vals'
      rmultinom(1, lib_size, gene_vals)
      
    } else {
      integer(0)
    }
  })) %>% as.numeric
  
  tsim_counts_t@x <- new_vals
  return (tsim_counts_t)
}


###################################### Sampling Cell Code ######################################

# NEED TO FIX CELL SAMPLING FOR CRISPR KNOCKOUT
#' @export
sample_cells <- function(model, ncells, ncells_per_grna, meta) {
    # need to add a test to make sure i have enough simulations to sample from
    if (!is.na(ncells)) {
        cell_info <- sampling_cells(meta, ncells)
    }

    if (!is.na(ncells_per_grna)) {
        cell_info <- sampling_cells_per_grna(model, meta, ncells_per_grna)
    }
    
    cell_info <-  cell_info %>%
                      mutate(target_perturbed = TRUE, cell_id = paste0("cell", row_number())) %>%
                      rowwise %>% 
                      mutate(target_perturbed=ifelse(crispr_type != 'Knockout', TRUE, 
                                              ifelse(grepl(paste0(target_gene, '_PRT'), sim_name), TRUE, FALSE))) %>%
                      select(cell_id, step_ix, target_gene, target_perturbed, grna, sim_name, phase, time) %>%
                      as.data.frame

    return (cell_info)
}

sampling_cells <- function(meta, num_cells) {
    sampled_cells <- NULL
    step_ixs <- meta %>% pull(step_ix)
    
    # sample individual cells
    for (row_i in 1:num_cells) {
            sampled_ix <- sample(step_ixs, size=1)
            row <- meta %>% filter(step_ix == sampled_ix)
            
            if (row_i == 1) {
                sampled_cells <- row
            } else {
                sampled_cells <- rbind(row, sampled_cells)
            }
    }
    
    sampled_cells <- sampled_cells %>% as.data.frame
    return (sampled_cells)
}

sampling_cells_per_grna <- function(model, meta, ncells_per_grna) {
    grna_sampled_cells <- list()
    sims <- model$sim_meta$sim_name

    for (sim in sims) {
        grna_meta <- meta %>% filter(grepl(sim, sim_name))
        
        if (model$simulation_system$crispr_param$crispr_type == 'Knockout') {
            perc <- sim_meta %>% filter(sim_name == sim) %>% pull(sample_percentage)
            ncells <- ceiling(ncells_per_grna * perc)
        } else {
            ncells <- ncells_per_grna
        }
        
        grna_sampled_cells[[sim]] <- sampling_cells(grna_meta, ncells)
    }
    
    grna_sampled_cells <- grna_sampled_cells %>% bind_rows()
    return (grna_sampled_cells)
}