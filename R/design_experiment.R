#' Set Parameters for In Silico gRNA Experiment
#' @param ncellpops null
#' @param edge_probs null
#' @param ntfs_per_cellpop null
#' @param nhk null
#' @param ngenes null
#' @param nsimulations null
#' @param ngrnas_per_target null
#' @param ncells_per_model null
#' @param census_interval null
#' @param tau null
#' @param ntarget null
#' @param target_tf_only null
#' @param grna_library null
#' @param experiment_type null
#' @param ctrl_label null
#' @param cache_dir null
#' @export
design_experiment <- function(
    ncellpops=1,
    edge_probs=0.2,
    ntfs_per_cellpop=sample(10:20, ncellpops, replace=TRUE),
    grna_to_ctrl_ratio=runif(1),
    nhk=9000,
    negenes=500,
    nsimulations=100,
    ngrnas_per_target=3,
    ncells_per_model=1000,
    census_interval=10,
    tau=100/3600,
    ntargets=sum(ntfs_per_cellpop),
    ncells_in_experiment=10000,
    target_tf_only=TRUE,
    grna_library='default',
    experiment_type='ko',
    ctrl_label='CTRL',
    cache_dir=file.path(getwd(), 'cache')) {
    
    cellpop_params <- config_cellpop_params(ncellpops, edge_probs, ntfs_per_cellpop)
    
    simulator_params <- config_simulator_params(nhk, negenes, nsimulations,
                                                ntfs_per_cellpop, ncells_per_model,
                                                census_interval, tau)
    
    experiment_params <- config_experiment_params(ntargets, target_tf_only, grna_library,
                                                  ngrnas_per_target, ncells_in_experiment,
                                                  grna_to_ctrl_ratio, experiment_type, ctrl_label,
                                                  cache_dir)
    
    crigen_obj <- list(cellpop_params=cellpop_params,
                       simulator_params=simulator_params,
                       experiment_params=experiment_params)
    
    return (crigen_obj)
}

config_cellpop_params <- function(ncellpops, edge_prob, ntfs_per_cellpop) {
    
    if (edge_prob >= 1){
        stop("Edge Probability must be less than 1...") 
    }
    
     if (!as.vector(ntfs_per_cellpop)){
         stop("ntfs_per_cellpop parameter must be a vector...")
     }
    
    if (ncellpops != length(ntfs_per_cellpop)){
        stop("ncellpops doesn't equal ntfs_per_cellpop...")
    }
    
    cellpop_params <- list(ncellpops=ncellpops, edge_prob=edge_prob, ntfs_per_cellpop=ntfs_per_cellpop)
    return (cellpop_params)
}

config_simulator_params <- function(nhk, negenes, nsimulations, ntfs_per_cellpop,
                                    ncells_per_model, census_interval, tau) {
    sim_params <- list(num_hks = nhk,
                       num_targets = negenes,
                       num_tfs = sum(ntfs_per_cellpop),
                       num_cells = ncells_per_model,
                       num_simulations = nsimulations,
                       census_intervals = census_interval,
                       tau = tau)
    
    return (sim_params)
}

config_experiment_params <- function(ntargets, target_tf_only, grna_library,
                                     ngrnas_per_target, ncells_in_experiment,
                                     grna_to_ctrl_ratio, experiment_type, ctrl_label,
                                     cache_dir) {
    if (ntargets <= 0) {
        stop("Number of Target Genes Specified is Equal to or less that 0...")
    }
    
    if (!target_tf_only %in% c(TRUE, FALSE)) {
        stop("target_tf_only Parameter is not a Boolean value...")
    }
    
    if (!grna_library %in% c("perfect", "default")) {
        stop("gRNA Library Specified is Currently Not Supported by Crigen...")
    }
    
    if (!experiment_type %in% c("ko", "interference")) {
        stop("Type of CRISPR Experimental Specified is Currently Not Supported by Crigen...")
    }
    
    ncells_to_sample = (ncells_in_experiment / ntargets)
    expr_params <- list(ntargets=ntargets,
                        grna_library=grna_library,
                        target_tf_only=target_tf_only,
                        ngrnas_per_target=ngrnas_per_target,
                        ncells_in_experiment=ncells_in_experiment,
                        ncells_to_sample=ncells_to_sample,
                        experiment_type=experiment_type,
                        grna_to_ctrl_ratio=grna_to_ctrl_ratio,
                        ctrl_label=ctrl_label,
                        cache_dir=cache_dir)
    
    return (expr_params)
}