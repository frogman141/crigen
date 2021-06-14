#' @export 
sample_simulated_cells <- function(crigen_obj) {
    # merge all of the simulated sces and downsample cells
    new_sces = list()
    dyngen_meta <- crigen_obj$experiment_meta$dyngen_meta
    
    for (row_i in 1:nrow(dyngen_meta)){
        row = dyngen_meta[row_i, ]
        sce = readRDS(row$sce_fp)
        
        sampled_cells = downsample_cells(colnames(sce), row$sample_percentage)
        new_sces[[row$model_name]] = sce[, sampled_cells]
    }
    
    sce <- merge_sces(new_sces)
    return (sce)
}

#' @export 
merge_sces <- function(sces) {
    metadata <- list(column=NULL, row=TRUE)
    assay_list <- list(counts=NULL, logcounts=NULL, counts_spliced=NULL, counts_unspliced=NULL, counts_protein=NULL)
    
    for (model_sce in sces) {
        metadata <- bind_metadata(model_sce, metadata)
        assay_list <- bind_count_data(model_sce, assay_list)
    }

    sce <- SingleCellExperiment(assay=assay_list, colData=metadata$column, rowData=metadata$row)
    return (sce)
}

bind_count_data <- function(current_sce, assays) {
    # merging all of the matrices by column and column metadata by row. 
    # because i am merging dyngen objects with the same kinetics row data is the same
    assays$counts <- cbind(assays$counts, assay(current_sce, 'logcounts'))
    assays$logcounts <- cbind(assays$logcounts, assay(current_sce, 'logcounts'))
    assays$counts_spliced <- cbind(assays$counts_spliced, assay(current_sce, 'logcounts'))
    assays$counts_unspliced <- cbind(assays$counts_unspliced, assay(current_sce, 'logcounts'))
    assays$counts_protein <- cbind(assays$counts_protein, assay(current_sce, 'logcounts'))
    
    return (assays)
}

bind_metadata <- function(current_sce, metadata){
    metadata$column <- rbind(metadata$column, colData(current_sce))

    if (isTRUE(metadata$row)) {
        metadata$row <- rowData(current_sce)
    }
    
    return (metadata)
}

# downsample_grnas <- function(sce, grna_to_ctrl_ratio){
#     cells_to_keep = c()
#     cell_meta <- colData(sce) %>% as.data.frame
#     grnas <- cell_meta$grna_name
    
#     for (grna in grnas) {
#         cells <- cell_meta %>% filter(grna_name == grna) %>% rownames()
        
#         if (grna == 'CTRL'){
#             cells_to_keep <- c(cells_to_keep, cells)
#             next()
#         }
        
#         sampled_cells <- downsample_cells(cells, grna_to_ctrl_ratio)
#         cells_to_keep <- c(cells_to_keep, sampled_cells)
#     }
    
#     sce <- sce[, cells_to_keep]
#     return (sce)
# }

#' @export 
downsample_cells <- function(cells, sample_percentage) {
    # random sample the cells we want to use 
    ncells <- round(length(cells) * sample_percentage)
    sampled_cells <- sample(cells, ncells)
        
    return (sampled_cells)
}
