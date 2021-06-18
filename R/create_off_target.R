#' @export
off_target <- function(off_target_genes, library='default'){
    num_off_targets <- sample_off_target_count(library)
    off_target_genes <- select_off_target_genes(off_target_genes, num_off_targets)
    num_of_genes <- num_off_targets + 1
    
    return (list(num_of_genes=num_of_genes, num_off_targets=num_off_targets, off_target_genes=off_target_genes))
}

sample_off_target_count = function(library) {
    # Select the number off target genes through random sampling of 
    # off-target probability mass functions. The default option is to sample
    # a Poisson Distribution with lambda of 1.
    num_off_targets <- 0
    grna_lib_meta <- readRDS('../metadata/grna_libraries_meta.Rds')   
    
    # sampling from the library distribution. Poisson with lambda=1 is the 
    # default distribution to sample from.
    if (library == 'perfect') {
        num_off_targets <- 0
    } else if (library == 'default'){
        num_off_targets <- rpois(1, lambda=1)
    } else if (library %in% colnames(grna_lib_meta)){
        num_off_targets <- sample_off_target_dist(grna_lib_meta[[library]]$pdf_off)
    } else {
        message("Specified Library is currently not covered by crigen...")
        exit()
    }
    
    return (num_off_targets)
}

# add code to determine which distribution to use for sampling off-target genes
select_off_target_genes <- function (genes, num_off_target) {
    
    if (num_off_target > length(genes)){
        message("Number of Off Targets is greater than the number of avaliable to sample. Sampling all genes..")
        num_off_target <- length(genes)
    }
    
    # sample from gene list provided the off target gene for the gRNA
    off_target_genes <- sample(genes, num_off_target)
    
    return (off_target_genes)
}

sample_off_target_dist <- function(pdf) {
    off_target <- sample(pdf$mismatch_up_to_2, size = 1, prob = pdf$prob)
    return (off_target)
}
