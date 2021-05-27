#' @export
on_target <- function(perturb_genes, num_of_genes, library='default', experiment_type='KO'){
    
    on_target_activity <- sample_on_target_activity(num_of_genes, library)
    
    if (experiment_type == 'Activation') {
        on_target_activity = on_target_activity + 1
    }
    
    return (on_target_activity)
}

sample_on_target_activity <- function(num_of_genes, library) {
    
    # generate randomized on target if specified
    if (library == 'default'){
        on_target <- runif(num_of_genes)
    } 
    
    # converting on_target score into a KD multiplier Dyngen acccepts
    on_target = 1 - on_target
    return (on_target)
}

####### Going to work on this in a different file #######

all_ko_combinations <- function(perturbed_genes){
    # generating all of the potiential combinations of KO effects for a given gRNA
    combinations = list()

    for (gene in target_genes){
        combinations[[gene]] = 0:1
    }
    
    ko_combinations = combinations %>% as.data.frame %>% expand.grid()
    return (ko_combinations)
}

calc_ko_sample_percentages <- function(on_target, ko_combinations){
    # calculate the percentage of cells that should be sampled from 
    # a given KO combination model.
    probs_mat = data.frame()
    
    for (i_row in 1:nrow(ko_combinations)) {
        row = ko_combinations[i_row, ]
        
        for (i_col in 1:ncol(ko_combinations)) { 
            gene_on_target = on_target[i_col]
            value = ko_combinations[i_row, i_col]

            if (value == 1) {
                # probability of grna edit
                row[i_col] = gene_on_target
            } else {
                # probability of no edit
                row[i_col] = 1 - gene_on_target
            }
        }

        probs_mat = rbind(probs_mat, row)
    }

    sample_percentages = rowProds(as.matrix(probs_mat))
    return (sample_percentages)
}