#' @export
initialize_kinetics <- function(model) {
    
    assert_that(!is.null(model$feature_info),
                !is.null(model$feature_network))
    
    # generate kinetics params and formulae
    message("Initializing Networks Kinetics...")
    model <- initializing_gene_kinetics(model)
    reactions <- generate_reactions(model)
  
    # create variables
    fid <- model$feature_info$feature_id
    model$feature_info$mol_premrna <- paste0("mol_premrna_", fid)
    model$feature_info$mol_mrna <- paste0("mol_mrna_", fid)
    model$feature_info$mol_protein <- paste0("mol_protein_", fid)
  
    molecule_ids <- c(model$feature_info$mol_premrna, 
                      model$feature_info$mol_mrna,
                      model$feature_info$mol_protein)
    
    initial_state <- set_names(rep(0, length(molecule_ids)), molecule_ids)
    
    # return system
    model$simulation_system$reactions <- reactions
    model$simulation_system$molecule_ids <- molecule_ids
    model$simulation_system$initial_state <- initial_state
    
    return (model)
}

#' @export
inject_sim_noise <- function(sim, model) {
    # get noise parameters
    sd <- model$simulation_system$noise_sd
    mean <- model$simulation_system$noise_mean
    
    # inject noise into both individual gene rates and their interactions
    info <- model$feature_info %>% rate_noise(mean, sd)
    network <- model$feature_network %>% interaction_noise(mean, sd)
    
    out <- calculate_dissociation(info, network)
    params <- extract_parameters(out$feature_info, out$feature_network)
    
    values <- params %>% select(value) %>% t %>% as.data.frame
    colnames(values) <- params$id
    rownames(values) <- sim
    
    return (values)
}

#' @export
calculate_dissociation <- function(feature_info, feature_network) {
  remove <- c("max_premrna", "max_mrna", "max_protein", "dissociation", "k", "max_protein")
  
  feature_info <- feature_info[, !colnames(feature_info) %in% remove]
  feature_network <- feature_network[, !colnames(feature_network) %in% remove]
  
  feature_info <- feature_info %>%
                    mutate(max_premrna = transcription_rate / (mrna_decay_rate + splicing_rate),
                           max_mrna = splicing_rate / mrna_decay_rate * max_premrna,
                           max_protein = translation_rate / protein_decay_rate * max_mrna)
  
  feature_network <- feature_network %>% 
                        left_join(feature_info %>% select(from = feature_id, max_protein), by = "from") %>% 
                        mutate(dissociation = max_protein / 2)
  
  return(lst(feature_info, feature_network))
}

#' @export
extract_parameters <- function(feature_info, feature_network) {
    # extract production / degradation rates, ind and bas
    feature_params <- feature_info %>% 
                        select(feature_id, transcription_rate, splicing_rate, translation_rate,
                               mrna_decay_rate, protein_decay_rate, bas = basal, ind = independence) %>% 
                        gather("param", "value", -feature_id) %>% 
                        mutate(id = paste0(param, "_", feature_id), type = "feature_info")

    # extract dis, hill, str
    edge_params <- feature_network %>% 
                        select(from, to, dis = dissociation, hill, str = strength) %>% 
                        gather("param", "value", -from, -to) %>% 
                        mutate(id = paste0(param, "_", from, "_", to), type = "feature_network")

    parameters <- bind_rows(feature_params, edge_params) %>% select(id, value)
    return (parameters)
}

#' @export
rate_noise <- function(feature_info, mean, std) {
    mrna_halflife <- protein_halflife <- NULL
    rate_cols <- c("transcription_rate", "splicing_rate", "translation_rate", "mrna_halflife", "protein_halflife")
    
    feature_info <- feature_info %>% 
                          mutate_at(c("basal"), ~ pmin(1, . * rnorm(length(.), mean = mean, sd = std))) %>% 
                          mutate_at(rate_cols, ~ . * rnorm(length(.), mean = mean, sd = std)) %>% 
                          mutate(mrna_decay_rate = log(2) / mrna_halflife, 
                                 protein_decay_rate = log(2) / protein_halflife)
    
    return (feature_info)
}    

#' @export
interaction_noise <- function(feature_network, mean, std) {
    # satisfy r cmd check
    feature_network <- feature_network %>%
                            mutate_at(c("strength", "hill"), ~ . * rnorm(length(.), mean = mean, sd = std))
    
    return(feature_network)
}

initializing_gene_kinetics <- function(model) {
  
    # fetch feature info and network
    rate_columns <- c("transcription_rate", "splicing_rate", "translation_rate", "mrna_halflife", "protein_halflife", "independence")
    interaction_columns <- c("effect", "strength", "hill")
    
    feature_info <- model$feature_info %>% add_columns(rate_columns, NA_real_)
    feature_network <- model$feature_network %>% add_columns(interaction_columns, NA_real_)
    
    # generate relatively stable kinetics for TFs
    feature_info <- feature_info %>% sample_rate_kinetics()
    feature_network <- feature_network %>% sample_interaction_kinetics()
    
    # calculate k
    dis_out <- calculate_dissociation(feature_info, feature_network)
    feature_info <- dis_out$feature_info
    feature_network <- dis_out$feature_network
    
    # calculate ba and a
    basal_df <- feature_network %>% 
                      rename(feature_id = .data$to) %>% 
                      group_by(feature_id) %>% 
                      summarise(basal = calculate_basal(.data$effect))
    
    feature_info <- feature_info %>%
                        left_join(basal_df, by = "feature_id") %>%
                        mutate(basal=ifelse(is.na(basal), 1, basal))
    
    feature_network <- qc_tf_interactions(feature_network, feature_info)

    model$feature_info <- feature_info
    model$feature_network <- feature_network

    return (model)
}

#' @importFrom GillespieSSA2 reaction
generate_reactions <- function(model) {
      
    # add helper information to feature info
    feature_info <- model$feature_info %>% 
                        left_join(model$feature_network %>% 
                                        group_by(feature_id = to) %>% 
                                        summarise(regulators = list(data.frame(from = .data$from,
                                                                               effect = .data$effect,
                                                                               strength = .data$strength)),
                                                  .groups = "drop"),
                                  by = "feature_id") %>% 
                        left_join(model$feature_network %>% 
                                    group_by(feature_id = from) %>% 
                                    summarise(num_targets = n()),
                                  by = "feature_id")

    # generate formula per feature
    out <- pbapply::pblapply(
        seq_len(nrow(feature_info)),
        function(i) {
            info <- feature_info %>% extract_row_to_list(i)
            fid <- info$feature_id

            w <- paste0("mol_premrna_", fid)
            x <- paste0("mol_mrna_", fid)
            y <- paste0("mol_protein_", fid)

            transcription_rate <- paste0("transcription_rate_", fid)
            splicing_rate <- paste0("splicing_rate_", fid)
            translation_rate <- paste0("translation_rate_", fid)
            mrna_decay_rate <- paste0("mrna_decay_rate_", fid)
            protein_decay_rate <- paste0("protein_decay_rate_", fid)

            basal <- paste0("bas_", fid)
            independence <- paste0("ind_", fid)

            if (!is.null(info$regulators)) {
                rid <- info$regulators$from
                eff <- info$regulators$effect
                str <- info$regulators$strength
                reg_ys <- paste0("mol_protein_", rid)
                reg_diss <- paste0("dis_", rid, "_", fid)
                reg_hills <- paste0("hill_", rid, "_", fid)
                reg_strs <- paste0("str_", rid, "_", fid)
                regulation_var <- paste0("chi_", rid, "_", fid)

                reg_affinity_calc <- paste(paste0(regulation_var, " = ", reg_strs, " * pow(", reg_ys, "/", reg_diss, ", ", reg_hills, "); "), collapse = "")
                
                numerator <-
                  if (sum(eff > 0) > 0) {
                    paste0(basal, " - pow(", independence, ",", sum(eff > 0), ") + ", paste("(", regulation_var[eff > 0], " + ", independence, ")", collapse = " * ", sep = ""))
                  } else {
                    basal
                  }
                denominator <- paste("(", regulation_var, " + 1)", collapse = " * ", sep = "")

                act_function <- paste0(reg_affinity_calc, transcription_rate, " * (", numerator, ")/(", denominator, ")")
             } else {
                act_function <- paste0(transcription_rate, " * ", basal)
                regulation_var <- character()
             }

             formulae <- list(
                # pre-mRNA production
                reaction(
                  name = paste0("transcription_", fid),
                  effect = set_names(1, w),
                  propensity = paste0(act_function)
                ),
                # splicing
                reaction(
                  name = paste0("splicing_", fid),
                  effect = set_names(c(1, -1), c(x, w)),
                  propensity = paste0(splicing_rate, " * ", w)
                ),
                # protein production
                reaction(
                  name = paste0("translation_", fid), 
                  effect = set_names(1, y),
                  propensity = paste0(translation_rate, " * ", x)
                ),
                # pre-mRNA degradation
                reaction(
                  name = paste0("premrna_degradation_", fid),
                  effect = set_names(-1, w),
                  propensity = paste0(mrna_decay_rate, " * ", w)
                ),
                # mRNA degradation
                reaction(
                  name = paste0("mrna_degradation_", fid),
                  effect = set_names(-1, x),
                  propensity = paste0(mrna_decay_rate, " * ", x)
                ),
                # protein degradation
                reaction(
                  name = paste0("protein_degradation_", fid),
                  effect = set_names(-1, y),
                  propensity = paste0(protein_decay_rate, " * ", y)
                )
              )

          formulae[[1]]$buffer_ids <- regulation_var
          return (formulae)
    })
  
    return(unlist(out, recursive = FALSE))
}

calculate_basal <- function(effects) {
    basal <- case_when(all(effects == -1) ~ 1,
                       all(effects == 1) ~ 0.0001,
                       TRUE ~ 0.5)
    
    return (basal)
}

sample_rate_kinetics <- function(feature_info) {
     feature_info <- feature_info %>%
                        mutate(independence = 1,
                               transcription_rate = runif(n(), 10, 20),
                               translation_rate = runif(n(), 100, 150),
                               mrna_halflife = runif(n(), 2.5, 5),
                               protein_halflife = runif(n(), 5, 10),
                               splicing_rate = log(2) / 2,
                               mrna_decay_rate = log(2) / .data$mrna_halflife,
                               protein_decay_rate = log(2) / .data$protein_halflife)
    
    return (feature_info)
}

#' @importFrom dyngen rnorm_bounded
sample_interaction_kinetics <- function(feature_network) {
    # this function needs to worked on so there is at least one condition where a TF is expressed
    feature_network <- feature_network %>%
                            mutate(effect = sample(c(-1L, 1L), n(), replace=TRUE, prob=c(.25, .75)),
                                   strength = 10 ^ runif(n(), log10(1), log10(100)),
                                   hill = rnorm_bounded(n(), 2, 2, min = 1, max = 10))
    
    return (feature_network)
}

qc_tf_interactions <- function(feature_network, feature_info) {
    tfs <- feature_info %>% filter(is_tf) %>% pull(feature_id)
    feature_network <- feature_network %>% mutate(row_i=row_number())
    
    for (tf in tfs) {
        tf_edges <- feature_network %>% filter(to == tf)
        all_edges <- tf_edges %>% nrow
        negative_edges <- tf_edges %>% filter(effect == -1L) %>% nrow
        
        if (all_edges == 0) {
            next()
        }
        
        if (negative_edges == all_edges) {
            edge_to_modify <- tf_edges %>% pull(row_i) %>% sample(1)
            feature_network[edge_to_modify, ]$effect <- 1
        }
    }
    
    return (feature_network)
}

simple_noise_simple <- function(feature_info, feature_network, mean, sd) {
    # satisfy r cmd check
    mrna_halflife <- protein_halflife <- NULL
    rate_cols <- c("transcription_rate", "splicing_rate", "translation_rate", "mrna_halflife", "protein_halflife")
    
    feature_info <- feature_info %>% 
                      mutate_at(c("basal"), ~ pmin(1, . * rnorm(length(.), mean = mean, sd = sd))) %>% 
                      mutate_at(rate_cols, ~ . * rnorm(length(.), mean = mean, sd = sd)) %>% 
                      mutate(mrna_decay_rate = log(2) / mrna_halflife, 
                             protein_decay_rate = log(2) / protein_halflife)
    
    feature_network <- feature_network %>%
                            mutate_at(c("strength", "hill"), ~ . * rnorm(length(.), mean = mean, sd = sd))
    
    return(lst(feature_info, feature_network))
}

add_columns <- function(df, colnames, fill = NA) {
    
    for (colname in colnames) {
        if (!colname %in% colnames(df)) {
            df[[colname]] <- fill
        }
    }
  
    return (df)
}