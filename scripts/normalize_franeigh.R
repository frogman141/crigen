#!/home/baker02/miniconda3/envs/scRNA_env/bin/Rscript

# execute the analysis pipeline for the ER data
# each step will write intermediate output files to disk

#! How many cores per task?
#SBATCH --cpus-per-task=4
#! How much memory do you need?
#SBATCH --mem=16G
#! How much wallclock time will be required?
#SBATCH --time=10-00:00:00
#! Specify your email address here otherwise you won't recieve emails!
#SBATCH --mail-user=alexander.baker@cruk.cam.ac.uk
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
#SBATCH --no-requeue
#! General partition
#SBATCH -p general

library(scran)
library(scater)
library(DESeq2)
library(argparse)
library(tidyverse)
library(BiocParallel)

multicore_param = MulticoreParam(4)
register(multicore_param, default=TRUE)

# Psuedo Bulk Differential Expression Utility Functions 
conduct_dea_with_deseq2 = function(dds, cache_fp){
  # Remove genes with low counts
  keep = rowSums(counts(dds)) >= 10
  dds = dds[keep,]
    
  # estimate size factors and train deseq2 model
  dds = DESeq(dds, parallel = TRUE)
  saveRDS(dds, cache_fp)

  return (dds)
}

fetch_dea_results = function(dds, dea_res_fp, gene_col='gene', control_label='Control'){
  # fetching the results for all shRNA and applying shrinkage to there log2 foldchange estimates.
  
  result = list()
  
  # fetching dea results for all target genes
  for(gene in unique(dds[[gene_col]])){
    print (gene)
    # skip control
    if (gene == control_label){
      next()
    }
    
    gene_contrast = c(gene_col, gene, control_label)
    gene_coefficient = paste0(gene_col, "_", gene, "_vs_", control_label)
      
    print(paste0("Fetching the Results for ", gene))
    res_unshrunken = results(dds, alpha = 0.05, contrast = gene_contrast, parallel = TRUE)
    result[[gene]] = lfcShrink(dds, gene_coefficient, res=res_unshrunken, type="normal", parallel = TRUE)
  }
  
  # save and returned fetched dea results
  saveRDS(result, dea_res_fp)
  return (result)
}

find_deg_of_target_genes = function(dea_results, col_of_interest='padj', threshold=0.05) {
  # looping through all target genes and extracting statistically significant differentially expressed genes 
  # with a median transcript count of at least 5. Return the union of all gene names. 
  deg_gene_list = c()
  
  for (target_gene in names(dea_results)){
    res = dea_results[[target_gene]]
    gene_degs = res[!is.na(res[[col_of_interest]]) & res[[col_of_interest]] <= threshold, ]
    deg_gene_list = union(deg_gene_list, rownames(gene_degs))
  }
  
  return(deg_gene_list)
}

extract_de_results = function(deg_results, list_of_degs, gene_metadata, padj_col='padj', logfc_col='log2FoldChange'){
  # Extracting Differential Expression Results of Target Genes
  df = data.frame()
  deg_target_genes = names(deg_results)
  
  for(target_gene in rownames(gene_metadata)) {
    
    if (!target_gene %in% deg_target_genes){
      next()
    }
    
    res = deg_results[[target_gene]][list_of_degs, ]
    res[[padj_col]][is.na(res[[padj_col]])] = 1
    
    target_df = data.frame(
      logfc = res[[logfc_col]],
      log_padj = -log(res[[padj_col]], 10),
      target_gene = rep(target_gene, length(degs)),
      deg = degs
    )
    
    rownames(target_df) = NULL
    df = rbind(df, target_df)
  }
  
  df$logfc[is.na(df$logfc)] = 0 
  
  return (df)
}


ctrl_label = 'Control'
franiegh_dir = '/scratcha/fmlab/baker02/Edwards/data/Franiegh'
############################### Loading the Franeigh et al 2021 dataset ###############################
parser = ArgumentParser(description='Conducting Differential Expression Analysis on Sub Groups of Franiegh Data')
parser$add_argument('group', help='sub group of Franiegh Data', type='character')

arguments =  parser$parse_args()
franiegh_group = arguments$group

################ Conducting Differential Expression Analysis on Franeigh et al. 2021 ################

franiegh_group = arguments$group
raw_deseq2_fp = paste0(franiegh_dir, '/DEGs/DESeq2/', franiegh_group, '.Rds')

if (!file.exists(raw_deseq2_fp)){
#     # because the data provided was preprocessed i just going to round
#     summed = aggregateAcrossCells(sce, id=colData(sce)[,c("gene", "random_group")])
#     meta = colData(summed)

#     gene_list = as.character(unique(meta$gene))
#     meta$random_group_label = factor(meta$random_group, levels=unique(meta$random_group))
#     meta$gene_label = factor(meta$gene, levels=c(ctrl_label, gene_list[-which(gene_list == ctrl_label)]))

#     dds = DESeqDataSetFromMatrix(countData = assay(summed, 'counts'), colData = meta, design = ~ random_group_label + gene_label)
#     saveRDS(dds, raw_deseq2_fp)
    dds = readRDS(paste0(franiegh_dir, '/DEGs/Prep_DDS/', franiegh_group, '.Rds'))
    dds = conduct_dea_with_deseq2(dds, raw_deseq2_fp)
    
} else {
    print ("Loading Raw Deseq2 Object...")
    dds = readRDS(raw_deseq2_fp)
}

########## Extracting the Results from Differential Expression Analysis ##########

dea_results_fp = paste0(franiegh_dir, '/DEGs/results/', franiegh_group, '.Rds')

if (!file.exists(dea_results_fp)){
    print ("Extracting Differential Expression Analysis Results from DESeq2 Object...")
    results = fetch_dea_results(dds, dea_results_fp, 'gene_label', ctrl_label)
} else {
    print ("Loading Differential Expression Analysis Results...")
    results = readRDS(dea_results_fp)
}
