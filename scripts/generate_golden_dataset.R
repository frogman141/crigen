#!/home/baker02/miniconda3/envs/tidyscreen/bin/Rscript

#! How many cores per task?
#SBATCH --cpus-per-task=1
#! How much memory do you need?
#SBATCH --mem=64G
#! How much wallclock time will be required?
#SBATCH --time=4-00:00:00
#! Specify your email address here otherwise you won't recieve emails!
#SBATCH --mail-user=alexander.baker@cruk.cam.ac.uk
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
#SBATCH --no-requeue
#! General partition
#SBATCH -p general

library(dyno)
library(scran)
library(scater)
library(slurmR)
library(dyngen)
library(igraph)
library(argparse)
library(argparse)
library(parallel)
library(dynutils)
library(tidyverse)
library(assertthat)
library(data.table)
library(SingleCellExperiment)
source("/scratcha/fmlab/baker02/tidyscreen/scripts/perturbation.R")
source("/scratcha/fmlab/baker02/tidyscreen/scripts/backbone_cellpop.R")


############################ Parsing Command Line Argument ############################
parser = ArgumentParser(description='Conducting Differential Expression Analysis on Sub Groups of Franiegh Data')
parser$add_argument('dataset_count', help='Ticker that is counting the number of datasets being generated', type='character')
parser$add_argument('--dataset-dir', help='There directory where Dyngen will save simulated data', type='character')
parser$add_argument('--num-cells', help='Number of Cells for Dyngen to Simulate', type='integer', default=1000)
parser$add_argument('--egenes-count', help='Number of Effect Genes to simulate', type='integer', default=250)
parser$add_argument('--hk-count', help='Number of House Keeping Genes to simulate', type='integer', default=4750)
arguments = parser$parse_args()

dataset_dir = arguments$dataset_dir
dataset_count = arguments$dataset_count
data_dir = file.path('/scratcha/fmlab/baker02/tidyscreen/data/', dataset_dir)

num_cells = arguments$num_cells

hk_count = arguments$hk_count
egenes_count = arguments$egenes_count

# generate the cell populations gene regulatory module network. Will requested 15 target genes
cellpop_backbone = backbone_single_cell_population(min_mod=15, max_mod=15)

# initialize the model and start the simulation
config = initialise_model(
    backbone = cellpop_backbone,
    num_cells = num_cells,
    num_targets = hk_count,
    num_hks = egenes_count,
    gold_standard_params = gold_standard_default(
      census_interval = 1,
      tau = 100 / 3600),
    simulation_params = simulation_default(
      census_interval = 10,
      ssa_algorithm = ssa_etl(tau = 300 / 3600),
      experiment_params = simulation_type_wild_type(
        num_simulations = 100
    )),
    experiment_params = experiment_snapshot(
      realcount = "GSE100866_CBMC_8K_13AB_10X-RNA_umi"
    )
)

config$num_cores = 5

# creating transcription factor network and generate kinetics
model_common = config %>%
                  generate_tf_network() %>%
                  generate_feature_network() %>% 
                  generate_kinetics() %>%
                  generate_gold_standard()

network_tfs = model_common$feature_info %>% filter(is_tf == TRUE & burn == FALSE) %>% pull(feature_id) %>% as.character

start.time = Sys.time()

perturb_expr = run_perturbation_experiment(model_common, network_tfs, cluster=TRUE)
sce = as_sce(perturb_expr)

end.time = Sys.time()
taken.time = end.time - start.time

message(taken.time)

dyngen_filename = paste0('dyngen_model_', dataset_count, '.Rds')
sce_filename = paste0('sce_of_dyngen_model_', dataset_count, '.Rds')

sce_fp = file.path(data_dir, sce_filename)
dyngen_fp = file.path(data_dir, dyngen_filename)

saveRDS(sce, sce_fp)
saveRDS(perturb_expr, dyngen_fp)
