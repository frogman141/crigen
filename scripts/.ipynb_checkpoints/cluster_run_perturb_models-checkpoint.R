#!/home/baker02/miniconda3/envs/tidyscreen/bin/Rscript

#! How many cores per task?
#SBATCH --cpus-per-task=1
#! How much memory do you need?
#SBATCH --mem=4G
#! How much wallclock time will be required?
#SBATCH --time=1-00:00:00
#! Specify your email address here otherwise you won't recieve emails!
#SBATCH --mail-user=alexander.baker@cruk.cam.ac.uk
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
#SBATCH --no-requeue
#! General partition
#SBATCH -p general

library(dyno)
library(dyngen)
library(argparse)

############################ Parsing Command Line Argument ############################
parser = ArgumentParser(description='Conducting Differential Expression Analysis on Sub Groups of Franiegh Data')
parser$add_argument('--model-fp', help='File Path to Dyngen Model', type='character')

arguments = parser$parse_args()
model_fp = arguments$model_fp

################## create knockdown simulation and generate cells #####################
perturb_model = readRDS(model_fp)
perturb_model =  perturb_model %>% generate_cells()

saveRDS(perturb_model, perturb_model_fp)