#!/home/baker02/miniconda3/envs/perturbverse/bin/Rscript
#! How many cores per task?
#SBATCH --cpus-per-task=1
#! How much memory do you need?
#SBATCH --mem=16G
#! How much wallclock time will be required?
#SBATCH --time=1-00:00:00
#! Specify your email address here otherwise you won't recieve emails!
#SBATCH --mail-user=alexander.baker@cruk.cam.ac.uk
#! Uncomment this to prevent the job from being requeued (e.g. if
#! interrupted by node failure or system downtime):
#SBATCH --no-requeue
#! General partition
#SBATCH -p general

library(crigen)

test_obj_fp = 'test_crigen.Rds'

if(!file.exists(test_obj_fp)) {
    crg <- design_experiment(, grna_library='brunello', grna_to_ctrl_ratio=0.3, nhk=500, negenes=500) %>% create_experiment()
    saveRDS(crg, test_obj_fp)
} else {
    crg <- readRDS(test_obj_fp)
}

crg <- crg %>% simulate_experiment(, slurm=TRUE)
saveRDS(crg, test_obj_fp)
