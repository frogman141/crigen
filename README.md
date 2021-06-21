# CRIGEN

CRIGEN is a Single Cell CRISPR Simulator.  It will generate SingleCellExperiment objects that simulate the effects of activating, inhibiting and knocking out specific genes, as well as CRISPR off target effects. Parameters controlling the underlying transcription factor networks are configurable.  

## Installation

This package installs from this repository with devtools, but first we need to install some dependencies manually.
```
 install.packages(devtools)
 devtools::install_github("dynverse/dyno")

 install.packages("BiocManager")
 BiocManager::install(version = "3.13")
 BiocManager::install(c("scater", "scran"))

 devtools::install_github("frogman141/crigen")
```

## Quick start guide


### Simulate a knock out experiment

Generate 10000 cells (ncells_in_experiment) from 1 cell population (ncellpops) with 1 of 5 transcription factors knocked out per cell (ntfs_per_cellpop), 500 downstream genes attached to each transcription factor, simulated from a network of 9000 housekeeping genes (nhk), and run 2 simulations (nsimulations).  

```
library(crigen)
crg <- design_experiment(,experiment_type='ko', ncells_in_experiment=10000, ncellpops=1, ntfs_per_cellpop=5, negenes=500, nhk=9000, nsimulations=2) %>%
  create_experiment() %>%
  simulate_experiment()
```

### Simulate off target effects

Offtarget effects can be added to other experiments or can be used alone. (TODO)


## Design experiment parameters and defaults

```
Parameters regarding generating cell populations
    ncellpops=1, # number of cell populations.  Each cell population will have a distinct network structure. 
    edge_probs=0.2, # probability that tfs are connected in initial network creation
    ntfs_per_cellpop=sample(10:20, ncellpops, replace=TRUE), # number of transcription factors that are knocked out per cell population, over the whole experiment
    grna_to_ctrl_ratio=runif(1), # number of cells transfected with guide rnas to control cells (generally 0.3)
    nhk=9000,                         # number of housekeeping genes in dyngen network.  This and negenes are used as a reference to simulate the transcription factor network from.
    negenes=500,                      # number of egenes attached to tf in dyngen network (note this is called target in sce object as in target of tf)
        
Network simulation parameters
    nsimulations=100,  #dyngen parameter, number of walks through graph; per model
    ngrnas_per_target=3, # number of guide rnas targeting one gene
    ncells_per_model=1000, # TODO -- is this calculated from the number of cells?
    census_interval=10, # dyngen parameter, time interval
    tau=100/3600,         # dyngen parameter




            #
    ntargets=sum(ntfs_per_cellpop), # number of targets
    ncells_in_experiment=10000, # number of cells
    target_tf_only=TRUE, # as opposed to being housekeeping or egenes
    grna_library='default', # uniform distribution for on target; poisson for off target **Alex**
    experiment_type='ko',             # ('ko', ‘interference’, ‘activation’)
    ctrl_label='CTRL',
```

The number of calls that crigen will make to dyngen is equal to the number of models * the number of simulations, where the number of models = Target genes * guides/target * max(1, # off targets).

## Crigen object

```
names(crg)
sce <- crg$sce
```

Parameters about knockdown effects are stored in crg2$experiment_meta$grna_meta.
Dyngen base model and run specific models are stored in object crg3$simulator_models.


