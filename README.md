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

Generate 10000 cells from 1 cell population with 2 transcription factors (ntf), 10 housekeeping genes (nhk)

Number of egenes

```
library(crigen)
crg <- design_experiment(,experiment_type='ko',ntfs_per_cellpop=2, nhk=10, negenes=2) %>%
  create_experiment() %>%
  simulate_experiment()
```

### Simulate off target effects

Offtarget effects can be added to other experiments or can be used alone.


## Design experiment parameters and defaults

```
    ncellpops=1,
    edge_probs=0.2,
    ntfs_per_cellpop=sample(10:20, ncellpops, replace=TRUE), # number of transcription factors
    grna_to_ctrl_ratio=runif(1),
    nhk=9000,                         # number of house keeping genes
    negenes=500,                      # number of egenes (downstream reporter genes)
    nsimulations=100,
    ngrnas_per_target=3,
    ncells_per_model=1000,
    census_interval=10,
    tau=100/3600,                     #
    ntargets=sum(ntfs_per_cellpop),
    ncells_in_experiment=10000,
    target_tf_only=TRUE,
    grna_library='default',
    experiment_type='ko',             # ('ko', )
    ctrl_label='CTRL',
```

