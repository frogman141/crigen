# CRIGEN

CRIGEN is a Single Cell CRISPR Simulator.

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
