# Support functions for MethylCIBERSORT

Custom scripts used within "Pediatric Pan-CNS Tumor Analysis of Immune-cell Infiltration Identifies Correlates of Antitumor Immunity" by Grabovska et al. (2019) 

- [Overview](#overview)
- [Documentation](#documentation)
- [System Requirements](#system-requirements)
- [Hardware requirements](#hardware-requirements)
- [Software Requirements](#software-requirements)
- [Usage Guide](#usage-guide)
- [Issues](https://github.com/dannlbol/mcibersort_scripts/issues)

# Overview
The functions provided serve two main purposes:

* facilitate and extend the 'feature selection' functionality of the `MethylCIBERSORT` R software package
* facilitate benchmarking of resulting signature matrix by generating *synthetic mixtures* of test populations with known input proportions 

# Documentation


# System Requirements
## Hardware Requirements
Functions provided here are compatible with any standard computer with enough RAM to support the in-memory operations.

## Software Requirements
Functions provided import and depend on a number of R packages. Functionality has been tested on *R 3.5.3*
### R Dependencies
`MethylCIBERSORT 0.2.1` was obtained by correspondence from Dr Ankur Chakravarthy.

```
magrittr
matrixStats
BiocGenerics
limma
dplyr
minfi
```

# Usage Guide:

```
source("./feature_select_WORKING_2019.R")
source("./synth_mix_2019.R")
```
