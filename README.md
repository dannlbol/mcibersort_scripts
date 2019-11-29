# Support functions for MethylCIBERSORT

Custom scripts used within "Pediatric Pan-CNS Tumor Analysis of Immune-cell Infiltration Identifies Correlates of Antitumor Immunity" by Grabovska et al. (2019) 

- [Overview](#overview)
- [Documentation](#documentation)
- [System Requirements](#system-requirements)
- [Hardware requirements](#hardware-requirements)
- [Software Requirements](#software-requirements)
- [Issues](https://github.com/dannlbol/mcibersort_scripts/issues)

# Overview
Text

# Documentation
Text

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
