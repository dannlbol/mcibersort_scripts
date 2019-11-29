# Support functions for MethylCIBERSORT

Custom scripts used within "Pediatric Pan-CNS Tumor Analysis of Immune-cell Infiltration Identifies Correlates of Antitumor Immunity" by Grabovska et al. (2019) 

- [Overview](#overview)
- [Documentation](#documentation)
- [feature.select.new)(#feature.select.new)
- [make.synth.mix)(#make.synth.mix)
- [System Requirements](#system-requirements)
- [Hardware requirements](#hardware-requirements)
- [Software Requirements](#software-requirements)
- [Usage Guide](#usage-guide)
- [Issues](https://github.com/dannlbol/mcibersort_scripts/issues)

# Overview
The functions provided serve two main purposes:

* `feature.select.new()`: facilitate and extend the 'feature selection' functionality of the `MethylCIBERSORT` R software package
* `make.synth.mix()`: facilitate benchmarking of resulting signature matrix by generating *synthetic mixtures* of test populations with known input proportions 

# Documentation
## feature.select.new

## make.synth.mix

# System Requirements
## Hardware Requirements
Functions provided here are compatible with any standard computer with enough RAM to support the in-memory operations.
The functions are supported on any operating system which supports `MethylCIBERSORT 0.2.1`

## Software Requirements
Functions provided import and depend on a number of R packages. Functionality has been tested on *R 3.5.3* *Ubuntu 16.04.5* LTS*
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

Feature selection example:
```
source("./feature_select_WORKING_2019.R")
require(FlowSorted.Blood.450k)
idx <- which(FlowSorted.Blood.450k$CellType %in% c("Bcell","CD4T","CD8T","NK"))
bvals <- getBeta(FlowSorted.Blood.450k[, idx])
getwd() ## outs stored in working directory
feature.select.new(Stroma.matrix = bvals,
                   Phenotype.stroma = as.factor(FlowSorted.Blood.450k$CellType[idx]),
                   sigName = "test")
```

Generating *synthetic mixtures* example:
```
source("./synth_mix_2019.R")
require(FlowSorted.Blood.450k)
idx <- which(FlowSorted.Blood.450k$CellType %in% c("Bcell","CD4T","CD8T","NK"))
bvals <- getBeta(FlowSorted.Blood.450k[, idx])
b.sig <- read.delim("OPT2_0.2_200_SigEdit.txt", row.names = 1, header = TRUE)
bvals <- bvals[rownames(b.sig), ]
all(rownames(bvals)==rownames(b.sig))

make.synth.mix(input.data = bvals,
               pop.factor = as.factor(FlowSorted.Blood.450k$CellType[idx]),
               pop.rows = 10,
               n.cores = 1)
```
