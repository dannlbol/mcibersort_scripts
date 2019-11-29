# Support functions for MethylCIBERSORT

Custom scripts used within "Pediatric Pan-CNS Tumor Analysis of Immune-cell Infiltration Identifies Correlates of Antitumor Immunity" by Grabovska et al. (2019) 

- [Overview](#overview)
- [Documentation](#documentation)
- [feature.select.new](#feature.select.new)
- [make.synth.mix](#make.synth.mix)
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
This function mimics the functionality of `FeatureSelect.V4()` from `MethylCIBERSORT 0.2.1`
### Default function call:
`feature.select.new(MaxDMPs = 100, deltaBeta = 0.2, useM = FALSE, CellLines.matrix = NULL, export = TRUE, export.fit = TRUE, export.cpg = TRUE, sigName = "methylCibersort", Stroma.matrix = NULL, Phenotype.stroma = NULL, FDR = 0.01, silent = TRUE)`

### Additional parmeters:
**useM**: specify whether conversion to M-values should be done before carrying out feature seletion
**export.fit**: whether to export the `limma` fit object during the feature selection
**export.cpg**: whether to export a table of the CpG probes selected alongside the name the name of the population which they were selected agaist

### Typical outputs:
* .txt file of the resulting signature matrix
* .rds file of the limma fit
* .txt file of CpGs selected and the comparison from which they were chosen
* .txt file of the Illumina/`minfi` CpG annotations for the relevant CpGs

## make.synth.mix
This function creates a matrix of proportions for a given set of data populations and generates weighted means based on those proportions for a given matrix of values, typically methylation array beta-values.
### Default function call:
`make.synth.mix(input.data = NULL, pop.factor = NULL, pop.rows = 100, output.dir = getwd(), output.name = gsub("-", "_", Sys.Date()), n.cores = 1)`

### Parameters:
**input.data**: a matrix of data, typically beta-values; rows = features (probes), columns = observations (samples)
**pop.factor**: a factor with levels describing the populations in input.data columns
**pop.rows**: numeric specifying how many rows of proportions to generate for each population
**output.dir**: location to save the resulting files
**output.name**: name to append to the resulting filenames
**n.cores**: set to >1 to run on multiple cores in parallel using `parallel::mcapply()`

### Typical outputs:
* .txt file of the resulting porportions matrix
* .txt file of the resulting mixtures

### Note on functionality:
The function populates the proportion table by column, with each population being assigned a set sequence of probabilities. A typical example for 3 populations, 6 rows per population:

<img src="/synth_mix_example.png" width="480">

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
parallel
```

# Usage Guide:

Feature selection example:
```
source("./feature_select_WORKING_2019.R")
require(FlowSorted.Blood.450k)
idx <- which(FlowSorted.Blood.450k$CellType %in% c("Bcell","CD4T","CD8T","NK"))
bvals <- getBeta(FlowSorted.Blood.450k[, idx])
getwd() ## outs stored in working directory
a <- Sys.time()
feature.select.new(Stroma.matrix = bvals,
                   Phenotype.stroma = as.factor(FlowSorted.Blood.450k$CellType[idx]),
                   sigName = "test")
b <- Sys.time()
b-a
# Time difference of 1.55212 mins
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
a <- Sys.time()
make.synth.mix(input.data = bvals,
               pop.factor = as.factor(FlowSorted.Blood.450k$CellType[idx]),
               pop.rows = 10,
               n.cores = 1)
b <- Sys.time()
b-a
## Time difference of 1.76159 secs
```
