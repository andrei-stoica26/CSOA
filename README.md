# CSOA
Cell Set Overlap Analysis (CSOA) is a tool for calculating per-cell gene 
signature scores in an scRNA-seq dataset. CSOA constructs a set for 
each gene in the signature, consisting of the cells that highly express the 
gene. Next, all overlaps of pairs of cell sets are computed, ranked, 
filtered and scored. The CSOA per-cell score is calculated by summing up all 
products of the overlap scores and the min-max-normalized expression of the 
two involved genes.

## Installation

To install CSOA, run the following R code:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("CSOA")
```

## Usage

The basic command to run CSOA is:

```
runCSOA(scObj, geneSets)
```

`scObj` must contain normalized and log-transformed gene expression data 
provided in one of the following formats:

- `Seurat`.
    - CSOA will use the expression matrix stored in the `data` layer.
- `SingleCellExperiment`
    - CSOA will use the expression matrix stored in the `logcounts` assay.
- `matrix`
- `dgCMatrix`

`geneSets` must be a named list of character vectors. 
