# CSOA
Cell Set Overlap Analysis (CSOA) is a tool for calculating per-cell gene 
signature scores in a scRNA-seq dataset. 

For each signature gene, CSOA constructs a set consisting of cells 
highly expressing the gene. Next, all overlaps of pairs of cell sets are 
computed, ranked, filtered and scored. The CSOA per-cell score is calculated by 
summing up, over all top cell set overlaps, the product of the overlap score 
and the min-max-normalized expression of each of the two genes. 

## Installation

To install CSOA, run the following R code:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("andrei-stoica26/CSOA")
```

## Usage

The basic command to run CSOA is:

```
runCSOA(scObj, genes)
```

`scObj` must contain normalized and log-transformed gene expression data provided 
in one of the following formats:

- `Seurat`.
    - CSOA will use the expression matrix stored in the `data` layer.
- `SingleCellExperiment`
    - CSOA will use the expression matrix stored in the `logcounts` assay.
- `matrix`
- `dgCMatrix`

`genes` must be a character vector. 

For scoring multiple gene signatures using the `runCSOAMultiple` function, 
`geneSets` must be a list of character vectors and `geneSetNames` must be a 
character vector of the same length:


```
runCSOAMultiple(scObj, geneSets, geneSetNames)
```
