# CSOA
Cell Set Overlap Analysis (CSOA) is a tool for calculating per-cell gene 
signature scores in a scRNA-seq dataset. For each gene in the signature, 
CSOA constructs a set consisting of cells highly expressing the gene, 
fter which all overlaps of pairs of cell sets are computed, ranked, 
filtered and scored. The CSOA per-cell score is calculated by summing up, 
over all top cell set overlaps, the product of the overlap score and the 
min-max-normalized expression of each of the two genes. CSOA can run on a Seurat object, a SingleCellExperiment object, a matrix and a dgCMatrix.

## Installation

To install CSOA, run the following R code:

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    
BiocManager::install("CSOA")
```
