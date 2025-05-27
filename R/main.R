#' @importClassesFrom Seurat Seurat
#' @importFrom kerntools minmax
#' @include quantile_sets.R
#' @include generics.R
#' @include generate_overlaps_tools.R
#' @include process_overlaps_tools.R
NULL

#' Generate overlaps of cell sets for input genes
#'
#' This function builds a set of cells for each input gene by dividing the cells
#' expressing the gene into quantiles based on their expression and retaining
#' the top-quantile cells. Subsequently, the significance of the overlaps of all
#' pairs of the constructed cell sets is computed, and the overlaps receive a
#' ranking.
#'
#' @param expObj A Seurat object or expression matrix of matrix class
#' @param genes Vector of genes. Must include at least two genes
#' @param nQuantiles An integer between 2 and 10
#' @param overlapFileName A character used to save the overlap data frame.
#' Default is NULL (the overlap data frame will not be saved)
#'
#' @return A data frame listing statistics for all cell set overlaps
#'
#' @export
#'
generateOverlaps <- function(expObj, genes, nQuantiles, overlapFileName=NULL){
  allCellSets <- fullQuantileSets(expObj, genes, nQuantiles)
  topCellSets <- topQuantileSets(allCellSets)
  overlapDF <- cellSetsOverlaps(topCellSets, dim(expObj)[2], overlapFileName)
  return(overlapDF)
}

#' Process data frame of overlaps of cell sets
#'
#' This function filters, ranks and scores previously generated overlaps of cell
#' sets
#'
#' @param overlapDF A data frame created with generateOverlaps
#' @param nPairs Number of overlaps that will be retained
#'
#' @return A data frame consisting of filtered, ranked and scored cell sets
#' overlaps
#'
#' @export
#'
processOverlaps <- function(overlapDF, nPairs=100){
  overlapDF <- subset(overlapDF, pval_adj < 0.05)
  overlapDF <- rankOverlaps(overlapDF)
  if (!is.null(nPairs)){
    if(nPairs > nrow(overlapDF))
      message(str_c('Will return only ', nrow(overlapDF), ' significant overlaps. More are not available.'))
    overlapDF <- overlapDF[1:min(nPairs, nrow(overlapDF)), ]
  }
  overlapDF <- scoreOverlaps(overlapDF)
  return(overlapDF)
}

#' Assign a per-cell gene set score to Seurat object
#'
#' This function uses the scored data frame of overlaps to compute a CSOA score
#' for each cell in a Seurat obj. It also requires
#'
#' @param seuratObj A Seurat object
#' @param overlapDF A data frame created with generateOverlaps
#' @param normExp A minmax-normalized by rows matrix of gene expression
#' @param colStr The name of the column where CSOA results will be stored
#'
#' @return A Seurat object with a CSOA score assigned for each cell
#'
#' @export
#'
scoreSeurat <- function(seuratObj, overlapDF, normExp, colStr='CSOA'){
  message('Computing per-cell scores for gene pairs...')
  pairsScores <- lapply(1:nrow(overlapDF), function(i) overlapDF[i, 'score'] *
                          normExp[overlapDF[i, 'gene1'], ] * normExp[overlapDF[i, 'gene2'], ])
  message('Computing per-cell pathway scores...')
  scores <- rowSums(data.frame(pairsScores))
  seuratObj@meta.data[[colStr]] <- scores / max(scores)
  return(seuratObj)
}

#' Run the CSOA pipeline
#'
#' This function generates cell set overlaps for an input gene set based on
#' quantiles of gene expression and retaining the top quantile genes.
#'
#' @param seuratObj A Seurat object
#' @param genes Vector of genes. Must include at least two genes
#' @param nQuantiles An integer between 2 and 10
#' @param nPairs Number of overlaps that will be retained
#' @param colStr The name of the column where CSOA results will be stored
#' @param overlapFileName The name of the file where the overlap data frame
#' @param overlapFileName The name of the file where the overlap data frame
#' will be saved. Default is NULL (the overlap data frame will not be saved)
#'
#' @return A Seurat object with a CSOA score assigned for each cell
#'
#' @export
#'
runCSOA <- function(seuratObj, genes, nQuantiles=10, nPairs=100, colStr='CSOA', overlapFileName=NULL){
  expObj <- expMat(seuratObj)
  overlapDF <- generateOverlaps(expObj, genes, nQuantiles, overlapFileName)
  overlapDF <- processOverlaps(overlapDF, nPairs)
  message('Normalizing expression matrix by rows...')
  normExp <- kerntools::minmax(expObj[genes, ], rows=T)
  seuratObj <- scoreSeurat(seuratObj, overlapDF, normExp, colStr)
  return(seuratObj)
}
