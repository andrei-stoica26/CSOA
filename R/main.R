#' @importClassesFrom Seurat Seurat
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SingleCellExperiment altExp
#' @importFrom SummarizedExperiment assay
#' @importFrom kerntools minmax
#' @importFrom stats setNames
#' @include percentile_sets.R
#' @include generics.R
#' @include generate_overlaps_tools.R
#' @include process_overlaps_tools.R
NULL

#' Generate overlaps of cell sets for input genes
#'
#' This function builds a set of cells for each input gene by first selecting all
#' the cells expressing the gene and then retaining only the cells expressing the
#' gene at the input percentile. Subsequently, overlaps of pairs of the constructed
#' cell sets are assessed for statistical significance.
#'
#'
#' @inheritParams percentileSets
#' @inheritParams cellSetsOverlaps
#'
#' @return A data frame listing statistics for all cell set overlaps
#'
#' @export
#'
generateOverlaps <- function(geneSetExp, percentile = 90, pairs = NULL, overlapFileName = NULL){
  cellSets <- percentileSets(geneSetExp, percentile)
  overlapDF <- cellSetsOverlaps(cellSets, dim(geneSetExp)[2], pairs, overlapFileName)
  return(overlapDF)
}

#' Process data frame of overlaps of cell sets
#'
#' This function filters, ranks and scores previously generated overlaps of cell
#' sets
#'
#' @inheritParams rankOverlaps
#' @param nPairs Number of overlaps that will be retained
#' @inheritParamas byCorrectDF
#'
#' @return A data frame consisting of filtered, ranked and scored cell sets
#' overlaps
#'
#' @export
#'
processOverlaps <- function(overlapDF, nPairs = 100, pvalThr = 0.05){
  if (nrow(overlapDF) > 1)
    overlapDF <- byCorrectDF(overlapDF, pvalThr) else overlapDF$pval_adj <- overlapDF$pval
  overlapDF <- rankOverlaps(overlapDF)
  if (!is.null(nPairs)){
    if(nPairs > nrow(overlapDF))
      message(paste0('Will return only ', nrow(overlapDF), ' significant overlaps. More are not available.'))
    overlapDF <- overlapDF[seq_len(min(nPairs, nrow(overlapDF))), ]
  }
  overlapDF <- scoreOverlaps(overlapDF)
  return(overlapDF)
}

#' Assign a per-cell gene set score to Seurat object
#'
#' This function uses the scored data frame of overlaps to compute a CSOA score
#' for each cell in a Seurat object
#'
#' @param overlapDF Processed overlap dataframe generated using processOverlaps
#' @param normExp A minmax-normalized by rows matrix of gene expression
#' @param cellNames Cell names
#' @param colStr The name of the column where CSOA results will be stored
#'
#' @return A Seurat object with a CSOA score assigned for each cell
#'
#' @export
#'
computeCellScores <- function(overlapDF, normExp, cellNames, colStr='CSOA'){
  message('Computing per-cell scores for gene pairs...')
  pairsScores <- lapply(1:nrow(overlapDF), function(i) overlapDF[i, 'score'] *
                          normExp[overlapDF[i, 'gene1'], ] * normExp[overlapDF[i, 'gene2'], ])
  message('Computing per-cell pathway scores...')
  scores <- rowSums(data.frame(pairsScores))
  scores <- kerntools::minmax(as.matrix(scores))
  scoreDF <- data.frame(setNames(list(scores), colStr))
  rownames(scoreDF) <- cellNames
  return(scoreDF)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for CSOA-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param scoreDF Dataframe of CSOA scores
#'
#' @rdname storeCellScores
#' @export
#'
storeCellScores.default <- function(scObj, scoreDF, ...)
  stop('Unrecognized input type: scObj must be a Seurat object with a data assay, a SingleCellExperiment with a logcounts assay or a matrix.')

#' @rdname storeCellScores
#' @export
#'
storeCellScores.Seurat <- function(scObj, scoreDF, ...){
  scObj@meta.data <- cbind(scObj@meta.data, scoreDF)
  return(scObj)
}

#' @param altExpName Name of the matrix storing CSOA scores
#'
#' @rdname storeCellScores
#' @export
#'
storeCellScores.SingleCellExperiment <- function(scObj, scoreDF, altExpName = 'CSOA', ...){
  colData(scObj) <- cbind(colData(scObj), scoreDF)
  return(scObj)
}

#' @rdname storeCellScores
#' @export
#'
storeCellScores.matrix <- function(scObj, scoreDF, ...)
  return(scoreDF)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Generate CSOA scores from overlap data frame
#'
#' This function computes per-cell CSOA scores after the overlap data frame has
#' been generated.
#'
#' @inheritParams generateOverlaps
#' @inheritParams processOverlaps
#' @inheritParams computeCellScores
#'
#' @return A data frame with a column corresponding to the CSOA scores
#'
#' @export
#'
scoreCells <- function(geneSetExp, overlapDF, colStr = 'CSOA', nPairs = 100){
  overlapDF <- processOverlaps(overlapDF, nPairs)
  message('Normalizing expression matrix by rows...')
  genes <- overlapGenes(overlapDF)
  normExp <- kerntools::minmax(geneSetExp[genes, ], rows=T)
  scoreDF <- computeCellScores(overlapDF, normExp, colnames(geneSetExp), colStr)
  return(scoreDF)
}

#' Run the CSOA pipeline
#'
#' This function generates cell set overlaps for an input gene set based on
#' percentiles of gene expression, computes the significance of these overlaps,
#' ranks, filters and scores the overlaps based on this significance, and builds
#' a per-cell score by summing the products of the scores of these overlaps and
#' the custom-normalized per-cell expressions of the corresponding pairs of genes.
#'
#' @inheritParams expMat
#' @param genes Vector of genes. Must include at least two genes
#' @inheritParams generateOverlaps
#' @inheritParams scoreCells
#'
#' @return An object of the same class as scObj with a CSOA score assigned for each cell
#'
#' @export
#'
runCSOA <- function(scObj, genes, colStr='CSOA', percentile = 90, nPairs = 100, overlapFileName = NULL){
  if (!min(is(genes)[1:2] == c('character', 'vector')) | length(genes) < 2)
    stop('genes must be a character vector of length >= 2')
  genes <- sort(genes)
  geneSetExp <- expMat(scObj, genes)
  overlapDF <- generateOverlaps(geneSetExp, percentile, pairs = NULL, overlapFileName)
  scoreDF <- scoreCells(geneSetExp, overlapDF, colStr, nPairs)
  return(storeCellScores(scObj, scoreDF))
}


