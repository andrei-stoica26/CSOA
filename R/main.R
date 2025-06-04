#' @importClassesFrom Seurat Seurat
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SingleCellExperiment altExp
#' @importFrom SummarizedExperiment assay
#' @importFrom kerntools minmax
#' @importFrom stats setNames
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
#' @inheritParams fullQuantileSets
#' @inheritParams cellSetsOverlaps
#'
#' @return A data frame listing statistics for all cell set overlaps
#'
#' @export
#'
generateOverlaps <- function(scObj, genes, nQuantiles, pairs=NULL, overlapFileName=NULL){
  allCellSets <- fullQuantileSets(scObj, genes, nQuantiles)
  topCellSets <- topQuantileSets(allCellSets)
  overlapDF <- cellSetsOverlaps(topCellSets, dim(scObj)[2], pairs, overlapFileName)
  return(overlapDF)
}

#' Process data frame of overlaps of cell sets
#'
#' This function filters, ranks and scores previously generated overlaps of cell
#' sets
#'
#' @inheritParams rankOverlaps
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
  stop('Unrecognized input type: scObj must be a Seurat object with a data assay, a SingleCellExperiment with a logcounts assay or an expression matrix.')

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
  SingleCellExperiment::altExp(scObj, altExpName) <- SummarizedExperiment::SummarizedExperiment(t(scoreDF))
  return(scObj)
}

#' @rdname storeCellScores
#' @export
#'
storeCellScores.matrix <- function(scObj, scoreDF, ...){
  return(scoreDF)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Generate CSOA scores from overlap data frame
#'
#' This function generates per-cell CSOA scores after the overlap data frame has
#' been generated.
#'
#' @param expression Expression matrix
#' @inheritParams processOverlaps
#' @inheritParams generateOverlaps
#' @inheritParams computeCellScores
#'
#' @return A data frame with a column corresponding to the CSOA scores
#'
#' @export
#'
scoreCells <- function(expression, overlapDF, nPairs=100, colStr='CSOA'){
  overlapDF <- processOverlaps(overlapDF, nPairs)
  message('Normalizing expression matrix by rows...')
  genes <- overlapGenes(overlapDF)
  normExp <- kerntools::minmax(expression[genes, ], rows=T)
  scoreDF <- computeCellScores(overlapDF, normExp, colnames(expression), colStr)
  return(scoreDF)
}

#' Run the CSOA pipeline
#'
#' This function generates cell set overlaps for an input gene set based on
#' quantiles of gene expression and retaining the top quantile genes.
#'
#' @inheritParams generateOverlaps
#' @inheritParams scoreCells
#'
#' @return An object of the same class as scObj with a CSOA score assigned for each cell
#'
#' @export
#'
runCSOA <- function(scObj, genes, colStr='CSOA', nQuantiles=10, nPairs=100, overlapFileName=NULL){
  expression <- expMat(scObj)
  overlapDF <- generateOverlaps(expression, genes, nQuantiles, overlapFileName)
  scoreDF <- scoreCells(expression, overlapDF, nPairs, colStr)
  return(storeCellScores(scObj, scoreDF))
}


