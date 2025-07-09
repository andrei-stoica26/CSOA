#' @importClassesFrom Seurat Seurat
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment assay
#' @importFrom stats setNames
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
generateOverlaps <- function(geneSetExp, percentile = 90, pairs = NULL,
                             overlapFileName = NULL){
  cellSets <- percentileSets(geneSetExp, percentile)
  overlapDF <- cellSetsOverlaps(cellSets, dim(geneSetExp)[2], pairs,
                                overlapFileName)
  return(overlapDF)
}

#' Process data frame of overlaps of cell sets
#'
#' This function filters, ranks and scores previously generated overlaps of cell
#' sets
#'
#' @inheritParams rankOverlaps
#' @inheritParams byCorrectDF
#' @inheritParams filterOverlaps
#' @inheritParams prepareFiltering
#' @param jaccardCutoff A cutoff used in the filtering of edges with low Jaccard
#' scores. NULL by default (no filtering of such edges will be performed).
#' @inheritParams scoreOverlaps
#'
#' @return A data frame consisting of filtered, ranked and scored cell sets
#' overlaps
#'
#' @export
#'
processOverlaps <- function(overlapDF, pvalThr = 0.05, saveCutoffPlot = FALSE, jaccardCutoff = NULL, osMethod = 'log'){

  if (nrow(overlapDF) > 1)
    overlapDF <- byCorrectDF(overlapDF, pvalThr) else overlapDF$pval_adj <- overlapDF$pval

  if (!nrow(overlapDF))
    return(overlapDF)

  overlapDF <- rankOverlaps(overlapDF)
  firstOutRawRank <- prepareFiltering(overlapDF, saveCutoffPlot)
  overlapDF <- filterOverlaps(overlapDF, firstOutRawRank)
  if (!is.null(jaccardCutoff)){
    overlapDF <- breakWeakTies(overlapDF, jaccardCutoff)
    firstOutRawRank <- NULL
  }

  overlapDF <- scoreOverlaps(overlapDF, osMethod, firstOutRawRank)
  return(overlapDF)
}

#' Assign a per-cell gene set score to Seurat object
#'
#' This function uses the scored data frame of overlaps to compute a CSOA score
#' for each cell in a Seurat object
#'
#' @inheritParams computePCPairScores
#' @inheritParams computePairScores
#' @inheritParams computePCSetScores
#'
#' @return A Seurat object with a CSOA score assigned for each cell
#'
#' @export
#'
computeCellScores <- function(overlapDF, normExp, cellNames, colStr='CSOA', pairFileName = NULL, keepOverlapOrder = FALSE){
  pcPairScores <- computePCPairScores(overlapDF, normExp)
  if(!is.null(pairFileName))
    pairScores <- computePairScores(overlapDF, pcPairScores, pairFileName, keepOverlapOrder)
  scoreDF <- computePCSetScores(pcPairScores, cellNames, colStr)
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
  for (colName in colnames(scoreDF))
    if (colName %in% colnames(scObj@meta.data))
      scObj@meta.data[[colName]] <- c()
  scObj@meta.data <- cbind(scObj@meta.data, scoreDF)
  return(scObj)
}

#' @param altExpName Name of the matrix storing CSOA scores
#'
#' @rdname storeCellScores
#' @export
#'
storeCellScores.SingleCellExperiment <- function(scObj, scoreDF, altExpName = 'CSOA', ...){
  for (colName in colnames(scoreDF))
    if (colName %in% colnames(colData(scObj)))
      colData(scObj)[[colName]] <- c()
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
scoreCells <- function(geneSetExp, overlapDF, colStr = 'CSOA', pvalThr = 0.05,
                       saveCutoffPlot = FALSE, jaccardCutoff = NULL, osMethod = 'log',
                       pairFileName = NULL, keepOverlapOrder = FALSE){
  overlapDF <- processOverlaps(overlapDF, pvalThr, saveCutoffPlot, jaccardCutoff,
                               osMethod)

  if(!nrow(overlapDF)){
    warning('No significant overlaps were identified. All cells will get a score of 0.')
    scoreDF <- data.frame(setNames(list(rep(0, dim(geneSetExp)[2])), colStr))
    return(scoreDF)
  }

  message('Normalizing expression matrix by rows...')
  genes <- overlapGenes(overlapDF)
  normExp <- kerntools::minmax(geneSetExp[genes, ], rows=TRUE)
  scoreDF <- computeCellScores(overlapDF, normExp, colnames(geneSetExp), colStr,
                               pairFileName, keepOverlapOrder)
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
#
#' @return An object of the same class as scObj with a CSOA score assigned for each cell
#'
#' @export
#'
runCSOA <- function(scObj, genes, colStr='CSOA', percentile = 90,
                    overlapFileName = NULL, pvalThr = 0.05,
                    saveCutoffPlot = FALSE, jaccardCutoff = NULL, osMethod = 'log',
                    pairFileName = NULL, keepOverlapOrder = FALSE){
  if (!min(is(genes)[1:2] == c('character', 'vector')) | length(genes) < 2)
    stop('genes must be a character vector of length >= 2.')
  geneSetExp <- expMat(scObj, genes)
  overlapDF <- generateOverlaps(geneSetExp, percentile, pairs=NULL,
                                overlapFileName)
  scoreDF <- scoreCells(geneSetExp, overlapDF, colStr, pvalThr, saveCutoffPlot,
                        jaccardCutoff, osMethod,
                        pairFileName, keepOverlapOrder)
  return(storeCellScores(scObj, scoreDF))
}
