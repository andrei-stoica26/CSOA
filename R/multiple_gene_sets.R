#' Generate CSOA scores from overlap data frame for multiple gene sets
#'
#' This function scores an overlap data frame generated using multiple
#' gene sets.
#' The overlap data frame is split based on the overlaps corresponding to
#' each gene set and scored, and the output is rejoined as a data frame.
#'
#' @inheritParams scoreCells
#' @param setPairs A list of overlaps corresponding to each input gene set.
#' @param geneSetNames Character vector of names of gene sets.
#' @param pairFileTemplate Character object used in the naming of the files
#' where the pair data frames will be saved. Default is NULL (the pair
#' data frames will not be saved).
#' @param keepOverlapOrder Keep the rank-based order of overlaps in the
#' pair score file, as opposed to changing it to a pair score-based order.
#' Ignored if pairFileTemplate is NULL.
#'
#' @return A data frame whose columns correspond to the CSOA scores of the
#' input gene sets.
#'
#' @export
#'
scoreCellsMultiple <- function(geneSetExp,
                               overlapDF,
                               setPairs,
                               geneSetNames,
                               pvalThr = 0.05,
                               saveCutoffPlot = FALSE,
                               jaccardCutoff = NULL,
                               osMethod = 'log',
                               pairFileTemplate = NULL,
                               keepOverlapOrder = FALSE){
  if(!is.null(pairFileTemplate))
    pairFileName <- paste0(pairFileTemplate, geneSetNames) else
      pairFileName <- NULL
  scoreDFList <- lapply(seq_along(setPairs), function(i) {
    setOverlapDF <- overlapSlice(overlapDF, setPairs[[i]])
    scoreDF <- scoreCells(geneSetExp, setOverlapDF, geneSetNames[i],
                          pvalThr, saveCutoffPlot, jaccardCutoff,
                          osMethod, pairFileName[i], keepOverlapOrder)
    return(scoreDF)
  })
  allScoresDF <- Reduce(cbind, scoreDFList)
  return(allScoresDF)
}

#' Run the CSOA pipeline for multiple gene sets
#'
#' This function generates cell set overlaps for input gene sets
#' based on percentiles of gene expression, computes the significance
#' of these overlaps, ranks, filters and scores the overlaps based on
#' significance, and builds per-cell score by summing the products of
#' the scores of these overlaps and the custom-normalized per-cell
#' expressions of the corresponding pairs of genes.
#'
#' @inheritParams runCSOA
#' @param geneSets List of character vectors
#' @inheritParams scoreCellsMultiple
#'
#' @return An object of the same class as scObj with per-gene-set CSOA scores
#' assigned for each cell
#'
#' @export
#'
runCSOAMultiple <- function(scObj,
                            geneSets,
                            geneSetNames,
                            percentile = 90,
                            overlapFileName = NULL,
                            pvalThr = 0.05,
                            saveCutoffPlot = FALSE,
                            jaccardCutoff = NULL,
                            osMethod = 'log',
                            pairFileTemplate = NULL,
                            keepOverlapOrder = FALSE){
  geneSets <- lapply(geneSets, sort)
  setPairs <- lapply(geneSets, getPairs)
  pairs <- Reduce(union, setPairs)
  genes <- Reduce(union, geneSets)
  geneSetExp <- expMat(scObj, genes)
  overlapDF <- generateOverlaps(geneSetExp, percentile,
                                pairs, overlapFileName)
  scoreDF <- scoreCellsMultiple(geneSetExp, overlapDF,
                                setPairs, geneSetNames,
                                pvalThr, saveCutoffPlot,
                                jaccardCutoff, osMethod,
                                pairFileTemplate, keepOverlapOrder)
  return(storeCellScores(scObj, scoreDF))
}
