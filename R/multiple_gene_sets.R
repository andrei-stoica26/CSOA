#' Join scores for multiple gene sets
#'
#' This function joins multiple data frames with CSOA scores and returns the
#' scored object
#'
#' @inheritParams runCSOA
#' @param scoreDFList A list of data frames with CSOA scores
#'
#' @return An object (Seurat, SingleCellExpression or matrix, depending on the
#' input) containing all the scores
#'
joinCellScores <- function(scObj, scoreDFList){
  allScoresDF <- Reduce(cbind, scoreDFList)
  return(storeCellScores(scObj, allScoresDF))
}

#' Run the CSOA pipeline for multiple gene sets
#'
#' This function generates cell set overlaps for input gene sets based on
#' quantiles of gene expression and retaining the top quantile genes.
#'
#' @inheritParams runCSOA
#' @param geneSets List of character vectors
#' @param geneSetsNames Character vector of names of gene sets
#'
#' @return An object of the same class as scObj with per-gene-set CSOA scores
#' assigned for each cell
#'
#' @export
#'
runCSOAMultiple <- function(scObj, geneSets, geneSetsNames, nQuantiles=10, nPairs=100, overlapFileName=NULL){
  expression <- expMat(scObj)
  geneSets <- lapply(geneSets, sort)
  genes <- unique(unlist(geneSets))
  setPairs <- lapply(geneSets, getPairs)
  pairs <- Reduce(union, setPairs)
  overlapDF <- generateOverlaps(expression, genes, nQuantiles, pairs, overlapFileName)
  scoreDFList <- lapply(seq_along(setPairs), function(i) {
    setOverlapDF <- overlapSlice(overlapDF, setPairs[[i]])
    setGenes <- geneSets[i]
    scoreDF <- scoreCells(expression, setOverlapDF, nPairs, geneSetsNames[i])
    return(scoreDF)
  })
  return(joinCellScores(scObj, scoreDFList))
}
