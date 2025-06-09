#' Compute per-cell gene pair scores
#'
#' This function scores each gene pair corresponding to a top overlap in each
#' cell
#'
#' @param overlapDF A scored overlap data frame
#' @param normExp A min-max normalized expression matrix of the genes involved in
#' top overlaps
#'
#' @return A data frame with the same dimensions as normExp
#'
#' @export
#'
computePCPairScores <- function(overlapDF, normExp){
  message('Computing per-cell scores for gene pairs...')
  df <- data.table::transpose(data.frame(lapply(1:nrow(overlapDF), function(i) overlapDF[i, 'score'] *
                         normExp[overlapDF[i, 'gene1'], ] * normExp[overlapDF[i, 'gene2'], ])))
  rownames(df) <- paste0(paste0(overlapDF$gene1, '_'), overlapDF$gene2)
  colnames(df) <- colnames(normExp)
  return(df)
}

#' Compute aggregate gene pair scores
#'
#' This function assesses the relative contribution of each gene pair to the CSOA
#' score
#'
#' @inheritParams computeCellScores
#' @param pcPairScores A date frame of pair scores in each cell for each pair
#' in the overlap data frame
#' @param pairFileName The name of the file where the pair data frame
#' will be saved. Default is NULL (the pair data frame will not be saved)
#' @param keepOverlapOrder Keep the rank-based order of overlaps in the pair score
#' file, as opposed to changing it to a pair score-based order. Ignored if
#' pairFileName is NULL
#'
#' @return A data frame with overlap and pair scores and ranks
#'
#' @export
#'
computePairScores <- function(overlapDF, pcPairScores, pairFileName = NULL, keepOverlapOrder = FALSE){
  df <- overlapDF[, c('gene1', 'gene2', 'score', 'rank')]
  colnames(df)[3:4] <- paste0('overlap', c('Score', 'Rank'))

  pairTotalScores <- rowSums(pcPairScores)
  totalScore <- sum(pairTotalScores)

  df$pairScore <- pairTotalScores / totalScore * 100
  df <- df[order(df$pairScore, decreasing = TRUE), ]

  df$pairRank <- seq_len(nrow(overlapDF))
  df$revCumsum <- spatstat.utils::revcumsum(df$pairScore)

  if (keepOverlapOrder)
    df <- df[order(df$overlapScore, decreasing = TRUE), ]

  if (!is.null(pairFileName)){
    pairFile <- paste0(pairFileName, '.qs')
    message(paste0('Saving pair scores file: ', pairFile, '...'))
    qsave(df, pairFile)
  }

  return(df)
}

#' Aggregate per-cell gene pair scores
#'
#' This function aggregates per-cell gene pair scores into per-cell gene
#' signature scores
#'
#' @inheritParams computePairScores
#' @param cellNames Cell names
#' @param colStr The name of the column where CSOA results will be stored
#'
#' @return A data frame with the per-cell gene signature score as a column
#'
#' @export
#'
computePCSetScores <- function(pcPairScores, cellNames, colStr = 'CSOA'){
  message('Computing per-cell gene signature scores...')
  scores <- colSums(pcPairScores)
  scores <- vMinmax(scores)
  scoreDF <- data.frame(setNames(list(scores), colStr))
  rownames(scoreDF) <- cellNames
  return(scoreDF)
}


