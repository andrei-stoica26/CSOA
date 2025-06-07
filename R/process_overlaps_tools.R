#' @importFrom spatstat.utils revcumsum
#'
NULL

#' Rank cell set overlaps
#'
#' This function ranks cell set overlaps by taking the average of the adjusted
#' p-value rank and the ratio of shared cells over expected shared cells rank.
#'
#' @param overlapDF A data frame created with generate_overlaps
#'
#' @return A data frame with ranked overlaps
#'
rankOverlaps <- function(overlapDF){
  overlapDF$pvalRank <- seq_len(nrow(overlapDF))
  overlapDF <- overlapDF[order(overlapDF$ratio, decreasing=TRUE), ]
  overlapDF$ratioRank <- seq_len(nrow(overlapDF))
  overlapDF$rank <- rowMeans(overlapDF[ , c('pvalRank', 'ratioRank')])
  overlapDF <- overlapDF[order(overlapDF$rank), ]
  return(overlapDF)
}

#' Score cell set overlaps
#'
#' This function scores cell set overlaps based on their ranks. The score of
#' the top overlap is set to 1, and the score decreases logarithmically towards
#' 0, which corresponds to the score of the first overlap not included in the
#' filtered data frame (i.e., if the filtered data frame contains 100 overlaps,
#' the 101st overlap corresponds to a score of 0)
#'
#' @inheritParams rankOverlaps
#'
#' @return A data frame with ranked overlaps
#'
scoreOverlaps <- function(overlapDF){
  overlapDF$rank <- seq_len(nrow(overlapDF))
  overlapDF$score <- log(seq(exp(1), 1, length.out = nrow(overlapDF) + 1)[overlapDF$rank])
  return(overlapDF)
}

#' Compute pair scores
#'
#' This function assesses the relative contribution of each gene pair to the CSOA
#' score
#'
#' @inheritParams computeCellScores
#' @param pairScores A list of pair scores in each cell for each pair in the
#' overlap data frame
#' @param pairFileName The name of the file where the pair data frame
#' will be saved. Default is NULL (the overlap data frame will not be saved)
#'
#' @return A data frame with overlap and pair scores and ranks
#'
#' @export
#'
computePairScores <- function(overlapDF, pairScores, pairFileName = NULL){
  df <- overlapDF[, c('gene1', 'gene2', 'score', 'rank')]
  colnames(df)[3:4] <- paste0('overlap', c('Score', 'Rank'))
  pairTotalScores <- colSums(data.frame(pairScores))
  totalScore <- sum(pairTotalScores)
  df$pairScore <- pairTotalScores / totalScore * 100
  df <- df[order(df$pairScore, decreasing = TRUE), ]
  df$pairRank <- seq_len(nrow(overlapDF))
  df$revCumsum <- spatstat.utils::revcumsum(df$pairScore)
  pairFile <- paste0(pairFileName, '.qs')
  message(paste0('Saving pair file: ', pairFile, '...'))
  qsave(df, pairFile)
  return(df)
}

#' Extract gene pairs from overlap matrix
#'
#' This function extracts the gene pairs from an overlap matrix
#'
#' @inheritParams rankOverlaps
#'
#' @return A list of gene pairs
#'
#' @export
#'
overlapPairs <- function(overlapDF)
  return(as.list(data.table::transpose(overlapDF[, c(1, 2)])))

#' Extract subset delineated using gene pairs from overlap matrix
#'
#' This function extracts the subset determined by input gene pairs from an
#' overlap matrix
#'
#' @inheritParams rankOverlaps
#' @param pairs Gene pairs corresponding to the extracted overlaps
#'
#' @return An overlap data frame corresponding to the selected gene pairs
#'
#' @export
#'
overlapSlice <- function(overlapDF, pairs)
  return(overlapDF[which(overlapPairs(overlapDF) %in% pairs),])
