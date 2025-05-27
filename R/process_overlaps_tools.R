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
  overlapDF$pvalRank <- 1:nrow(overlapDF)
  overlapDF <- overlapDF[order(overlapDF$ratio, decreasing=T), ]
  overlapDF$ratioRank <- 1:nrow(overlapDF)
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
#' @param overlapDF A data frame created with generate_overlaps
#'
#' @return A data frame with ranked overlaps
#'
scoreOverlaps <- function(overlapDF){
  overlapDF$rank <- 1:nrow(overlapDF)
  overlapDF$score <- log(seq(exp(1), 1, length.out = nrow(overlapDF) + 1)[1:nrow(overlapDF)])
  return(overlapDF)
}
