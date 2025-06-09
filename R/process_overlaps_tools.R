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
#' @export
#'
rankOverlaps <- function(overlapDF){
  overlapDF$pvalRank <- seq_len(nrow(overlapDF))
  overlapDF <- overlapDF[order(overlapDF$ratio, decreasing=TRUE), ]
  overlapDF$ratioRank <- seq_len(nrow(overlapDF))
  overlapDF$rawAggRank <- rowMeans(overlapDF[ , c('pvalRank', 'ratioRank')])
  overlapDF <- overlapDF[order(overlapDF$rawAggRank), ]
  overlapDF$rank <- seq_len(nrow(overlapDF))
  return(overlapDF)
}

#' Find the raw aggregate rank of the highest non-top overlap
#'
#' This function finds the raw aggregate rank of the highest non-top overlap.
#'
#' @param overlapDF A ranked overlap data frame
#' @param nPairs Number of overlaps that will be retained
#'
#' @return A numeric value
#'
#' @export
#'
firstExcluded <- function(overlapDF, nPairs=100){
  if (nrow(overlapDF) > nPairs)
    return(overlapDF$rawAggRank[nPairs + 1])
  return(NULL)
}

#' Filter cell set overlaps
#'
#' This function filters cell set overlaps after the overlap data frame
#' has been ranked
#'
#' @inheritParams firstExcluded
#'
#' @return A filtered overlap data frame
#'
#' @export
#'
filterOverlaps <- function(overlapDF, nPairs = 100){
  if (!is.null(nPairs)){
    if(nPairs > nrow(overlapDF))
      message(paste0('Will return only ', nrow(overlapDF), ' significant overlaps. More are not available.'))
    overlapDF <- overlapDF[seq_len(min(nPairs, nrow(overlapDF))), ]
  }
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
#' @param overlapDF A filtered overlap data frame
#' @param osMethod The method used to compute overlap scores. Options are 'log'
#' and 'minmax'
#' @param firstOutRank The raw rank of the highest-ranked overlap that was not
#' recorded among top overlaps. Ignored if osMethod is set to 'log'
#'
#' @return A data frame with ranked overlaps
#'
#' @export
#'
scoreOverlaps <- function(overlapDF, osMethod = 'log', firstOutRank = NULL){
  if (nrow(overlapDF) == 1){
    overlapDF$score <- 1
    return(overlapDF)
  }
  if(!osMethod %in% c('log', 'minmax'))
    stop('Unrecognized overlap scoring method. See ?CSOA::scoreOverlaps for the accepted methods')

  if (osMethod == 'log')
    overlapDF$score <- log(seq(exp(1), 1, length.out = nrow(overlapDF) + 1)[overlapDF$rank])
  if (osMethod == 'minmax'){
    rawRank <- overlapDF$rawAggRank
    if (!is.null(firstOutRank))
      rawRank <- c(rawRank, firstOutRank) else
        rawRank <- c(rawRank, 2 * rawRank[nrow(overlapDF)] - rawRank[nrow(overlapDF) - 1])
    overlapDF$score <- (1 - vMinmax(rawRank))[overlapDF$rank]
  }
  return(overlapDF)
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
