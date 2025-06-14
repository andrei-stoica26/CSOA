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
  overlapDF <- overlapDF[order(overlapDF$pval), ]
  overlapDF$pvalRank <- rank(overlapDF$pval, ties.method = 'min')
  overlapDF <- overlapDF[order(overlapDF$ratio, decreasing = T), ]
  overlapDF$ratioRank <- rank(-overlapDF$ratio, ties.method = 'min')

  geneConn <- geneBestEdgeRank(overlapDF)
  overlapDF$pvalRank <- (geneConn[overlapDF$gene1, 1] + geneConn[overlapDF$gene2, 1]) / 2
  overlapDF$ratioRank <- (geneConn[overlapDF$gene1, 2] + geneConn[overlapDF$gene2, 2]) / 2
  overlapDF$rawAggRank <- (overlapDF$pvalRank + overlapDF$ratioRank) / 2

  overlapDF$rank <- rank(overlapDF$rawAggRank, ties.method = 'min')
  overlapDF <- overlapDF[order(overlapDF$rank), ]
  return(overlapDF)
}

#' Find the raw aggregate rank of the highest non-top overlap
#'
#' This function finds the raw aggregate rank of the highest non-top overlap.
#'
#' @param overlapDF A ranked overlap data frame
#'
#' @return A numeric value
#'
#' @export
#'
firstExcluded <- function(overlapDF){
  freqDF <- data.frame(rank = unique(overlapDF$rank), freq = as.numeric(table(overlapDF$rank)))
  freqSub <- subset(freqDF, freq == max(freq))
  rankCutoff <- mean(c(max(freqSub$rank), min(freqSub$rank)))

  outDF <- subset(overlapDF, rank > rankCutoff)
  if (nrow(outDF))
    return(outDF$rawAggRank[1])
  return(NULL)
}

#' Filter cell set overlaps
#'
#' This function filters cell set overlaps after the overlap data frame
#' has been ranked
#'
#' @inheritParams firstExcluded
#' @param firstOutRawRank The raw aggregate rank of the first overlap that will be
#' excluded
#'
#' @return A filtered overlap data frame
#'
#' @export
#'
filterOverlaps <- function(overlapDF, firstOutRawRank = NULL){
  if(is.null(firstOutRawRank))
    return(overlapDF)
  return(subset(overlapDF, rawAggRank < firstOutRawRank))
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
#' @param firstOutRawRank The raw rank of the highest-ranked overlap that was not
#' recorded among top overlaps. Ignored if osMethod is set to 'log'
#'
#' @return A data frame with ranked overlaps
#'
#' @export
#'
scoreOverlaps <- function(overlapDF, osMethod = 'log', firstOutRawRank = NULL){
  if (nrow(overlapDF) == 1){
    overlapDF$score <- 1
    return(overlapDF)
  }
  if(!osMethod %in% c('log', 'minmax'))
    stop('Unrecognized overlap scoring method. See ?CSOA::scoreOverlaps for the accepted methods')

  if (osMethod == 'log'){
    rankVals <- unique(overlapDF$rank)
    logVals <- log(seq(exp(1), 1, length.out = length(rankVals) + 1))[seq_along(rankVals)]
    names(logVals) <- rankVals
    overlapDF$score <- logVals[as.character(overlapDF$rank)]
  }

  if (osMethod == 'minmax'){
    rawRank <- overlapDF$rawAggRank
    if (!is.null(firstOutRawRank))
      rawRank <- c(rawRank, firstOutRawRank) else
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
