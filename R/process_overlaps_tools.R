#' @importFrom spatstat.utils revcumsum
#'
NULL

#' Find the connectivity of each gene
#'
#' This function finds the connectivity of each gene
#' from an overlap data frame ranked using p-values
#' and recorded-over-expected ratios.
#'
#' @param overlapDF Overlap data frame with the pvalRank
#' and ratioRank columns
#' @param asRanks Whether to replace connectivity
#' scores by ranks
#'
#' @return A data frame with genes involved in the overlaps
#' as row names, and two columns, corresponding to connectivity
#' ranks (by default) or scores (if asRanks is set to FALSE)
#' for both p-value and recorded-over-observed size ratio.
#'
#' @noRd
#'
geneBestEdgeRank <- function(overlapDF, asRanks = TRUE){
  genes <- overlapGenes(overlapDF)
  mat <- do.call(rbind, lapply(genes, function(gene){
    geneDF <- overlapDF[overlapDF$gene1 == gene |
                          overlapDF$gene2 == gene, ]
    return(c(min(geneDF$pvalRank), min(geneDF$ratioRank)))
  }))
  rownames(mat) <- genes
  colnames(mat) <- c('connPvalRank', 'connRatioRank')
  if (asRanks){
    mat <- rankReplace(mat, 'connPvalRank')
    mat <- rankReplace(mat, 'connRatioRank')
  }
  return(mat)
}

#' Rank cell set overlaps
#'
#' This function ranks cell set overlaps by taking the average of
#' the adjusted p-value rank and the ratio of shared cells over
#' expected shared cells rank.
#'
#' @param overlapDF A data frame created with generate_overlaps
#'
#' @return A data frame with ranked overlaps
#'
#' @export
#'
#' @examples
#' overlapDF <- data.frame(gene1=paste0('G', c(1, 3, 7, 6, 8, 2, 4, 3, 4, 5)),
#' gene2=paste0('G', c(2, 7, 2, 5, 4, 5, 1, 2, 2, 8)),
#' ratio=runif(10, 2, 10),
#' pval=runif(10, 0, 1e-10))
#' rankOverlaps(overlapDF)
#'
rankOverlaps <- function(overlapDF){
  if (!nrow(overlapDF))
    return(overlapDF)
  overlapDF <- overlapDF[order(overlapDF$pval), ]
  overlapDF$pvalRank <- rankFun(overlapDF$pval)
  overlapDF <- overlapDF[order(overlapDF$ratio,
                               decreasing=TRUE), ]
  overlapDF$ratioRank <- rankFun(-overlapDF$ratio)

  geneConn <- geneBestEdgeRank(overlapDF)
  overlapDF$pvalRank <- (geneConn[overlapDF$gene1, 1] +
                           geneConn[overlapDF$gene2, 1]) / 2
  overlapDF$ratioRank <- (geneConn[overlapDF$gene1, 2] +
                            geneConn[overlapDF$gene2, 2]) / 2
  overlapDF$rawAggRank <- (overlapDF$pvalRank +
                             overlapDF$ratioRank) / 2

  overlapDF$rank <- rankFun(overlapDF$rawAggRank)
  overlapDF <- overlapDF[order(overlapDF$rank), ]
  return(overlapDF)
}

#' Find overlap rank cutoff
#'
#' This function finds the cutoff for rank-based
#' filtering of overlaps.
#'
#' @param freqDF A frequency data frame of overlap ranks.
#'
#' @return Rank cutoff.
#'
#' @export
#'
#' @examples
#' freqDF <- data.frame(rank = c(1, 2, 4, 7),
#' n = c(1, 3, 3, 2))
#' findRankCutoff(freqDF)
#'
findRankCutoff <- function(freqDF){
  freqSub <- subset(freqDF, n == max(n))
  rankCutoff <- mean(c(max(freqSub$rank), min(freqSub$rank)))
  return(rankCutoff)
}

#' Find the raw aggregate rank of the highest non-top overlap
#'
#' This function finds the raw aggregate rank of the
#' highest non-top overlap.
#'
#' @param overlapDF A ranked overlap data frame
#' @param saveCutoffPlot Whether to save overlap cutoff plot
#'
#' @return A numeric value
#'
#' @export
#'
#' @examples
#' overlapDF <- data.frame(gene1 = paste0('G', c(1, 2, 3, 4, 7, 7)),
#' gene2 = paste0('G', c(2, 5, 1, 8, 4, 9)),
#' rawAggRank = c(7, 9, 9, 11.5, 11.5, 13),
#' rank = c(1, 2, 2, 4, 4, 6))
#' prepareFiltering(overlapDF)
#'
#'
prepareFiltering <- function(overlapDF, saveCutoffPlot = FALSE){
  if (!nrow(overlapDF))
    return(NULL)
  freqDF <- dplyr::count(overlapDF, rank)
  rankCutoff <- findRankCutoff(freqDF)

  if(saveCutoffPlot){
    if (nrow(freqDF) < 2)
      stop('overlapCutoffPlot requires at least two points.')

    freqs <- freqDF$n
    nFreq <- length(unique(freqs))

    if (nFreq < 2)
      stop('overlapCutoffPlot requires at least two',
           'distinct rank frequencies.')

    maxFreq <- max(freqs)
    maxApps <- which(freqs %in% maxFreq)
    if (length(maxApps) == 1)
      if (maxApps %in% c(1, nrow(freqDF)))
        stop('overlapCutoffPlot requires that the',
        'maximum rank frequency is reached at a non-extremal point.')

    devPlot(overlapCutoffPlot, freqDF, rankCutoff)
  }

  outDF <- subset(overlapDF, rank > rankCutoff)
  if (nrow(outDF))
    firstOutRawRank <- outDF$rawAggRank[1] else
      if (nrow(overlapDF) > 1)
        firstOutRawRank <- 2 * overlapDF$rawAggRank[nrow(overlapDF)] -
    overlapDF$rawAggRank[nrow(overlapDF) - 1] else
          firstOutRawRank <- NULL

  return(firstOutRawRank)
}

#' Filter cell set overlaps
#'
#' This function filters cell set overlaps after the overlap data frame
#' has been ranked
#'
#' @inheritParams prepareFiltering
#' @param firstOutRawRank The raw aggregate rank of the
#' first overlap that will be excluded
#'
#' @return A filtered overlap data frame
#'
#' @examples
#' overlapDF <- data.frame(gene1 = paste0('G', c(1, 2, 3, 4, 7)),
#' gene2 = paste0('G', c(2, 5, 1, 8, 4)),
#' rawAggRank = c(1, 2, 2, 4, 4))
#' filterOverlaps(overlapDF, 2)
#'
#'
#' @export
#'
filterOverlaps <- function(overlapDF, firstOutRawRank = NULL){
  if(is.null(firstOutRawRank) | !nrow(overlapDF))
    return(overlapDF)
  return(subset(overlapDF, rawAggRank < firstOutRawRank))
}

#' Score cell set overlaps
#'
#' This function scores cell set overlaps based on their ranks.
#' The score of the top overlap is set to 1, and the score decreases
#' logarithmically towards 0, which corresponds to the score of the
#' first overlap not included in the filtered data frame (i.e., if
#' the filtered data frame contains 100 overlaps, the 101st overlap
#' corresponds to a score of 0).
#'
#' @param overlapDF A filtered overlap data frame.
#' @param osMethod The method used to compute overlap scores.
#' Options are 'log' and 'minmax'.
#' @param firstOutRawRank The raw rank of the highest-ranked
#' overlap that was not recorded among top overlaps. Ignored
#' if osMethod is set to 'log'.
#'
#' @return A data frame with ranked overlaps.
#'
#' @export
#'
scoreOverlaps <- function(overlapDF,
                          osMethod = 'log',
                          firstOutRawRank = NULL){
  message(nrow(overlapDF), ' overlap', rep('s', nrow(overlapDF) != 1),
          ' will be used in the calculation of CSOA scores.')
  if (nrow(overlapDF) == 1){
    overlapDF$score <- 1
    return(overlapDF)
  }
  if(!osMethod %in% c('log', 'minmax'))
    stop('Unrecognized overlap scoring method.',
         'See ?CSOA::scoreOverlaps for the accepted methods.')

  if (osMethod == 'log'){
    rankVals <- unique(overlapDF$rank)
    logVals <- log(seq(exp(1), 1, length.out =
                         length(rankVals) + 1))[seq_along(rankVals)]
    names(logVals) <- rankVals
    overlapDF$score <- logVals[as.character(overlapDF$rank)]
  }

  if (osMethod == 'minmax'){
    rawRank <- c(overlapDF$rawAggRank, firstOutRawRank)
    overlapDF$score <- 1 - vMinmax(rawRank)[seq_len(nrow(overlapDF))]
  }

  return(overlapDF)
}
