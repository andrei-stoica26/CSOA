#' @importFrom spatstat.utils revcumsum
#'
NULL

#' Find the connectivity of each gene
#'
#' This function finds the connectivity of each gene
#' from an overlap data frame ranked using p-values
#' and recorded-over-expected ratios.
#'
#' @param overlapDF Overlap data frame.
#' @param asRanks Whether to replace connectivity
#' scores by ranks.
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
#' This function ranks cell set overlaps.
#'
#' @details Overlap ranks are calculated as follows:
#'
#' 1. Two preliminary overlap ranks are computed based on:
#' 1) adjusted p-value;
#' 2) recorded-over-expected size ratio.
#'
#' 2. The preliminary ranks are replaced by \strong{connectivity-based ranks}
#' (\code{pvalRank}, \code{ratioRank}). For each edge corresponding to an
#' overlap, the p-value rank is replaced by the minimum p-value rank of the
#' vertices neighboring the edge—genes with a significant cell set overlap
#' with any of the two overlap genes. The ratio rank is replaced similarly.
#'
#' 3. A raw aggregate rank (\code{rawAggRank}) is computed as the average of
#' \code{pvalRank} and \code{ratioRank}.
#'
#' 4. A final aggregate rank (\code{rank}) defined as the order imposed by the
#' raw aggregate rank is constructed.
#'
#'
#' @param overlapDF An overlap data frame.
#'
#' @return A data frame with ranked overlaps.
#'
#' @noRd
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
#' This function finds the rank cutoff for overlap filtering.
#'
#' @details The rank cutoff is the average of the minimum and maximum
#' highest-frequency rank.
#'
#' @param freqDF A data frame of overlap rank frequencies.
#'
#' @return Rank cutoff.
#'
#' @noRd
#'
findRankCutoff <- function(freqDF){
    if(length(setdiff(c('rank', 'n'), colnames(freqDF))))
        stop('Columns `rank` and `n` must exist in the data frame.')
    freqSub <- freqDF[freqDF$n == max(freqDF$n), , drop=FALSE]
    rankCutoff <- mean(c(max(freqSub$rank), min(freqSub$rank)))
    return(rankCutoff)
}

#' Find the raw aggregate rank of the highest non-top overlap
#'
#' This function finds the raw aggregate rank of the
#' highest non-top overlap.
#'
#' @details This function calls \code{findRankCutoff} to determine the rank
#' cutoff and uses this cutoff to find the best-ranking overlap excluded by it.
#'
#' @param overlapDF A ranked overlap data frame.
#'
#' @return A numeric value
#'
#' @noRd
#'
prepareFiltering <- function(overlapDF){
    if (!nrow(overlapDF))
        return(NULL)

    freqDF <- dplyr::count(overlapDF, rank)
    rankCutoff <- findRankCutoff(freqDF)
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
#' has been ranked.
#'
#' @details If \code{firstOutRawRank} is \code{NULL}, the data frame will be
#' returned unchanged.
#'
#' @param overlapDF A ranked overlap data frame.
#' @param firstOutRawRank The raw aggregate rank of the
#' first overlap that will be excluded.
#'
#' @return A filtered overlap data frame.
#'
#' @noRd
#'
filterOverlaps <- function(overlapDF, firstOutRawRank = NULL){
    if(!'rawAggRank' %in% colnames(overlapDF))
        stop('Column `rawAggRank` must exist in the data frame.')
    if(is.null(firstOutRawRank) | !nrow(overlapDF))
        return(overlapDF)
    return(overlapDF[overlapDF$rawAggRank < firstOutRawRank, , drop=FALSE])
}

#' Score cell set overlaps
#'
#' This function scores cell set overlaps based on their ranks.
#'
#' @details Both the log and the minmax method assign higher scores to better
#' ranks and equal scores to equal ranks. All scores are in (0, 1].
#'
#' The log method maps each of the \eqn{n} unique overlap raw ranks sorted
#' increasingly to \eqn{\log\left(e - \frac{(e - 1) k}{n}\right)},
#'  \eqn{k \in \{0, \ldots n - 1\}}.
#'
#' The minmax method performs min-max normalization on the set consisting of
#' the unique overlap raw ranks and a minimum raw rank provided by
#' \code{prepareFiltering} in the standard pipeline. The latter is introduced
#' to guarantee that all top overlaps will get a positive score.
#'
#' @param overlapDF A ranked and filtered overlap data frame.
#' @param osMethod Method used to compute overlap scores.
#' Options are "log" and "minmax".
#' @param firstOutRawRank Raw rank of the highest-ranked non-top overlap.
#' Ignored if \code{osMethod} is "log".
#'
#' @return An overlap data frame with overlap scores.
#'
#' @noRd
#'
scoreOverlaps <- function(overlapDF,
                          osMethod = c('log', 'minmax'),
                          firstOutRawRank = NULL){
    message(nrow(overlapDF), ' overlap', rep('s', nrow(overlapDF) != 1),
            ' will be used in the calculation of CSOA scores.')
    if (nrow(overlapDF) == 1){
        overlapDF$score <- 1
        return(overlapDF)
    }

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

#' Process data frame of overlaps of cell sets
#'
#' This function filters, ranks and scores previously generated
#' overlaps of cell sets.
#'
#' @details Wrapper around \code{byCorrectDF}, \code{rankOverlaps},
#' \code{prepareFiltering}, \code{filterOverlaps} and \code{scoreOverlaps}.
#'
#' If \code{jaccardCutoff} is not \code{NULL}, it also calls
#' \code{breakWeakTies} between \code{filterOverlaps} and \code{scoreOverlaps}.
#'
#' @param overlapDF Overlap data frame.
#' @param mtMethod Multiple testing correction method. Options are
#' Bonferroni ('bf'), Benjamini-Hochberg('bh'), and the default
#' Benjamini-Yekutieli ('by').
#' @param jaccardCutoff A cutoff used in the filtering of edges with low
#' Jaccard scores. If \code{NULL} (as default), no filtering of such edges
#' will be performed.
#' @param osMethod Method used to compute overlap scores.
#' Options are "log" and "minmax".
#' @param ... Additional arguments passed to \code{mtCorrectDF}.
#'
#' @return A data frame consisting of filtered, ranked and scored cell sets
#' overlaps
#'
#' @examples
#' overlapDF <- data.frame(gene1=paste0('G',
#' c(1, 3, 7, 6, 8, 2, 4, 3, 4, 5)),
#' gene2=paste0('G',
#' c(2, 7, 2, 5, 4, 5, 1, 2, 2, 8)),
#' ratio=runif(10, 2, 10),
#' pval=runif(10, 0, 1e-10))
#' processOverlaps(overlapDF)
#'
#' @export
#'
processOverlaps <- function(overlapDF,
                            mtMethod = c('by', 'bh', 'bf'),
                            jaccardCutoff = NULL,
                            osMethod = c('log', 'minmax'),
                            ...){

    mtMethod <- match.arg(mtMethod, c('by', 'bh', 'bf'))
    osMethod <- match.arg(osMethod, c('log', 'minmax'))

    if (nrow(overlapDF) > 1)
        overlapDF <- mtCorrectDF(overlapDF, mtMethod, ...) else
            overlapDF$pval_adj <- overlapDF$pval

    if (!nrow(overlapDF))
        return(overlapDF)

    overlapDF <- rankOverlaps(overlapDF)
    firstOutRawRank <- prepareFiltering(overlapDF)
    overlapDF <- filterOverlaps(overlapDF, firstOutRawRank)
    if (!is.null(jaccardCutoff)){
        overlapDF <- breakWeakTies(overlapDF, jaccardCutoff)
        firstOutRawRank <- NULL
    }

    overlapDF <- scoreOverlaps(overlapDF, osMethod, firstOutRawRank)
    return(overlapDF)
}
