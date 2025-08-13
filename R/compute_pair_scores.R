#' Compute per-cell gene pair scores
#'
#' This function scores each gene pair corresponding to a top overlap in each
#' cell.
#'
#' @details The score is calculated by multiplying the overlap score with the
#' min-max-normalized expression of the two corresponding genes.
#'
#' @param overlapDF An overlap data frame.
#' @param normExp A min-max normalized expression matrix of the genes involved in
#' top overlaps.
#'
#' @return A data frame with gene pairs as rows and cells as columns.
#'
#' @keywords internal
#'
computePCPairScores <- function(overlapDF, normExp){
    if(length(setdiff(c('gene1', 'gene2', 'score'), colnames(overlapDF))))
        stop('Columns gene1, gene2 and score must exist in overlapDF.')
    if(max(overlapDF$score) != 1 | min(overlapDF$score) <= 0)
        stop('The maximum of the score column in overlapDF must be 1',
             ' and its minimum must be positive.')
    if(!is.numeric(normExp) | !is.matrix(normExp))
        stop('normExp must be a numeric matrix.')
    if(max(normExp) != 1 | min(normExp) != 0)
        stop('The maximum value of normExp must be 1',
             'and its minimum must be 0.')

    message('Computing per-cell scores for gene pairs...')
    gene1 <- overlapDF$gene1
    gene2 <- overlapDF$gene2
    scores <- overlapDF$score

    pcPairScores <- normExp[gene1, , drop=FALSE] * normExp[gene2, , drop=FALSE]
    pcPairScores <- pcPairScores * scores
    pcPairScores <- as.data.frame(pcPairScores)
    rownames(pcPairScores) <- paste0(gene1, "_", gene2)
    return(pcPairScores)
}

#' Compute aggregate gene pair scores
#'
#' This function assesses the relative contribution of each gene pair to the
#' CSOA score
#'
#' @inheritParams computePCPairScores
#' @param pcPairScores A date frame of pair scores in each cell for each pair
#' in the overlap data frame.
#' @param pairFileName The name of the file where the pair data frame
#' will be saved.
#' @param keepOverlapOrder Whether to keep the rank-based order of overlaps
#' in the pair score file, as opposed to changing it to a
#' pair score-based order.
#'
#' @return A data frame with overlap and pair scores and ranks.
#'
#' @keywords internal
#'
computePairScores <- function(overlapDF, pcPairScores, pairFileName = 'pairs',
                              keepOverlapOrder = FALSE){

    df <- overlapDF[, c('gene1', 'gene2', 'score', 'rank')]
    colnames(df)[c(3, 4)] <- paste0('overlap', c('Score', 'Rank'))

    pairTotalScores <- rowSums(pcPairScores)
    totalScore <- sum(pairTotalScores)

    df$pairScore <- pairTotalScores / totalScore * 100
    df <- df[order(df$pairScore, decreasing=TRUE), ]

    df$pairRank <- rankFun(-df$pairScore)
    df$revCumsum <- spatstat.utils::revcumsum(df$pairScore)

    if (keepOverlapOrder)
        df <- df[order(df$overlapScore, decreasing=TRUE), ]

    pairFile <- paste0(pairFileName, '.qs')
    message('Saving pair scores file: ', pairFile, '...')
    qsave(df, pairFile)

    return(df)
}
