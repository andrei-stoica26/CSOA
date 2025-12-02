#' Generate CSOA scores from overlap data frame for a single gene set
#'
#' This function computes per-cell CSOA scores from overlap data frame for
#' a single gene set.
#'
#' @inheritParams generateOverlaps
#' @inheritParams processOverlaps
#' @inheritParams computePairScores
#' @param colStr Name of column where CSOA scores will be stored.
#'
#' @return A data frame with a column corresponding to the CSOA scores.
#'
#' @keywords internal
#'
scoreCellsCore <- function(geneSetExp,
                           overlapDF,
                           colStr = 'CSOA',
                           mtMethod = c('BY', 'BH'),
                           jaccardCutoff = NULL,
                           osMethod = c('log', 'minmax'),
                           pairFileName = NULL,
                           keepOverlapOrder = FALSE,
                           ...){

    overlapDF <- processOverlaps(overlapDF, mtMethod,
                                 jaccardCutoff, osMethod, ...)
    if(!nrow(overlapDF)){
        warning('No significant overlaps were identified.',
                ' All cells will get a score of 0.')
        scoreDF <- data.frame(rep(0, dim(geneSetExp)[2]))
        colnames(scoreDF) <- colStr
        rownames(scoreDF) <- colnames(geneSetExp)
        return(scoreDF)
    }
    genes <- overlapGenes(overlapDF)
    message('Normalizing expression matrix by rows...')
    normExp <- kerntools::minmax(geneSetExp[genes, ], rows=TRUE)

    pcPairScores <- computePCPairScores(overlapDF, normExp)
    if(!is.null(pairFileName))
        pairScores <- computePairScores(overlapDF, pcPairScores,
                                        pairFileName, keepOverlapOrder)

    message('Computing per-cell gene signature scores...')
    scores <- colSums(pcPairScores)
    scores <- vMinmax(scores)
    scoreDF <- data.frame(setNames(list(scores), colStr))
    rownames(scoreDF) <- colnames(pcPairScores)
    return(scoreDF)
}

#' Generate CSOA scores from overlap data frame and list of pairs
#'
#' This function scores an overlap data frame using its associated list of
#' pairs. The overlap data frame is split based on the overlaps corresponding
#' to each gene set and scored, and the output is rejoined as a data frame.
#'
#' @details This function calls \code{scoreCells} to score each gene set
#' data frame split from the full overlap data frame.
#'
#' @inheritParams scoreCellsCore
#' @param setPairs A list of overlaps corresponding to each input gene set.
#' @param geneSetNames Character vector of names of gene sets.
#' @param pairFileTemplate Character object used in the naming of the files
#' where the pair data frames will be saved. Default is \code{NULL} (the pair
#' data frames will not be saved).
#' @param keepOverlapOrder Keep the rank-based order of overlaps in the
#' pair score file, as opposed to changing it to a pair score-based order.
#' Ignored if pairFileTemplate is \code{NULL}.
#'
#' @return A data frame whose columns correspond to the CSOA scores of the
#' input gene sets.
#'
#' @examples
#' mat <- matrix(0, 500, 300)
#' rownames(mat) <- paste0('G', seq(500))
#' colnames(mat) <- paste0('C', seq(300))
#' mat[sample(8000)] <- runif(8000, max=13)
#' genes <- paste0('G', seq(200))
#' mat[genes, 20:50] <- matrix(runif(200 * 31, min=14, max=15),
#' nrow=200, ncol=31)
#' geneSet1 <- paste0('G', seq(1, 150))
#' geneSet2 <- paste0('G', seq(50, 200))
#' geneSets <- list(geneSet1, geneSet2)
#' geneSets <- lapply(geneSets, sort)
#' setPairs <- lapply(geneSets, getPairs)
#' pairs <- Reduce(union, setPairs)
#' genes <- union(geneSet1, geneSet2)
#' mat <- mat[genes, ]
#' overlapDF <- generateOverlaps(mat, pairs=pairs)
#' scoreDF <- scoreCells(mat, overlapDF, setPairs, c('set1', 'set2'))
#' head(scoreDF)
#'
#' @export
#'
scoreCells <- function(geneSetExp,
                       overlapDF,
                       setPairs,
                       geneSetNames,
                       mtMethod = c('BY', 'BH'),
                       jaccardCutoff = NULL,
                       osMethod = c('log', 'minmax'),
                       pairFileTemplate = NULL,
                       keepOverlapOrder = FALSE,
                       ...){

    mtMethod <- match.arg(mtMethod, c('BY', 'BH'))
    osMethod <- match.arg(osMethod, c('log', 'minmax'))

    message('Processing overlaps...')

    if(!is.null(pairFileTemplate))
        pairFileName <- paste0(pairFileTemplate, geneSetNames) else
            pairFileName <- NULL

    scoreDFList <- lapply(seq_along(setPairs), function(i) {
        setOverlapDF <- overlapSlice(overlapDF, setPairs[[i]])
        scoreDF <- scoreCellsCore(geneSetExp, setOverlapDF, geneSetNames[i],
                                  mtMethod, jaccardCutoff, osMethod,
                                  pairFileName[i], keepOverlapOrder, ...)
        return(scoreDF)
    })

    allScoresDF <- Reduce(cbind, scoreDFList)
    return(allScoresDF)
}
