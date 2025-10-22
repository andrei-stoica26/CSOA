#' Run the CSOA pipeline
#'
#' This function generates cell set overlaps for input gene sets
#' based on percentiles of gene expression, computes the significance
#' of these overlaps, ranks, filters and scores the overlaps, and builds a
#' per-cell score by summing the products of overlap scores and the
#' min-max-normalized expression of the corresponding pairs of genes.
#'
#' @details Wrapper around \code{expMat}, \code{generateOverlaps},
#' \code{scoreCells} and \code{attachCellScores}.
#'
#' @inheritParams expMat
#' @inheritParams generateOverlaps
#' @param geneSets Named list of character vectors of which each must contain
#' at least two genes.
#' @inheritParams scoreCells
#' @param pairFileTemplate Character object used in the naming of the files
#' where the pair data frames will be saved. Default is \code{NULL} (the pair
#' data frames will not be saved).
#' @param keepOverlapOrder Keep the rank-based order of overlaps in the
#' pair score file, as opposed to changing it to a pair score-based order.
#' Ignored if \code{pairFileTemplate} is \code{NULL}.
#'
#' @return An object of the same class as scObj with per-gene-set CSOA scores
#' assigned for each cell.
#'
#' @examples
#' mat <- matrix(0, 500, 300)
#' rownames(mat) <- paste0('G', seq(500))
#' colnames(mat) <- paste0('C', seq(300))
#' mat[sample(8000)] <- runif(8000, max=13)
#' genes <- paste0('G', seq(200))
#' mat[genes, 20:50] <- matrix(runif(200 * 31, min = 14, max = 15),
#' nrow = 200, ncol = 31)
#' geneSet1 <- paste0('G', seq(1, 150))
#' geneSet2 <- paste0('G', seq(50, 200))
#' df <- runCSOA(mat, list(a = geneSet1, b = geneSet2))
#' head(df)
#'
#' @export
#'
runCSOA <- function(scObj,
                    geneSets,
                    percentile = 90,
                    mtMethod = c('by', 'bh', 'bf'),
                    jaccardCutoff = NULL,
                    osMethod = c('log', 'minmax'),
                    overlapFileName = NULL,
                    pairFileTemplate = NULL,
                    keepOverlapOrder = FALSE,
                    ...){

    mtMethod <- match.arg(mtMethod, c('by', 'bh', 'bf'))
    osMethod <- match.arg(osMethod, c('log', 'minmax'))
    if (is.null(names(geneSets)))
        stop('The gene sets must have names.')
    geneSets <- lapply(geneSets, sort)
    setPairs <- lapply(geneSets, getPairs)
    pairs <- Reduce(union, setPairs)
    genes <- Reduce(union, geneSets)
    geneSetExp <- expMat(scObj, genes)
    overlapDF <- generateOverlaps(geneSetExp, percentile, pairs,
                                  overlapFileName)
    scoreDF <- scoreCells(geneSetExp, overlapDF, setPairs, names(geneSets),
                          mtMethod, jaccardCutoff, osMethod, pairFileTemplate,
                          keepOverlapOrder, ...)
    return(attachCellScores(scObj, scoreDF))
}

