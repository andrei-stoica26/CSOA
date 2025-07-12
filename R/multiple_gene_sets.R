#' Generate CSOA scores from overlap data frame for multiple gene sets
#'
#' This function scores an overlap data frame generated using multiple
#' gene sets.
#' The overlap data frame is split based on the overlaps corresponding to
#' each gene set and scored, and the output is rejoined as a data frame.
#'
#' @details This function calls \code{scoreCells} to score each gene set
#' data frame split from the full overlap data frame.
#'
#' @inheritParams scoreCells
#' @param setPairs A list of overlaps corresponding to each input gene set.
#' @param geneSetNames Character vector of names of gene sets.
#' @param pairFileTemplate Character object used in the naming of the files
#' where the pair data frames will be saved. Default is NULL (the pair
#' data frames will not be saved).
#' @param keepOverlapOrder Keep the rank-based order of overlaps in the
#' pair score file, as opposed to changing it to a pair score-based order.
#' Ignored if pairFileTemplate is NULL.
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
#' mat[genes, 20:50] <- matrix(runif(200 * 31, min = 14, max = 15), nrow = 200, ncol = 31)
#' geneSet1 <- paste0('G', seq(1, 150))
#' geneSet2 <- paste0('G', seq(50, 200))
#' geneSets <- list(geneSet1, geneSet2)
#' geneSets <- lapply(geneSets, sort)
#' setPairs <- lapply(geneSets, getPairs)
#' pairs <- Reduce(union, setPairs)
#' genes <- union(geneSet1, geneSet2)
#' mat <- mat[genes, ]
#' overlapDF <- generateOverlaps(mat, pairs = pairs)
#' scoreDF <- scoreCellsMultiple(mat, overlapDF, setPairs, c('set1', 'set2'))
#' head(scoreDF)
#'
#' @export
#'
scoreCellsMultiple <- function(geneSetExp,
                               overlapDF,
                               setPairs,
                               geneSetNames,
                               pvalThr = 0.05,
                               saveCutoffPlot = FALSE,
                               jaccardCutoff = NULL,
                               osMethod = 'log',
                               pairFileTemplate = NULL,
                               keepOverlapOrder = FALSE){
    if(!is.null(pairFileTemplate))
        pairFileName <- paste0(pairFileTemplate, geneSetNames) else
            pairFileName <- NULL
    scoreDFList <- lapply(seq_along(setPairs), function(i) {
        setOverlapDF <- overlapSlice(overlapDF, setPairs[[i]])
        scoreDF <- scoreCells(geneSetExp, setOverlapDF, geneSetNames[i],
                              pvalThr, saveCutoffPlot, jaccardCutoff,
                              osMethod, pairFileName[i], keepOverlapOrder)
        return(scoreDF)
    })
    allScoresDF <- Reduce(cbind, scoreDFList)
    return(allScoresDF)
}

#' Run the CSOA pipeline for multiple gene sets
#'
#' This function generates cell set overlaps for input gene sets
#' based on percentiles of gene expression, computes the significance
#' of these overlaps, ranks, filters and scores the overlaps, and builds a
#' per-cell score by summing the products of overlap scores and the
#' min-max-normalized expression of the corresponding pairs of genes.
#'
#' #' @details Wrapper around \code{expMat}, \code{generateOverlaps},
#' \code{scoreCellsMultiple} and \code{attachCellScores}.
#'
#' @inheritParams runCSOA
#' @param geneSets List of character vectors
#' @inheritParams scoreCellsMultiple
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
#' df <- runCSOAMultiple(mat, list(geneSet1, geneSet2), c('set1', 'set2'))
#' head(df)
#'
#' @export
#'
runCSOAMultiple <- function(scObj,
                            geneSets,
                            geneSetNames,
                            percentile = 90,
                            overlapFileName = NULL,
                            pvalThr = 0.05,
                            saveCutoffPlot = FALSE,
                            jaccardCutoff = NULL,
                            osMethod = 'log',
                            pairFileTemplate = NULL,
                            keepOverlapOrder = FALSE){
    geneSets <- lapply(geneSets, sort)
    setPairs <- lapply(geneSets, getPairs)
    pairs <- Reduce(union, setPairs)
    genes <- Reduce(union, geneSets)
    geneSetExp <- expMat(scObj, genes)
    overlapDF <- generateOverlaps(geneSetExp, percentile,
                                  pairs, overlapFileName)
    scoreDF <- scoreCellsMultiple(geneSetExp, overlapDF,
                                  setPairs, geneSetNames,
                                  pvalThr, saveCutoffPlot,
                                  jaccardCutoff, osMethod,
                                  pairFileTemplate, keepOverlapOrder)
    return(attachCellScores(scObj, scoreDF))
}
