#' @importClassesFrom Seurat Seurat
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment assay
#' @importFrom stats setNames
NULL

#' Generate overlaps of cell sets for input genes
#'
#' This function constructs, for each gene in the expression matrix, a set of
#' cells expressing the gene at or above the input percentile.
#' Subsequently, overlaps of pairs of the constructed cell sets are assessed
#' for statistical significance.
#'
#' @details Wrapper around \code{percentileSets} and \code{cellSetsOverlaps}.
#' @inheritParams percentileSets
#' @inheritParams cellSetsOverlaps
#'
#' @return A data frame listing statistics for all cell set overlaps
#'
#' @examples
#' mat <- matrix(0, 2000, 500)
#' rownames(mat) <- paste0('G', seq(2000))
#' colnames(mat) <- paste0('C', seq(500))
#' mat[sample(length(mat), 270000)] <- sample(50, 270000, TRUE)
#' mat <- mat[paste0('G', sample(2000, 5)), ]
#' generateOverlaps(mat)
#'
#' @export
#'
generateOverlaps <- function(geneSetExp, percentile = 90, pairs = NULL,
                             overlapFileName = NULL){
    cellSets <- percentileSets(geneSetExp, percentile)
    if(!length(cellSets))
        return(data.frame())
    overlapDF <- cellSetsOverlaps(cellSets, dim(geneSetExp)[2], pairs,
                                  overlapFileName)
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
#' @inheritParams byCorrectDF
#' @inheritParams rankOverlaps
#' @inheritParams prepareFiltering
#' @inheritParams filterOverlaps
#' @param jaccardCutoff A cutoff used in the filtering of edges with low
#' Jaccard scores. If \code{NULL} (as default), no filtering of such edges
#' will be performed.
#' @inheritParams scoreOverlaps
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
                            pvalThr = 0.05,
                            saveCutoffPlot = FALSE,
                            jaccardCutoff = NULL,
                            osMethod = 'log'){

    if (nrow(overlapDF) > 1)
        overlapDF <- byCorrectDF(overlapDF, pvalThr) else
            overlapDF$pval_adj <- overlapDF$pval

    if (!nrow(overlapDF))
        return(overlapDF)

    overlapDF <- rankOverlaps(overlapDF)
    firstOutRawRank <- prepareFiltering(overlapDF, saveCutoffPlot)
    overlapDF <- filterOverlaps(overlapDF, firstOutRawRank)
    if (!is.null(jaccardCutoff)){
        overlapDF <- breakWeakTies(overlapDF, jaccardCutoff)
        firstOutRawRank <- NULL
    }

    overlapDF <- scoreOverlaps(overlapDF, osMethod, firstOutRawRank)
    return(overlapDF)
}

#' Assign a per-cell gene set score to Seurat object
#'
#' This function uses the scored data frame of overlaps to
#' compute a CSOA score for each cell in a Seurat object.
#'
#' @details The per-cell score is the sum of all products of overlap scores and
#' the min-max-normalized gene expression of each of the two overlap genes.
#'
#' Wrapper around \code{computePCPairScores} and \code{computePCSetScores}.
#'
#' @inheritParams computePCPairScores
#' @inheritParams computePairScores
#' @inheritParams computePCSetScores
#'
#' @return A Seurat object with a CSOA score assigned for each cell.
#'
#' @examples
#' overlapDF <- data.frame(gene1 = paste0('G', c(1, 2, 7, 8)),
#' gene2 = paste0('G', c(3, 7, 1, 2)),
#' score = c(0.22, 0.98, 1, 0.76))
#' normExp <- matrix(0, 5, 26)
#' normExp[sample(length(normExp), 50)] <- runif(50)
#' normExp[1, 2] <- 1
#' rownames(normExp) <- union(overlapDF$gene1, overlapDF$gene2)
#' computeCellScores(overlapDF, normExp)
#'
#' @export
#'
computeCellScores <- function(overlapDF,
                              normExp,
                              colStr = 'CSOA',
                              pairFileName = NULL,
                              keepOverlapOrder = FALSE){
    pcPairScores <- computePCPairScores(overlapDF, normExp)
    if(!is.null(pairFileName))
        pairScores <- computePairScores(overlapDF, pcPairScores,
                                        pairFileName, keepOverlapOrder)
    scoreDF <- computePCSetScores(pcPairScores, colStr)
    return(scoreDF)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for CSOA-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @param scoreDF Dataframe of CSOA scores
#'
#' @rdname attachCellScores
#' @export
#'
attachCellScores.default <- function(scObj, scoreDF, ...)
    stop('Unrecognized input type: scObj must be a Seurat object with a',
       ' data assay, a SingleCellExperiment with a logcounts assay',
       ' a matrix or a dgCMatrix.')

#' @rdname attachCellScores

#' @return A Seurat object with CSOA scores added to metadata.
#'
#' @export
#'
attachCellScores.Seurat <- function(scObj, scoreDF, ...){
    for (colName in colnames(scoreDF))
        if (colName %in% colnames(scObj[[]]))
            scObj[[]][[colName]] <- c()
        scObj[[]] <- cbind(scObj[[]], scoreDF)
    return(scObj)
}

#' @rdname attachCellScores
#'
#' @return A SingleCellExperiment object with CSOA scores added to
#' \code{colData}.
#'
#' @export
#'
attachCellScores.SingleCellExperiment <- function(scObj, scoreDF, ...){
    for (colName in colnames(scoreDF))
        if (colName %in% colnames(colData(scObj)))
            colData(scObj)[[colName]] <- c()
    colData(scObj) <- cbind(colData(scObj), scoreDF)
    return(scObj)
}

#' @rdname attachCellScores
#'
#' @return A list containing the expression matrix and the CSOA scores data
#' frame.
#'
#' @export
#'
attachCellScores.matrix <- function(scObj, scoreDF, ...)
    return(list(object = scObj, scores = scoreDF))

#' @rdname attachCellScores
#'
#' @return A list containing the expression matrix and the CSOA scores data
#' frame.
#'
#' @export
#'
attachCellScores.dgCMatrix <- function(scObj, scoreDF, ...)
        return(list(object = scObj, scores = scoreDF))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Generate CSOA scores from overlap data frame
#'
#' This function computes per-cell CSOA scores after the overlap data frame has
#' been generated.
#'
#' @details Wrapper around \code{processOverlaps} and \code{computeCellScores}.
#'
#' @inheritParams generateOverlaps
#' @inheritParams processOverlaps
#' @inheritParams computeCellScores
#'
#' @return A data frame with a column corresponding to the CSOA scores
#'
#' @examples
#' mat <- matrix(0, 500, 300)
#' rownames(mat) <- paste0('G', seq(500))
#' colnames(mat) <- paste0('C', seq(300))
#' mat[sample(8000)] <- runif(8000, max=13)
#' genes <- paste0('G', seq(200))
#' mat[genes, 20:50] <- matrix(runif(200 * 31, min = 14, max = 15), nrow = 200, ncol = 31)
#' genes <- paste0('G', seq(1, 200))
#' mat <- mat[genes, ]
#' overlapDF <- generateOverlaps(mat)
#' scoreDF <- scoreCells(mat, overlapDF)
#' head(scoreDF)
#'
#' @export
#'
scoreCells <- function(geneSetExp,
                       overlapDF,
                       colStr = 'CSOA',
                       pvalThr = 0.05,
                       saveCutoffPlot = FALSE,
                       jaccardCutoff = NULL,
                       osMethod = 'log',
                       pairFileName = NULL,
                       keepOverlapOrder = FALSE){
    overlapDF <- processOverlaps(overlapDF, pvalThr, saveCutoffPlot,
                                 jaccardCutoff, osMethod)
    if(!nrow(overlapDF)){
        warning('No significant overlaps were identified.',
                ' All cells will get a score of 0.')
        scoreDF <- data.frame(setNames(rep(0, dim(geneSetExp)[2]),
                                       colStr))
        rownames(scoreDF) <- colnames(geneSetExp)
        return(scoreDF)
    }

    message('Normalizing expression matrix by rows...')
    genes <- overlapGenes(overlapDF)
    normExp <- kerntools::minmax(geneSetExp[genes, ], rows=TRUE)
    scoreDF <- computeCellScores(overlapDF, normExp, colStr,
                                 pairFileName, keepOverlapOrder)
    return(scoreDF)
}

#' Run the CSOA pipeline
#'
#' This function generates cell set overlaps for an input gene set
#' based on percentiles of gene expression, computes the significance
#' of these overlaps, ranks, filters and scores the overlaps, and builds a
#' per-cell score by summing the products of overlap scores and the
#' min-max-normalized expression of the corresponding pairs of genes.
#'
#' @details Wrapper around \code{expMat}, \code{generateOverlaps},
#' \code{scoreCells} and \code{attachCellScores}.
#'
#' @inheritParams expMat
#' @param genes Vector of genes. Must include at least two genes.
#' @inheritParams generateOverlaps
#' @inheritParams scoreCells
#
#' @return An object of the same class as scObj with a CSOA score assigned
#' for each cell.
#'
#' @examples
#' mat <- matrix(0, 500, 300)
#' rownames(mat) <- paste0('G', seq(500))
#' colnames(mat) <- paste0('C', seq(300))
#' mat[sample(8000)] <- runif(8000, max=15)
#' genes <- paste0('G', seq(200))
#' mat[genes, 20:50] <- matrix(runif(200 * 31, min = 14, max = 15), nrow = 200, ncol = 31)
#' df <- runCSOA(mat, genes)
#' head(df)
#'
#' @export
#'
runCSOA <- function(scObj, genes, colStr='CSOA', percentile = 90,
                    overlapFileName = NULL, pvalThr = 0.05,
                    saveCutoffPlot = FALSE, jaccardCutoff = NULL,
                    osMethod = 'log', pairFileName = NULL,
                    keepOverlapOrder = FALSE){
    if (!min(is(genes)[c(1, 2)] == c('character', 'vector')) | length(genes) < 2)
        stop('genes must be a character vector of length >= 2.')
    geneSetExp <- expMat(scObj, genes)
    overlapDF <- generateOverlaps(geneSetExp, percentile, pairs=NULL,
                                  overlapFileName)
    scoreDF <- scoreCells(geneSetExp, overlapDF, colStr, pvalThr,
                          saveCutoffPlot, jaccardCutoff, osMethod,
                          pairFileName, keepOverlapOrder)
    return(attachCellScores(scObj, scoreDF))
}
