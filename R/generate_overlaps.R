#' Compute the pairwise overlaps of a list of cell sets
#'
#' This function performs the pairwise overlaps of a list of cell sets.
#'
#' @param cellSets A named list of two character arrays.
#' @param pairs Pairs.
#' @param nCells An integer.
#'
#' @return A vector comprising the names of genes, the cell counts,
#' the recorded and expected shared cells, the recorded-over-expected ratio,
#' and the hypergeometric p-value.
#'
#' @noRd
#'
pairOverlap <- function(cellSets, pairs, nCells){
    gene1 <- vapply(pairs, `[[`, character(1), 1)
    gene2 <- vapply(pairs, `[[`, character(1), 2)

    set1 <- cellSets[gene1]
    set2 <- cellSets[gene2]

    xCount <- lengths(set1)
    yCount <- lengths(set2)

    recorded <- mapply(function(a, b) length(intersect(a, b)), set1, set2)
    expected <- xCount * yCount / nCells
    ratio <- recorded / expected
    pval <- phyper(recorded - 1, xCount, nCells - xCount, yCount,
                   lower.tail = FALSE)

    df <- data.frame(gene1 = gene1, gene2 = gene2, ncells1 = xCount,
                     ncells2 = yCount, shared_cells = recorded,
                     exp_shared_cells = expected, ratio = ratio,
                     pval = pval)
    return(df)
}

#' Calculates the significance of overlaps of pairs of cells sets
#'
#' This function computes the statistical significance of overlaps of pairs of
#' cell sets.
#'
#' @param cellSets A list of character arrays.
#' @param nCells The total number of cells in the Seurat object.
#' @param pairs Pairs of cell sets to be assessed. If \code{NULL} (as default),
#' all pairs will be assessed.
#' @param overlapFileName The name of the file where the overlap data frame
#' will be saved. This option can be used to save time when performing
#' exploratory analyses such as trying different \code{jaccardCutoff} parameters
#' in \code{breakWeakTies}. Default is \code{NULL} (the overlap data frame will
#' not be saved).
#'
#' @return A data frame listing statistics for all cell set overlaps: cell set
#' sizes, recorded and expected shared cells, the recorded-over-expected ratio
#' and the hypergeometric p-value.
#'
#' @examples
#' cellSets <- list(G1 = c('A', 'H', 'J'),
#' G2 = c('B', 'D', 'E', 'F', 'J'),
#' G3 = c('C', 'I', 'L'))
#' cellSetsOverlaps(cellSets, 40)
#'
#' @export
#'
cellSetsOverlaps <- function(cellSets,
                             nCells,
                             pairs = NULL,
                             overlapFileName = NULL){
    message('Generating cell set overlaps...')
    genes <- names(cellSets)
    if(!length(genes))
        stop('The cell sets must have names.')
    if(is.null(pairs))
        pairs <- getPairs(genes)
    df <- pairOverlap(cellSets, pairs, nCells)
    if (!is.null(overlapFileName)){
        overlapFile <- paste0(overlapFileName, '.qs')
        message('Saving overlap file: ', overlapFile, '...')
        qsave(df, overlapFile)
    }
    return(df)
}

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
generateOverlaps <- function(geneSetExp,
                             percentile = 90,
                             pairs = NULL,
                             overlapFileName = NULL){
    cellSets <- percentileSets(geneSetExp, percentile)
    if(!length(cellSets))
        return(data.frame())
    overlapDF <- cellSetsOverlaps(cellSets,
                                  dim(geneSetExp)[2],
                                  pairs,
                                  overlapFileName)
    return(overlapDF)
}
