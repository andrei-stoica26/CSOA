#' @importFrom stats p.adjust phyper quantile runif setNames
#'
NULL

#' Generates cell expressing input genes at an input percentile
#'
#' This function constructs, for each gene in the expression matrix, a set of
#' cells expressing the gene at or above the input percentile.
#'
#' @param geneSetExp A gene expression non-sparse matrix with the rows
#' restricted to the genes for which cell sets will be computed.
#' @param percentile A positive number under 100.
#'
#' @return A named list of character vectors of length
#' equaling the number of input genes. Each vector stores the cells expressing
#' the gene at or above the input percentile.
#'
#' @examples
#' mat <- matrix(0, 1000, 500)
#' rownames(mat) <- paste0('G', seq(1000))
#' colnames(mat) <- paste0('C', seq(500))
#' mat[sample(length(mat), 70000)] <- sample(50, 70000, TRUE)
#' mat <- mat[paste0('G', sample(1000, 3)), ]
#' percentileSets(mat)
#'
#' @export
#'
percentileSets <- function(geneSetExp, percentile=90){
    if (!is.numeric(percentile) | length(percentile) > 2 |
        percentile < 0 | percentile >= 100)
        stop('percentile must be a non-negative number lower than 100.')
    if (is.null(rownames(geneSetExp)))
        stop('geneSetExp has no row names.')
    if (is.null(colnames(geneSetExp)))
        stop('geneSetExp has no column names.')
    genes <- rownames(geneSetExp)
    fraction <- percentile / 100
    message('Computing percentile sets...')
    expList <- lapply(genes, function(x){
        geneExp <- geneSetExp[x, ]
        geneExp <- geneExp[geneExp > 0]
        thresh <- as.numeric(quantile(geneExp, fraction))
        return(names(geneExp[geneExp > thresh]))
    })
    names(expList) <- genes
    expList <- expList[vapply(expList, length, numeric(1)) > 0]
    if (!length(expList))
        warning('No cell sets can be constructed at the indicated percentile',
                ' for the input genes.')
    if (length(expList) < length(genes))
        warning(length(genes) - length(expList), ' gene(s) had no',
                ' top cells at the indicated percentile.',
                ' These are now excluded from the gene signature.')
    return(expList)
}
