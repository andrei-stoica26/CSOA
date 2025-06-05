#' @importFrom stats quantile
#'
NULL


#' Generates cell expressing input genes at an input percentile
#'
#' This function first finds, for each input gene, the cells with a non-zero
#' expression of the gene. Subsequently, each cell set is filtered as to contain
#' only the cells showcasing the input percentile of the expression of the gene
#' defining the set.
#'
#' @inheritParams expMat
#' @param genes Vector of genes. Must include at least two genes
#' @param percentile A non-negative number under 100
#'
#' @return A named list of character vectors of length equaling the number of
#' input genes, storing, for each gene, the cells showing the input percentile in
#' terms of their expression of the gene
#'
#' @export
#'
percentileSets <- function(scObj, genes, percentile=90){
  if (!min(is(genes)[1:2] == c('character', 'vector')) | length(genes) < 2)
    stop('genes must be a character vector of length >= 2')
  if (!is.numeric(percentile) | length(percentile) > 2 | percentile < 0 | percentile >= 100)
    stop('percentile must be a non-negative number lower than 100')
  expression <- expMat(scObj)
  fraction <- percentile / 100
  message('Computing percentile sets...')
  expList <- lapply(genes, function(x){
    geneExp <- expression[x, ]
    geneExp <- geneExp[geneExp > 0]
    thresh <- as.numeric(quantile(geneExp, fraction))
    return(names(geneExp[geneExp > thresh]))
  })
  names(expList) <- genes
  return(expList)
}
