#' @importFrom qs qsave
#' @importFrom stringr str_c
#' @importFrom stats phyper
#'
NULL

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
#' @param pairs Pairs of cell sets to be assessed. If NULL (as default), all
#' pairs will be assessed.
#' @param overlapFileName The name of the file where the overlap data frame
#' will be saved. Default is \code{NULL} (the overlap data frame will not
#' be saved).
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
cellSetsOverlaps <- function(cellSets, nCells, pairs = NULL,
                             overlapFileName = NULL){
  message('Assessing gene overlaps...')
  genes <- names(cellSets)
  if(!length(genes))
    stop('The cell sets must be named')
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
