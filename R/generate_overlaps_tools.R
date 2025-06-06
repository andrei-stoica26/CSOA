#' @importFrom data.table transpose
#' @importFrom qs qsave
#' @importFrom stringr str_c
#' @importFrom stats phyper
#'
NULL

#' Compute the overlap of two cell sets
#'
#' This function performs the overlap
#' column of p-values. Additionally, the function filters the dataframe based on
#' p-values, unless pvalThr is set as NULL.
#'
#' @param pairCellSets A named list of two character arrays
#' @param nCells An integer
#'
#' @return A vector comprising the names of genes, the cell counts, the recorded and
#' expected shared cells, the recorded-over-expected ratio, and the hypergeometric
#' p-value
#'
pairOverlap <- function(pairCellSets, nCells){
  xCount <- length(pairCellSets[[1]])
  yCount <- length(pairCellSets[[2]])
  recorded <- length(intersect(pairCellSets[[1]], pairCellSets[[2]]))
  expected <- xCount * yCount / nCells
  ratio <-  recorded / expected
  pval <- phyper(recorded - 1, xCount, nCells - xCount, yCount, F)
  res <- c(names(pairCellSets)[1], names(pairCellSets)[2], xCount, yCount, recorded, expected, ratio, pval)
  return(res)
}

#' Calculates the significance of overlaps of pairs of cells sets
#'
#' This function computes the statistical significance of overlaps of pairs of
#' cell sets
#'
#' @param cellSets A list of character arrays
#' @param nCells The total number of cells in the Seurat object
#' @param pairs Pairs of cell sets to be assessed. If NULL (as default), all
#' pairs will be assessed
#' @param overlapFileName The name of the file where the overlap data frame
#' will be saved. Default is NULL (the overlap data frame will not be saved)
#'
#' @return A data frame listing statistics for all cell set overlaps
#'
cellSetsOverlaps <- function(cellSets, nCells, pairs = NULL, overlapFileName = NULL){
  message('Assessing gene overlaps...')
  genes <- names(cellSets)
  if(is.null(pairs))
    pairs <- getPairs(genes)
  df <- lapply(pairs, function(x) pairOverlap(cellSets[x], nCells))
  df <- data.frame(Reduce(rbind, df))
  if (ncol(df) == 1)
    df <- data.table::transpose(df)
  df[, c(3:8)] <- apply(df[, c(3:8)], 2, as.numeric)
  colnames(df) <- c('gene1', 'gene2', 'ncells1', 'ncells2', 'shared_cells', 'exp_shared_cells', 'ratio', 'pval')
  rownames(df) <- 1:length(rownames(df))
  if (!is.null(overlapFileName))
    qsave(df, paste0(overlapFileName, '.qs'))
  return(df)
}

overlapPairs <- function(overlapDF)
  return(as.list(data.table::transpose(overlapDF[, c(1, 2)])))

overlapSlice <- function(overlapDF, pairs)
  return(overlapDF[which(overlapPairs(overlapDF) %in% pairs),])

