#'@importFrom dplyr count
#'@include utils.R
#'
NULL

#' Convert a matrix to a data frame suitable for basicHeatmap
#'
#' This function converts a matrix to a long data frame suitable for basicHeatmap
#'
#' @param mat A matrix
#' @param colNames A character vector of size 3 representing the column names
#' of the output data frame
#'
#' @return A data frame suitable for basicHeatmap
#'
#' @export
#'
heatmapDF <- function(mat, colNames = c('x', 'y', 'Fill')){
  if (!is.matrix(mat))
    stop('mat must be a matrix')
  mat <- cluster_matrix(mat)
  df <- reshape2::melt(mat, varnames = colNames[1:2], value.name = colNames[3])
  return(df)
}

#' Prepare overlap data frame for network plot
#'
#' This function prepares a ranked and filtered overlap data frame for network
#' plot.
#'
#' @inheritParams warnUnfiltered
#' @param rankCol Name of the rank column
#' @param edgeScale Scaling factor used in generating edge weights
#'
#' @return A data frame ready to serve as input to networkPlot
#'
#' @export
#'
networkPlotDF <- function(overlapDF, rankCol = 'rank', edgeScale = 2){
  preWeight <- log(max(overlapDF[[rankCol]]) / overlapDF[[rankCol]] + 0.01)
  overlapDF$weight <- edgeScale * preWeight / max(preWeight)
  overlapDF <- overlapDF[, c('gene1', 'gene2', 'weight')]
  return(overlapDF)
}


#' Build a data frame showing the distribution of cells in sets among all cells
#'
#' This function builds a binary data frame showing the distribution of cells in
#' input sets among all cells in the dataset
#'
#' @param cellSets Cell sets
#' @param allCells Names of all cells in the dataset. If not specified, assumed
#' to be the union of the cell sets
#'
#' @return A data frame representing the distribution of cells from input sets
#' among all the sets in the dataset
#'
#' @export
#'
cellDistribution <- function(cellSets, allCells = Reduce(union, cellSets)){
  df <- data.table::transpose(data.frame(lapply(cellSets, function(x) allCells %in% x)))
  rownames(df) <- names(cellSets)
  colnames(df) <- allCells
  return(as.matrix(df))
}

