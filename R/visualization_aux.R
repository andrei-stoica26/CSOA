#'@importFrom dplyr count
#'@include utils.R
#'
NULL

#' Convert a matrix to a data frame suitable for basicHeatmap
#'
#' This function converts a matrix to a long data frame
#' suitable for basicHeatmap
#'
#' @param mat A matrix
#' @param colNames A character vector of size 3 representing
#' the column names
#' of the output data frame
#'
#' @return A data frame suitable for basicHeatmap
#'
#' @noRd
#'
heatmapDF <- function(mat, colNames = c('x', 'y', 'Fill')){
  if (!is.matrix(mat))
    stop('mat must be a matrix')
  mat <- cluster_matrix(mat)
  df <- reshape2::melt(mat, varnames=colNames[c(1, 2)],
                       value.name=colNames[3])
  return(df)
}

#' Prepare overlap data frame for network plot
#'
#' This function prepares a ranked and filtered overlap
#' data frame for network plot.
#'
#' @param overlapDF Overlap data frame.
#' @param rankCol Name of the rank column.
#' @param edgeScale Scaling factor used in generating
#' edge weights.
#'
#' @return A data frame ready to serve as input to
#' networkPlot.
#'
#' @keywords internal
#'
networkPlotDF <- function(overlapDF,
                          rankCol = 'rank',
                          edgeScale = 2){
  preWeight <- log(max(overlapDF[[rankCol]]) /
                     overlapDF[[rankCol]] + 0.01)
  overlapDF$weight <- edgeScale * preWeight /
    max(preWeight)
  overlapDF <- overlapDF[,
                         c('gene1', 'gene2', 'weight')]
  return(overlapDF)
}

#' Show the distribution of cell sets among cells
#'
#' This function returns a matrix that shows the presence of cell sets among
#' cells.
#'
#' @param cellSets A list of character vectors.
#' @param allCells Names of all cells in the dataset.
#'
#' @return A logical matrix with genes as rows and cells as columns.
#'
#' @examples
#' cellSets <- list(c('A', 'H', 'J'),
#' c('B', 'D', 'E', 'F', 'J'),
#' c('C', 'I', 'L'))
#' allCells <- LETTERS[seq(15)]
#' cellDistribution(cellSets, allCells)
#'
#' @export
#'
cellDistribution <- function(cellSets, allCells){
  res <- do.call(rbind, lapply(cellSets, function(x) allCells %in% x))
  rownames(res) <- names(cellSets)
  colnames(res) <- allCells
  return(res)
}

