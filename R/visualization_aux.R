#'@importFrom dplyr count
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
#' This function returns a logical matrix that shows the representation of
#' cell sets among all cells.
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

#' Adds a gradient color scale using two wesanderson colors
#'
#' This function a gradient color scale to a ggplot object using a wesanderson
#' palette, an index marking low values, and an index marking high values. The
#' indices are used to select colors from the wesanderson palette of choice.
#'
#' @param p A ggplot object.
#' @param palType Palette type: color or fill, continuous or discrete. Accepted
#' values are 'colorCont', 'fillCont', 'colDis' and 'fillDis'. The function shows
#' a warning and does not change the color scheme if a different value is passed
#' here.
#' @param wesPal A wesanderson palette.
#' @param wesLow Index of color marking low values.
#' @param wesHigh Index of color marking high values.
#' @param ... Arguments passed to other functions.
#'
#' @return A ggplot object with a new color scheme.
#'
#' @keywords internal
#'
wesBinaryGradient <- function(p,
                              palType,
                              wesPal = 'Royal1',
                              wesLow = 3,
                              wesHigh = 2,
                              ...){
  if(!palType %in% c('colorCont', 'fillCont', 'colDis', 'fillDis')){
    warning('Unrecognized palette type (see ?CSOA::wesBinaryGradient for the available palette types). The color scheme will
            not be changed')
    return(p)
  }
  colorPair <- wes_palette(wesPal)[c(wesLow, wesHigh)]
  if(palType == 'colorCont')
    p <- p + scale_color_gradientn(colours = colorPair, ...)
  if(palType == 'fillCont')
    p <- p + scale_fill_gradientn(colours = colorPair, ...)
  if(palType == 'colDis')
    p <- p + scale_color_manual(values = colorPair, ...)
  if(palType == 'fillDis')
    p <- p + scale_fill_manual(values = colorPair, ...)
  return(p)
}
