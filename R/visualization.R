#'@importFrom Seurat FeaturePlot
#'@importFrom ggplot2 ggtitle
#'@importFrom ggplot2 theme
#'@importFrom ggplot2 element_text
#'@importFrom ggplot2 scale_colour_gradientn
#'@importFrom wesanderson wes_palette
#'
NULL

#' Add an aesthetic title to a plot
#'
#' This function adds an aesthetic title to a ggplot object
#'
#' @param p A ggplot object
#' @param title Plot title
#' @param titleSize Title font size
#'
#' @return A ggplot object
#'
#' @export
#'
titlePlot <- function(p, title, titleSize=12)
  return(p + ggtitle(title) + theme(plot.title=element_text(hjust=0.5, size=titleSize, color='slateblue', face='bold')))

#' Improved feature plot with a highly distinctive color scheme
#'
#' This function customizes the appearance of Seurat's FeaturePlot for improved
#' distinctiveness and aesthetics.
#'
#' @param seuratObj A SeuratObj
#' @param feature Seurat feature
#' @param title Plot title
#' @param titleSize Title font size
#' @param wesPal Wes Anderson palette
#' @param wesLow Index of color marking low expression
#' @param wesHigh Index of color marking high expression
#'
#' @return A ggplot object
#'
#' @export
#'
featureWes <- function(seuratObj, feature, title, titleSize=12, wesPal='Royal1', wesLow=3, wesHigh=2){
  p <- FeaturePlot(seuratObj, feature)
  p <- titlePlot(p, title, titleSize)
  p <- p + scale_colour_gradientn(colours = wes_palette(wesPal)[c(wesLow, wesHigh)])
  return(p)
}
