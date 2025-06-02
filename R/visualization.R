#'@importFrom Seurat FeaturePlot
#'@importFrom ggeasy easy_remove_axes
#'@importFrom ggforce geom_circle
#'@importFrom ggnewscale new_scale_color new_scale_fill
#'@importFrom ggplot2 aes coord_fixed element_text geom_point ggplot ggtitle labs margin scale_color_discrete scale_color_gradientn theme theme_classic theme_void
#'@importFrom ggraph geom_edge_link geom_node_point geom_node_text ggraph scale_edge_width
#'@importFrom ggrepel geom_text_repel
#'@importFrom graphics par
#'@importFrom tidygraph as_tbl_graph
#'@importFrom viridis scale_color_viridis scale_fill_viridis
#'@importFrom wesanderson wes_palette
#'@include utils.R
#'@include visualization_aux.R
#'
NULL

#' Add an aesthetic title to a plot
#'
#' This function adds an aesthetic title to a ggplot object
#'
#' @param p A ggplot object
#' @param title Plot title
#' @param ... Other arguments passed to element_text
#'
#' @return A ggplot object
#'
#' @export
#'
titlePlot <- function(p, title, ...)
  return(p + ggtitle(title) + theme(plot.title=element_text(hjust=0.5, ...)))

#' Improved feature plot with a highly distinctive color scheme
#'
#' This function customizes the appearance of Seurat's FeaturePlot for improved
#' distinctiveness and aesthetics.
#'
#' @param seuratObj A SeuratObj
#' @param feature Seurat feature
#' @param title Plot title
#' @param wesPal Wes Anderson palette
#' @param wesLow Index of color marking low expression
#' @param wesHigh Index of color marking high expression
#' @param ... Other arguments passed to FeaturePlot
#'
#' @return A ggplot object
#'
#' @export
#'
featureWes <- function(seuratObj, feature, title = feature, wesPal='Royal1', wesLow = 3, wesHigh = 2, ...){
  p <- FeaturePlot(seuratObj, feature, ...)
  p <- titlePlot(p, title, ...)
  p <- p + scale_color_gradientn(colours = wes_palette(wesPal)[c(wesLow, wesHigh)])
  return(p)
}

#' Plot the overlaps as a network
#'
#' This function interprets the overlaps as edges in the graph and plots the
#' graph
#'
#' @inheritParams warnUnfiltered
#' @param title Plot title
#' @param nodePointSize Point size of graph nodes
#' @param nodeTextSize Text size of graph nodes
#' @param ... Additional parameters passed to other functions

#' @return A network plot
#'
#' @export
#'
networkPlot <- function(overlapDF, title = 'Top overlaps network plot', nodePointSize = 10, nodeTextSize = 2.3,
                        ...){
  df <- networkPlotDF(overlapDF)
  tblGraph <- tidygraph::as_tbl_graph(df, directed = FALSE)
  p <- ggraph(tblGraph, layout = "nicely") +
    geom_edge_link(aes(width = weight), color = 'green4') +
    scale_edge_width(range = c(0.1, 1)) +
    geom_node_point(size = nodePointSize, color = 'orange') +
    geom_node_text(aes(label = name), color = 'black', size = nodeTextSize) +
    theme_void() +
    theme(legend.position = 'none')
  p <- titlePlot(p, title, ...)
  return(p)
}

#' Gene circle plot for an overlap data frames
#'
#' This functions draws a gene-circle plot for an overlap data frame.
#'
#' @inheritParams edgeLists.list
#' @param title Plot title
#' @param groupStr Column used for grouping
#'
#' @return A ggplot object
#' @export
#'
geneCirclePlot <- function(overlapObj, groupStr = NULL, groupNames = NULL, cutoff = NULL, title = 'Top overlap genes plot'){
  geneCoordsDF <- geneCoords(overlapObj, groupNames, cutoff)
  circleCoordsDF <- circleCoords(geneCoordsDF)
  message('Plotting genes...')
  legendStep <- as.integer(geneCoordsDF$nEdges[1] / 6) + 1
  p <- ggplot() +
    geom_circle(aes(x0 = x, y0 = y, r = r, fill = nEdges, color = nEdges), data = circleCoordsDF) +
    scale_fill_viridis(option = 'viridis', begin = 0.4, breaks = seq(geneCoordsDF$nEdges[1], 1, -legendStep)) +
    scale_color_viridis(option = 'viridis', begin = 0.4, breaks = seq(geneCoordsDF$nEdges[1], 1, -legendStep), guide = 'none') +
    labs(fill = 'Degree') +
    theme_classic() + easy_remove_axes() + coord_fixed() +
    theme(plot.margin = margin(0, 0, 0, 0), legend.title = element_text(size = 10), legend.text = element_text(size = 10)) +
    geom_text_repel(aes(x, y, label = gene), data = geneCoordsDF, size = 4)
  if (!is.null(groupStr))
    p <- p + new_scale_color() +
    new_scale_fill() +
    geom_point(aes(x, y, color = group), data = geneCoordsDF, size = 1) +
    scale_color_discrete(type = c('red', 'purple1', 'olivedrab1','darkorange1', 'snow', 'thistle1', 'green1','violetred4',
                                  'goldenrod1', 'firebrick4')) +
    labs(color = groupStr) else p <- p + geom_point(aes(x, y), data = geneCoordsDF, color = 'red', size = 0.8)
  p <- titlePlot(p, title)
  return(p)
}
