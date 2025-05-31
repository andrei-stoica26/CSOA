#'@importFrom Seurat FeaturePlot
#'@importFrom ggeasy easy_remove_axes
#'@importFrom ggforce geom_circle
#'@importFrom ggnewscale new_scale_color new_scale_fill
#'@importFrom ggplot2 aes coord_fixed element_text geom_point ggplot ggtitle labs margin scale_color_discrete scale_color_gradientn theme theme_classic
#'@importFrom ggrepel geom_text_repel
#'@importFrom graphics par
#'@importFrom viridis scale_color_viridis scale_fill_viridis
#'@importFrom wesanderson wes_palette
#'@importFrom igraph E graph_from_data_frame layout_nicely
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
#' @param ... Other arguments passed to FeaturePlot
#'
#' @return A ggplot object
#'
#' @export
#'
featureWes <- function(seuratObj, feature, title = feature, titleSize = 12, wesPal='Royal1', wesLow = 3, wesHigh = 2, ...){
  p <- FeaturePlot(seuratObj, feature, ...)
  p <- titlePlot(p, title, titleSize)
  p <- p + scale_color_gradientn(colours = wes_palette(wesPal)[c(wesLow, wesHigh)])
  return(p)
}


#' Plot the overlaps as a network
#'
#' This function interprets the overlaps as edges in the graph and plots the
#' graph
#'
#' @inheritParams warnUnfiltered
#' @param ... Additional parameters passed to the networkPlotDF helper

#' @return A network plot
#'
#'
#'
networkPlot <- function(overlapDF, ...){
  df <- networkPlotDF(overlapDF, ...)
  g <- graph_from_data_frame(df, directed = F)
  par(mar=c(0,0,0,0))
  plot(g, edge.width = E(g)$weight, vertex.size = 15,
       vertex.color = 'orange',
       vertex.frame.color = 'orange',
       vertex.shape = 'circle',
       vertex.label.color = 'black',
       vertex.label.cex = 0.65,
       vertex.label.font = 1,
       edge.color = 'green4',
       layout = layout_nicely(g))
}

#' Gene circle plot for an overlap data frames
#'
#' This functions draws a gene-circle plot for an overlap data frame.
#'
#' @inheritParams circlesInfo
#' @param circlesDF A data frame storing the details of the circles generated
#' using circleInfo
#' @param title Plot title
#' @param groupStr groupStr Column used for grouping
#'
#' @return A ggplot object
#' @export
#'
geneCirclePlot <- function(coordsDF, circlesDF, title, groupStr = NULL){
  legendStep <- as.integer(coordsDF$nEdges[1] / 6) + 1
  p <- ggplot() +
    geom_circle(aes(x0 = x, y0 = y, r = r, fill = nEdges, color = nEdges), data = circlesDF) +
    scale_fill_viridis(option = 'viridis', begin = 0.4, breaks = seq(coordsDF$nEdges[1], 1, -legendStep)) +
    scale_color_viridis(option = 'viridis', begin = 0.4, breaks = seq(coordsDF$nEdges[1], 1, -legendStep), guide = 'none') +
    labs(fill = 'Degree') +
    theme_classic() + easy_remove_axes() + coord_fixed() +
    theme(plot.margin = margin(0, 0, 0, 0), legend.title = element_text(size = 10), legend.text = element_text(size = 10)) +
    geom_text_repel(aes(x, y, label = gene), data = coordsDF, size = 4)
  if (!is.null(groupStr))
    p <- p + new_scale_color() +
    new_scale_fill() +
    geom_point(aes(x, y, color = group), data = coordsDF, size = 1) +
    scale_color_discrete(type = c('red', 'purple1', 'olivedrab1','darkorange1', 'snow', 'thistle1', 'green1','violetred4',
                                  'goldenrod1', 'firebrick4')) +
    labs(color = groupStr) else p <- p + geom_point(aes(x, y), data = coordsDF, color = 'red', size = 0.8)
  p <- titlePlot(p, title)
  return(p)
}

#' Gene circle plot for multiple overlap data frames
#'
#' This functions draws a gene-circle plot for multiple overlap data frames.
#'
#' @inheritParams geneDegrees
#' @inheritParams geneCirclePlot
#'
#' @return A ggplot object
#' @export
#'
plotDegrees <- function(overlapDFs, title, groupStr, datasets, cutoffs = NULL){
  coordsDF <- geneCoords(overlapDFs, datasets, cutoffs)
  colnames(coordsDF)[5] <- 'group'
  circlesDF <- circlesInfo(coordsDF)
  return(geneCirclePlot(coordsDF, circlesDF, title, groupStr))
}
