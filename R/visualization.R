#'@importFrom Seurat FeaturePlot
#'@importFrom ggeasy easy_remove_axes
#'@importFrom ggforce geom_circle
#'@importFrom ggnewscale new_scale_color new_scale_fill
#'@importFrom ggplot2 aes coord_fixed element_blank element_text geom_line geom_point geom_tile ggplot ggtitle labs margin scale_color_discrete scale_color_gradientn scale_color_manual scale_fill_gradientn scale_fill_manual scale_x_continuous scale_y_continuous theme theme_classic theme_minimal theme_void
#'@importFrom ggraph geom_edge_link geom_node_point geom_node_text ggraph scale_edge_width
#'@importFrom ggrepel geom_text_repel
#'@importFrom graphics par
#'@importFrom grDevices dev.new dev.off
#'@importFrom reshape2 melt
#'@importFrom SeuratObject Idents
#'@importFrom textshape cluster_matrix
#'@importFrom tidygraph as_tbl_graph
#'@importFrom viridis scale_color_viridis scale_fill_viridis
#'@importFrom wesanderson wes_palette
#'@include utils.R
#'@include visualization_aux.R
#'
NULL

#' @rdname devPlot
#' @export
#'
devPlot.default <- function(plotObject, ...)
  stop('Unrecognized input type: plotObject must be a function, a ggplot object or a list of ggplot objects')

#' @rdname devPlot
#' @export
#'
devPlot.function <- function(plotObject, ...){
  dev.new(noRStudioGD = TRUE)
  print(plotObject(...))
  dev.off()
}

#' @rdname devPlot
#' @export
#'
devPlot.ggplot <- function(plotObject, ...)
  devPlot.function(identity, plotObject)

#' @rdname devPlot
#' @export
#'
devPlot.list <- function(plotObject, ...)
  invisible(lapply(plotObject, devPlot.ggplot))

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

#' Adds a gradient color scale using two Wes Anderson colors
#'
#' This function a gradient color scale to a ggplot object using a Wes Anderson
#' palette, an index marking low values, and an index marking high values. The
#' indices are used to select colors from the Wes Anderson palette of choice.
#'
#' @param p A ggplot object
#' @param wesPal Wes Anderson palette
#' @param wesLow Index of color marking low values
#' @param wesHigh Index of color marking high values
#' @param palType Palette type: color or fill, continuous or discrete. Accepted
#' values are 'colorCont', 'fillCont', 'colDis' and 'fillDis'. The function shows
#' a warning and does not change the color scheme if a different value is passed
#' here
#' @param ... Arguments passed to other functions
#'
#' @return A ggplot object with a new color scheme
#'
#' @export
#'
wesBinaryGradient <- function(p, wesPal = 'Royal1', wesLow = 3, wesHigh = 2, palType = 'colorCont', ...){
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

#' Improved feature plot with a highly distinctive color scheme
#'
#' This function customizes the appearance of Seurat's FeaturePlot for improved
#' distinctiveness and aesthetics.
#'
#' @param seuratObj A SeuratObj
#' @param feature Seurat feature
#' @param title Plot title
#' @param idClass Column to be used for labelling
#' @param label Whether to label the identity classes
#' @param labelSize Size of labels
#' @param repel Whether to make names of labels repel
#' @param titleSize Title size
#' @param ... Other arguments passed to wesBinaryGradient
#'
#' @return A ggplot object
#'
#' @export
#'
featureWes <- function(seuratObj, feature, title = feature, idClass = 'orig.ident',
                       label = TRUE, labelSize = 3, repel = TRUE, titleSize = 12, ...){
  Idents(seuratObj) <- idClass
  p <- FeaturePlot(seuratObj, feature, label = label, label.size = labelSize, repel = repel)
  p <- titlePlot(p, title, size = titleSize)
  p <- wesBinaryGradient(p, ...)
  return(p)
}

#' Plot the overlaps as a network
#'
#' This function interprets the overlaps as edges in the graph and plots the
#' graph
#'
#' @inheritParams warnUnfiltered
#' @param title Plot title
#' @inheritParams networkPlotDF
#' @param nodePointSize Point size of graph nodes
#' @param nodeTextSize Text size of graph nodes
#' @param ... Additional parameters passed to other functions
#'
#' @return A network plot
#'
#' @export
#'
networkPlot <- function(overlapDF, title = 'Top overlaps network plot', rankCol = 'rank', edgeScale = 2,
                        nodePointSize = 10, nodeTextSize = 2.3, ...){
  df <- networkPlotDF(overlapDF, rankCol, edgeScale)
  tblGraph <- tidygraph::as_tbl_graph(df, directed = FALSE)
  p <- ggraph(tblGraph, layout = "nicely") +
    geom_edge_link(aes(width = weight), color = 'green4') +
    scale_edge_width(range = c(0.1, 0.3)) +
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
#' @inheritParams circleCoords
#'
#' @return A ggplot object
#'
#' @export
#'
geneCirclePlot <- function(overlapObj, groupStr = NULL, groupNames = NULL, cutoff = NULL, title = 'Top overlap genes plot',
                           extraCircles = 0){
  geneCoordsDF <- geneCoords(overlapObj, groupNames, cutoff)
  circleCoordsDF <- circleCoords(geneCoordsDF, extraCircles)
  message('Plotting genes...')
  legendStep <- as.integer(geneCoordsDF$nEdges[1] / 6) + 1
  p <- ggplot() +
    geom_circle(aes(x0 = x, y0 = y, r = r, fill = nEdges, color = nEdges), data = circleCoordsDF) +
    scale_fill_viridis(option = 'viridis', begin = 0.4, breaks = seq(geneCoordsDF$nEdges[1], 1, -legendStep)) +
    scale_color_viridis(option = 'viridis', begin = 0.4, breaks = seq(geneCoordsDF$nEdges[1], 1, -legendStep), guide = 'none') +
    labs(fill = 'Number of top overlaps') +
    theme_classic() + easy_remove_axes() + coord_fixed() +
    theme(plot.margin = margin(0, 0, 0, 0), legend.title = element_text(size = 10), legend.text = element_text(size = 10)) +
    geom_text_repel(aes(x, y, label = gene), data = geneCoordsDF, size = 3)
  if (!is.null(groupStr))
    p <- p + new_scale_color() +
    new_scale_fill() +
    geom_point(aes(x, y, color = group), data = geneCoordsDF, size = 0.8) +
    scale_color_discrete(type = c('red', 'purple1', 'olivedrab1','darkorange1', 'lavender', 'thistle1', 'green1','violetred4',
                                  'goldenrod1', 'firebrick4')) +
    labs(color = groupStr) else p <- p + geom_point(aes(x, y), data = geneCoordsDF, color = 'red', size = 0.8)
  p <- titlePlot(p, title)
  return(p)
}

#' Plot the gene pair rank versus the overlap rank
#'
#' This functions plots the gene pair rank versus the overlap rank
#'
#' @param pairScoreDF A dataframe with columns gene1, gene2 (character),
#' overlapRank and pairRank (numeric)
#' @param title Plot title
#' @param pointColor Point color
#' @param pointSize Point size
#' @param labelSize Label size
#'
#' @return A ggplot object
#'
#' @export
#'
birankPlot <- function(pairScoreDF, title = 'Overlap and gene pair ranks', pointColor = 'deeppink3',
                       pointSize = 1.5, labelSize = 3.5){
  pairScoreDF$overlaps <- paste0(pairScoreDF$gene1, '_', pairScoreDF$gene2)
  p <- ggplot(pairScoreDF, aes(overlapRank, pairRank)) + theme_classic() +
    geom_point(color = pointColor, size =  pointSize) +
    geom_text_repel(aes(label = overlaps), size = labelSize) +
    scale_x_continuous(trans = 'reverse') +
    scale_y_continuous(trans = 'reverse') +
    labs(x='Overlap rank', y='Gene pair rank')
  p <- titlePlot(p, title)
  return(p)
}

#' Plot a simple heatmap
#'
#' This function plots a simple heatmap with clustering but no dendograms
#'
#' @param mat A matrix
#' @param aesNames A character vector of size 3 representing the y, x and fill
#' aes elements
#' @param title Plot title
#' @param axisTextSize Axis text size
#' @param palType Palette type
#' @param ... Arguments passed to other functions
#'
#' @return A ggplot object
#'
#' @export
#'
basicHeatmap <- function(mat, aesNames, title = 'Gene expression heatmap', axisTextSize = 7, palType = 'fillCont', ...){
  df <- heatmapDF(mat, aesNames)
  p <- ggplot(df, aes(x = .data[[aesNames[2]]], y = .data[[aesNames[1]]], fill = .data[[aesNames[3]]])) +
    geom_tile() +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = axisTextSize),
          axis.title.y = element_blank())
  p <- wesBinaryGradient(p, palType = palType, ...)
  p <- titlePlot(p, title)
  return(p)
}


#' Plot a simple heatmap
#'
#' This function plots a simple heatmap with clustering but no dendograms
#'
#' @param df A data frame with rank and score columns
#' @param title Plot title
#' @param lineColor Line color
#' @param pointColor Point color
#' @param pointSize Point size
#' @param axisTitleSize Axis title size
#'
#' @return A ggplot object
#'
#' @export
#'
rankScorePlot <- function(df, title, lineColor = 'mediumpurple4', pointColor = 'red', pointSize = 1, axisTitleSize = 10){
  p <- ggplot(df, aes(x = rank, y = score)) + geom_line(color = lineColor) + geom_point(color = pointColor, size = pointSize) +
    labs(x = 'Overlap rank', y = 'Overlap score')
  p <-  p <- p + theme_classic() + theme(axis.text.x = element_text(size = axisTitleSize - 1),
                                         axis.text.y = element_text(size = axisTitleSize - 1),
                                         axis.title = element_text(size = axisTitleSize))
  p <- titlePlot(p, title)
  return(p)
}

