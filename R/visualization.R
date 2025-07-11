#'@importFrom Seurat FeaturePlot
#'@importFrom ggeasy easy_remove_axes
#'@importFrom ggforce geom_circle
#'@importFrom ggnewscale new_scale_color new_scale_fill
#'@importFrom ggplot2 aes annotate coord_fixed element_blank element_text geom_hline geom_line geom_point geom_polygon geom_segment geom_tile geom_vline ggplot ggtitle labs margin scale_color_discrete scale_color_gradientn scale_color_manual scale_fill_gradientn scale_fill_manual scale_size_continuous scale_x_continuous scale_y_continuous theme theme_classic theme_minimal theme_void xlim ylim
#'@importFrom ggraph geom_edge_link geom_node_point geom_node_text ggraph scale_edge_width
#'@importFrom ggrepel geom_text_repel
#'@importFrom graphics par
#'@importFrom grDevices chull dev.new dev.off
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

#' Add a centered title to a plot
#'
#' This function adds a centered title to a ggplot object
#'
#' @param p A ggplot object.
#' @param title Plot title.
#' @param ... Other arguments passed to element_text.
#'
#' @return A ggplot object.
#'
#' @noRd
#'
titlePlot <- function(p, title, ...)
  return(p + ggtitle(title) + theme(plot.title=element_text(hjust=0.5, ...)))

#' Improved feature plot with a highly distinctive color scheme
#'
#' This function customizes the appearance of Seurat's FeaturePlot for improved
#' distinctiveness and aesthetics.
#'
#' @param seuratObj A Seurat object.
#' @param feature Seurat feature.
#' @param title Plot title.
#' @param idClass Column to be used for labelling. If NULL, no column-based labels
#' will be generated.
#' @param labelSize Size of labels. Ignored if idClass is NULL.
#' @param titleSize Title size.
#' @inheritParams wesBinaryGradient
#' @param ... Additional arguments passed to FeaturePlot.
#'
#' @return A ggplot object.
#'
#' @examples
#' library(Seurat)
#' mat <- matrix(0, 3000, 800)
#' mat[sample(length(mat), 90000)] <- sample(8, 90000, TRUE)
#' seuratObj <- CreateSeuratObject(counts = mat)
#' seuratObj <- FindVariableFeatures(seuratObj, nfeatures=200)
#' seuratObj <- NormalizeData(seuratObj)
#' seuratObj <- ScaleData(seuratObj)
#' seuratObj <- RunPCA(seuratObj, verbose=FALSE)
#' seuratObj <- RunUMAP(seuratObj, dims=1:20, verbose=FALSE)
#' featureWes(seuratObj, 'Feature3')
#'
#' @export
#'
featureWes <- function(seuratObj, feature,
                       title = feature,
                       idClass = NULL,
                       labelSize = 3,
                       titleSize = 12,
                       wesPal = 'Royal1',
                       wesLow = 3,
                       wesHigh = 2,
                       ...){
  if(is.null(idClass))
    p <- FeaturePlot(seuratObj, feature) else{
      Idents(seuratObj) <- idClass
      p <- FeaturePlot(seuratObj, feature, label=TRUE, label.size=labelSize, ...)
    }
  p <- titlePlot(p, title, size = titleSize)
  p <- wesBinaryGradient(p, 'colorCont', wesPal, wesLow, wesHigh)
  return(p)
}

#' Plot the overlaps as a network
#'
#' This function plots the graph of the overlap data frame, with genes as vertices
#' and overlaps as edges.
#'
#' @inheritParams networkPlotDF
#' @param title Plot title.
#' @param nodePointSize Point size of graph nodes.
#' @param nodeTextSize Text size of graph nodes.
#' @param ... Additional parameters passed to titlePlot.
#'
#' @return A network plot.
#'
#' @examples
#' overlapDF <- data.frame(gene1 = paste0('G', c(1, 2, 5, 6, 7, 17)),
#' gene2 = paste0('G', c(2, 5, 8, 11, 11, 11)),
#' rank = c(1, 1, 3, 3, 3, 3))
#' networkPlot(overlapDF)
#'
#' @export
#'
networkPlot <- function(overlapDF, title = 'Top overlaps network plot', rankCol = 'rank',
                        edgeScale = 2, nodePointSize = 10, nodeTextSize = 2.3, ...){
  df <- networkPlotDF(overlapDF, rankCol, edgeScale)
  tblGraph <- tidygraph::as_tbl_graph(df, directed=FALSE)
  p <- ggraph(tblGraph, layout="nicely") +
    geom_edge_link(aes(width=weight), color='green4') +
    scale_edge_width(range=c(0.1, 0.3)) +
    geom_node_point(size=nodePointSize, color='orange') +
    geom_node_text(aes(label=name), color='black', size=nodeTextSize) +
    theme_void() +
    theme(legend.position='none')
  p <- titlePlot(p, title, ...)
  return(p)
}

#' Radial plot for an overlap data frame
#'
#' This functions draws a radial plot for an overlap data frame to illustrate
#' gene participation in top overlaps.
#'
#' @inheritParams edgeLists.list
#' @param title Plot title.
#' @param groupLegendName The title of the group legend; if NULL, no groups will
#' be distinguished.
#' @inheritParams circleCoords
#'
#' @return A ggplot object.
#'
#' @examples
#' edgesDF <- data.frame(gene1 = paste0('G', c(1, 2, 3, 4, 7, 8, 10,
#' 11, 11, 10, 10, 10)),
#' gene2 = paste0('G', c(2, 5, 1, 8, 4, 9, 12,
#' 13, 14, 13, 16, 14)))
#' edgesDF <- connectedComponents(edgesDF, 'group')
#' geneRadialPlot(edgesDF, 'component', extraCircles=1)
#'
#' @export
#'
geneRadialPlot <- function(overlapObj, groupLegendName = NULL, groupNames = NULL, cutoff = NULL,
                           title = 'Top overlap genes plot',
                           extraCircles = 2){

  geneCoordsDF <- geneCoords(overlapObj, groupNames, cutoff)
  colors <- c('red', 'purple1', 'olivedrab1','darkorange1', 'lavender', 'thistle1',
              'green1','violetred4','goldenrod1', 'firebrick4')
  nColors <- length(colors)
  if(!is.null(groupLegendName) & length(unique(geneCoordsDF$group)) > nColors)
    stop('Only up to ', nColors, ' groups are allowed.')

  circleCoordsDF <- circleCoords(geneCoordsDF, extraCircles)

  legendStep <- as.integer(geneCoordsDF$nEdges[1] / 6) + 1
  p <- ggplot() +
    geom_circle(aes(x0=x, y0=y, r=r, fill=nEdges, color=nEdges),
                data=circleCoordsDF) +
    scale_fill_viridis(option='viridis', begin=0.4,
                       breaks=seq(geneCoordsDF$nEdges[1], 1, -legendStep)) +
    scale_color_viridis(option='viridis', begin=0.4,
                        breaks=seq(geneCoordsDF$nEdges[1], 1, -legendStep), guide='none') +
    labs(fill='Number of top overlaps') +
    theme_classic() + easy_remove_axes() + coord_fixed() +
    theme(plot.margin=margin(0, 0, 0, 0),
          legend.title=element_text(size=10),
          legend.text=element_text(size=10)) +
    geom_text_repel(aes(x, y, label=gene), data=geneCoordsDF, size=3)
  if (!is.null(groupLegendName))
    p <- p + new_scale_color() +
    new_scale_fill() +
    geom_point(aes(x, y, color=group), data = geneCoordsDF, size=0.8) +
    scale_color_discrete(type=colors) +
    labs(color=groupLegendName) else p <- p + geom_point(aes(x, y),
                                                  data=geneCoordsDF,
                                                  color='red',
                                                  size=0.8)
  p <- titlePlot(p, title)
  return(p)
}

#' Plot a simple heatmap
#'
#' This function plots a simple heatmap, with clustering but no dendograms.
#'
#' @param mat A matrix.
#' @param aesNames A character vector of size 3 representing the y, x and fill
#' aes elements.
#' @param title Plot title.
#' @param axisTextSize Axis text size.
#' @inheritParams wesBinaryGradient
#'
#' @return A ggplot object
#'
#' @examples
#' mat <- matrix(0, 10, 20)
#' mat[sample(length(mat), 50)] <- runif(50, max = 2.5)
#' basicHeatmap(mat)
#'
#' @export
#'
basicHeatmap <- function(mat,
                         aesNames = c('x', 'y', 'Score'),
                         title = 'Heatmap',
                         axisTextSize = 7,
                         palType = 'fillCont',
                         wesPal = 'Royal1',
                         wesLow = 3,
                         wesHigh = 2){
  df <- heatmapDF(mat, aesNames)
  p <- ggplot(df, aes(x=.data[[aesNames[2]]],
                      y=.data[[aesNames[1]]],
                      fill=.data[[aesNames[3]]])) +
    geom_tile() +
    theme_minimal() +
    theme(axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_text(size = axisTextSize),
          axis.title=element_blank())
  p <- wesBinaryGradient(p, palType, wesPal, wesLow, wesHigh)
  p <- titlePlot(p, title)
  return(p)
}


#' Plot the selection of overlaps
#'
#' This plots shows the process of selecting the overlap rank cutoff, showcasing
#' the convex hull of the points representing the frequencies of each rank.
#'
#' @inheritParams findRankCutoff
#' @param rankCutoff Rank cutoff
#' @param title Plot title
#'
#' @return A ggplot object
#' @examples
#' freqDF <- data.frame(rank = c(1, 2, 4, 7, 10,
#' 12, 13, 15, 16),
#' n = c(1, 1, 2, 3, 3, 2,
#' 1, 2, 1))
#' overlapCutoffPlot(freqDF, 8.5)
#'
#'
#' @export
#'
overlapCutoffPlot <- function(freqDF, rankCutoff, title = 'Overlap cutoff plot'){
  if(length(setdiff(c('rank', 'n'), colnames(freqDF))))
    stop('freqDF must have columns rank and n.')
  xMin <- min(freqDF$rank)
  xMax <- max(freqDF$rank)
  yMin <- min(freqDF$n)
  yMax <- max(freqDF$n)

  if (nrow(freqDF) < 2)
    stop('overlapCutoffPlot requires at least two points.')

  colors <- c('purple', 'gold')

  hull <- freqDF[sort(chull(freqDF$rank, freqDF$n)), c('rank', 'n')]
  colnames(hull) <- c('x', 'y')
  hullSegments <- pointsToSegments(hull)
  plg <- hullToPolygon(hull, rankCutoff)
  plgOut <- hullToPolygon(hull, rankCutoff, 'out')

  p <- ggplot() +  theme_classic() + xlim(xMin, xMax) + ylim(yMin, yMax) +
    labs(x='Overlap rank', y='Frequency') +
    geom_polygon(data = plg, aes(x, y, fill='Accepted overlaps'), alpha=0.2, ) +
    geom_polygon(data = plgOut, aes(x, y, fill='Discarded overlaps'), alpha=0.2) +
    geom_segment(data = hullSegments, aes(x, y, xend=xEnd, yend=yEnd),
                 color='black', linewidth=0.8) +
    geom_point(data = freqDF, aes(rank, n), color = 'black', size=1, shape=24) +
    scale_fill_manual(values=colors,
                      labels=c('Accepted overlaps', 'Discarded overlaps')) +
    geom_vline(xintercept=rankCutoff,
               color='blue',
               linewidth=0.3,
               linetype='dashed') +
    theme(legend.title=element_blank(),
          legend.position='bottom')
  p <- titlePlot(p, title)
  return(p)
}
