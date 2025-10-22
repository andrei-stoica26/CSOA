#' @importFrom Seurat FeaturePlot Idents<-
#' @importFrom ggplot2 aes element_blank element_text geom_tile ggplot scale_color_gradientn scale_color_manual scale_fill_gradientn scale_fill_manual theme theme_minimal
#' @importFrom henna centerTitle connectedComponents hullPlot networkPlot radialPlot
#' @importFrom reshape2 melt
#' @importFrom rlang .data
#' @importFrom textshape cluster_matrix
#' @importFrom wesanderson wes_palette
#'
NULL

#' Plot a simple heatmap
#'
#' This function plots a simple heatmap, with clustering but no dendograms.
#'
#' @param mat A matrix.
#' @param aesNames A character vector of length 3 representing the y, x and fill
#' aes elements.
#' @param title Plot title.
#' @param axisTextSize Axis text size.
#' @inheritParams wesBinaryGradient
#' @param ... Other arguments passed to \code{henna::centerTitle}.
#'
#' @return A ggplot object.
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
                         wesHigh = 2,
                         ...){
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
    p <- centerTitle(p, title, ...)
    return(p)
}

#' A feature plot with a more distinctive color scheme.
#'
#' This function customizes the appearance of \code{Seurat::FeaturePlot} for
#' improved distinctiveness and aesthetics.
#'
#' @param seuratObj A Seurat object.
#' @param feature Seurat feature.
#' @param title Plot title.
#' @param idClass Column to be used for labelling. If NULL, no column-based
#' labels will be generated.
#' @param labelSize Size of labels. Ignored if idClass is NULL.
#' @param titleSize Title size.
#' @inheritParams wesBinaryGradient
#' @param ... Additional arguments passed to \code{Seurat::FeaturePlot}.
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
                       labelSize = 3.5,
                       titleSize = 12,
                       wesPal = 'Royal1',
                       wesLow = 3,
                       wesHigh = 2,
                       ...){
    if(is.null(idClass))
        p <- FeaturePlot(seuratObj, feature) else{
            if(!idClass %in% colnames(seuratObj[[]]))
                stop('Grouping variable `', idClass,
                     '` not found in the metadata of the Seurat object.')
            Idents(seuratObj) <- idClass
            p <- FeaturePlot(seuratObj, feature, label=TRUE,
                             label.size=labelSize, ...)
        }
    p <- centerTitle(p, title, size=titleSize)
    p <- wesBinaryGradient(p, 'colorCont', wesPal, wesLow, wesHigh)
    return(p)
}

#' Radial plot for an overlap data frame
#'
#' This function draws a radial plot for an overlap data frame to illustrate
#' gene participation in top overlaps.
#'
#' @details The function can separate genes by groups. The groups can be, for
#' instance, different gene sets, or different connected components of the same
#' overlap data frame. A wrapper around \code{henna::radialPlot}
#'
#' @inheritParams edgeLists.list
#' @param title Plot title.
#' @param degreeLegendTitle The title of the degree legend.
#' @param groupLegendTitle The title of the group legend. If \code{NULL},
#' no groups will be distinguished.
#' @param extraCircles Number of extra circles to be displayed on the plot.
#' @inheritParams edgeLists.list
#' @param ... Additional parameters passed to \code{henna::radialPlot}.
#'
#' @return A ggplot object.
#'
#' @examples
#' edgesDF <- data.frame(gene1 = paste0('G', c(1, 2, 3, 4, 7, 8, 10,
#' 11, 11, 10, 10, 10)),
#' gene2 = paste0('G', c(2, 5, 1, 8, 4, 9, 12,
#' 13, 14, 13, 16, 14)))
#' edgesDF <- henna::connectedComponents(edgesDF, 'group')
#' geneRadialPlot(edgesDF, groupLegendTitle='Component', extraCircles=1)
#'
#' @export
#'
geneRadialPlot <- function(overlapObj,
                           title = 'Top overlap genes plot',
                           degreeLegendTitle = 'Number of top overlaps',
                           groupLegendTitle = 'Group',
                           extraCircles = 2,
                           groupNames = NULL,
                           cutoff = NULL,
                           ...){

    edgesDFs <- edgeLists(overlapObj, groupNames, cutoff)
    degreesDF <- geneDegrees(edgesDFs)
    p <- radialPlot(degreesDF,
                    title,
                    degreeLegendTitle,
                    groupLegendTitle,
                    extraCircles,
                    ...)
    return(p)
}

#' Plot the selection of overlaps
#'
#' This function illustrates the process of selecting the overlap rank cutoff
#' by plotting rank frequencies against ranks and showcasing the convex hull of
#' the rank-frequency points.
#'
#' @details A wrapper around \code{henna::hullPlot}.
#'
#' @param overlapDF Processed overlap data frame created
#' with \code{processOverlaps}.
#' @param title Plot title.
#' @param palette Color palette. Must have two colors, the first one
#' representing accepted overlaps and the other representing discarded overlaps.
#' @param hullWidth Width of the convex hull.
#' @param xLab x axis label.
#' @param yLab y axis label.
#' @param legendLabs Legend labels.
#' @param pointShape Point shape.
#' @param ... Additional arguments passed to \code{henna::hullPlot}.
#'
#' @return A ggplot object.
#'
#' @examples
#' overlapDF <- data.frame(gene1=paste0('G', c(1, 3, 7, 6, 8, 2, 4, 3, 4, 5)),
#' gene2=paste0('G', c(2, 7, 2, 5, 4, 5, 1, 2, 2, 8)),
#' rank=c(1, 2, 3, 4, 4, 6, 7, 7, 7, 10))
#' overlapCutoffPlot(overlapDF)
#'
#' @export
#'
overlapCutoffPlot <- function(overlapDF,
                              title = 'Overlap cutoff plot',
                              palette = c('purple', 'yellow'),
                              hullWidth = 0.8,
                              xLab = 'Overlap rank',
                              yLab = 'Frequency',
                              legendLabs=c('Accepted overlaps',
                                           'Discarded overlaps'),
                              pointShape = 24,
                              ...
){
    if(!'rank' %in% colnames(overlapDF))
        stop('rank column not found in `overlapDF`.',
             ' Use `process_overlaps` to generate it.')
    freqDF <- dplyr::count(overlapDF, rank)

    if (nrow(freqDF) < 2)
        stop('`overlapCutoffPlot` requires at least two points.')

    freqs <- freqDF$n
    if (length(unique(freqs)) < 2)
        stop('`overlapCutoffPlot` requires at least two',
             ' distinct rank frequencies.')

    maxApps <- which(freqs %in% max(freqs))
    if (identical(maxApps, 1)  | identical(maxApps, nrow(freqDF)))
        stop('`overlapCutoffPlot` requires that the',
             ' maximum rank frequency is reached at a non-extremal point.')

        rankCutoff <- findRankCutoff(freqDF)

        p <- hullPlot(freqDF,
                      title,
                      rankCutoff,
                      palette=palette,
                      hullWidth=hullWidth,
                      xLab=xLab,
                      yLab=yLab,
                      legendLabs=legendLabs,
                      pointShape=pointShape,
                      ...)
        return(p)
}

#' Plot the overlaps as a network
#'
#' This function plots the graph of the overlap data frame, with genes as
#' vertices and overlaps as edges.
#'
#' @details A thin wrapper around \code{henna::networkPlot}.
#'
#' @param overlapDF Overlap data frame.
#' @param title Plot title.
#' @param nodeColor The color of nodes. If \code{NULL}, the default
#' \code{henna::networkPlot} color scheme will be used, which uses different
#' colors for nodes belonging to different connected components.
#' @param edgeColor The color of edges.
#' @param ... Additional parameters passed to \code{henna::networkPlot}.
#'
#' @return An overlap network plot.
#'
#' @examples
#' overlapDF <- data.frame(gene1 = paste0('G', c(1, 2, 5, 6, 7, 17)),
#' gene2 = paste0('G', c(2, 5, 8, 11, 11, 11)),
#' rank = c(1, 1, 3, 3, 3, 3))
#' overlapNetworkPlot(overlapDF)
#'
#' @export
#'
overlapNetworkPlot <- function(overlapDF,
                               title = 'Top overlaps network plot',
                               nodeColor = 'orange',
                               edgeColor = 'green4',
                               ...)
    return(networkPlot(overlapDF,
                       title,
                       nodeColor=nodeColor,
                       edgeColor=edgeColor,
                       ...))
