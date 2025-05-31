#'@importFrom dplyr count
#'@include utils.R
#'
NULL

#' Prepare overlap data frame for network plot
#'
#' This function prepares a ranked and filtered overlap data frame for network
#' plot.
#'
#' @inheritParams warnUnfiltered
#' @param weightFactor A scaling factor used for generating edge weights
#'
#' @return A data frame ready to serve as input to networkPlot
#'
#' @export
#'
networkPlotDF <- function(overlapDF, weightFactor = 2, raiseWarning = 1000){
  warnUnfiltered(overlapDF, raiseWarning)
  preWeight <- log(max(overlapDF$rank) / overlapDF$rank + 0.01)
  df$weight <- weightFactor * preWeight / max(preWeight)
  df <- df[, c('gene1', 'gene2', 'weight')]
  return(df)
}

#' Calculate gene degrees from edges data frame
#'
#' This function calculates gene degrees from a data frame with columns gene1,
#' gene2 and group
#' @param edgesDF A data frame of edges generated with edgeLists
#'
#' @return A gene degrees data frame
#'
#' @export
#'
geneDegreesCore <- function(edgesDF){
  genes <- union(edgesDF$gene1, edgesDF$gene2)
  df <- data.table::transpose(data.frame(sapply(genes, function(x){
    geneGraph <- subset(edgesDF, gene1 == x | gene2 == x)
    return(c(x, nrow(geneGraph), geneGraph$group[1]))
    })))
  colnames(df) <- c('gene', 'nEdges', 'group')
  df$nEdges <- as.integer(df$nEdges)
  df <- df[order(df$nEdges, decreasing = TRUE), ]
  return(df)
}

#' Calculate gene degrees from multiple data frames of edges
#'
#' This function calculates gene degrees from the list of data frames of edges
#' generated with edgeLists
#' @param edgesDFs A list of data frames of edges generated with edgeLists
#'
#' @return A gene degrees data frame
#'
#' @export
#'
geneDegrees <- function(edgesDFs){
  dfList <- lapply(edgesDFs, geneDegreesCore)
  df <- Reduce(rbind, dfList)
  df <- df[order(df$nEdges, decreasing = TRUE), ]
  return(df)
}

#' Map gene degrees to distances from the center and find the frequency of these
#' distances
#'
#' This function interprets gene degrees as distances from a center (high
#' degrees = low distances) and calculates the frequencies of these distances.
#' Used later to draw concentric circles with the frequencies representing the
#' number of points on a circle of the same radius
#'
#' @param degreesDF Gene degree data frame
#'
#' @return A data frame of gene distances
#'
#' @export
#'
distFreq <- function(degreesDF){
  message('Finding frequencies of gene degrees...')
  center <- degreesDF$nEdges[1] + 1
  if (degreesDF$nEdges[1] != degreesDF$nEdges[2])
    center <- degreesDF$nEdges[1]
  df <- dplyr::count(degreesDF, nEdges)
  df <- df[order(df$nEdges, decreasing = T), ]
  df$nEdges <- center - df$nEdges
  colnames(df) <- c('Dist', 'Freq')
  return(df)
}

#' Compute the coordinates of genes on the figure made from concentric circles
#'
#' This function computes the coordinates of genes on the figure made from
#' concentric circles
#'
#' @inheritParams edgeLists.list
#'
#' @return A data frame containing the coordinates of the genes
#'
#' @export
#'
geneCoords <- function(overlapObj, groupNames = NULL, cutoff = NULL){
  edgesDFs <- edgeLists(overlapObj, groupNames, cutoff)
  degreesDF <- geneDegrees(edgesDFs)
  distFreqDF <- distFreq(degreesDF)
  message('Finding gene coordinates...')
  df <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in seq_len(nrow(distFreqDF)))
    df <- rbind(df, pointsOnCircle(distFreqDF$Dist[i], distFreqDF$Freq[i]))
  df <- Reduce(cbind, list(degreesDF[, 1, drop = F], df, degreesDF[, 2:3]))
  df[, 5] <- as.factor(df[, 5])
  return(df)
}

#' Store the radii of the circles and the corresponding number of edges
#'
#' This function store the radii of the circles and the corresponding number
#' of edges
#'
#' @param geneCoordsDF Dataframe wih gene coordinates
#' @param extraCircles Number of circles drawn beyond those needed to include the
#' points representing the genes. Default is 0
#'
#' @return A data frame containing the radius and the number of edges for each
#' circle
#'
#' @export
#'
circleCoords <- function(geneCoordsDF, extraCircles = 0){
  minDegree <- geneCoordsDF$nEdges[nrow(geneCoordsDF)] - extraCircles
  maxDegree <- geneCoordsDF$nEdges[1]
  nCircles <- maxDegree - minDegree + 1 + extraCircles
  hasSharedMax <- geneCoordsDF$nEdges[1] == geneCoordsDF$nEdges[2]
  df <- data.frame(
    x = rep(0, nCircles),
    y = rep(0, nCircles),
    r = seq(nCircles + hasSharedMax - 0.5, hasSharedMax + 0.5, -1),
    nEdges = seq(minDegree, maxDegree)
  )
  return(df)
}
