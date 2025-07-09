#' Calculate gene degrees from edges data frame
#'
#' This function calculates gene degrees from a data frame with columns gene1,
#' gene2 and group
#' @param edgesDF A data frame of edges generated with edgeLists
#'
#' @return A gene degrees data frame
#'
#' @noRd
#'
geneDegreesCore <- function(edgesDF){
  genes <- union(edgesDF$gene1, edgesDF$gene2)
  df <- as.data.frame(table(c(edgesDF$gene1, edgesDF$gene2)))
  colnames(df) <- c('gene', 'nEdges')
  df$group <- as.factor(edgesDF$group[1])
  df <- df[order(df$nEdges, decreasing=TRUE), ]
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
#' @noRd
#'
geneDegrees <- function(edgesDFs){
  dfList <- lapply(edgesDFs, geneDegreesCore)
  df <- do.call(rbind, dfList)
  df <- df[order(df$nEdges, decreasing=TRUE), ]
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
#' @noRd
#'
distFreq <- function(degreesDF){
  message('Finding frequencies of gene degrees...')
  center <- degreesDF$nEdges[1] + 1
  if (degreesDF$nEdges[1] != degreesDF$nEdges[2])
    center <- degreesDF$nEdges[1]
  df <- dplyr::count(degreesDF, nEdges)
  df <- df[order(df$nEdges, decreasing=TRUE), ]
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
#' @noRd
#'
geneCoords <- function(overlapObj, groupNames = NULL, cutoff = NULL){
  edgesDFs <- edgeLists(overlapObj, groupNames, cutoff)
  degreesDF <- geneDegrees(edgesDFs)
  distFreqDF <- distFreq(degreesDF)
  message('Finding gene coordinates...')
  circlePoints <- do.call(rbind, lapply(seq_len(nrow(distFreqDF)), function(i)
    pointsOnCircle(distFreqDF$Dist[i], distFreqDF$Freq[i])))
  df <- cbind(degreesDF[, 1, drop=FALSE], circlePoints, degreesDF[, c(2, 3)])
  df[, 5] <- as.factor(df[, 5])
  return(df)
}

#' Store the radii of the circles and the corresponding number of edges
#'
#' This function store the radii of the circles and the corresponding number
#' of edges
#'
#' @param geneCoordsDF Dataframe wih gene coordinates
#' @param extraCircles Number of circles drawn beyond those needed to include
#' the points representing the genes. Default is 0
#'
#' @return A data frame containing the radius and the number of edges for each
#' circle
#'
#' @keywords internal
#'
circleCoords <- function(geneCoordsDF, extraCircles = 0){
  minDegree <- geneCoordsDF$nEdges[nrow(geneCoordsDF)] - extraCircles
  maxDegree <- geneCoordsDF$nEdges[1]
  nCircles <- maxDegree - minDegree + 1
  hasSharedMax <- geneCoordsDF$nEdges[1] == geneCoordsDF$nEdges[2]
  df <- data.frame(
    x = rep(0, nCircles),
    y = rep(0, nCircles),
    r = seq(nCircles + hasSharedMax - 0.5, hasSharedMax + 0.5, -1),
    nEdges = seq(minDegree, maxDegree)
  )
  return(df)
}
