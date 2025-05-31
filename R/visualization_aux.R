#'@importFrom dplyr count
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
#'
#'
networkPlotDF <- function(overlapDF, weightFactor = 2, raiseWarning = 1000){
  warnUnfiltered(overlapDF, raiseWarning)
  preWeight <- log(max(overlapDF$rank) / overlapDF$rank + 0.01)
  df$weight <- weightFactor * preWeight / max(preWeight)
  df <- df[, c('gene1', 'gene2', 'weight')]
  return(df)
}

#' Compute the gene degrees in a filtered overlap data frame
#'
#' This function computes the gene degrees in a filtered overlap data frame.
#'
#' @inheritParams warnUnfiltered
#' @param dataset Used when analyzing multiple datasets
#' @param cutoff Allows further filtering of the overlap data frame
#'
#' @return A data frame showcasing the gene degrees
#'
#' @export
#'
geneDegreesCore <- function(overlapDF, dataset = NULL, cutoff = NULL, raiseWarning = 1000){
  warnUnfiltered(overlapDF, raiseWarning)
  if(!is.null(cutoff))
    overlapDF <- overlapDF[1:cutoff, ]
  genes <- union(overlapDF$gene1, overlapDF$gene2)
  df <- data.table::transpose(data.frame(sapply(genes, function(x) {
    geneGraph <- subset(overlapDF, gene1 == x | gene2 == x)
    return(c(nrow(geneGraph), geneGraph$component[1]))
  })))
  df <- cbind(genes, df)
  colnames(df) <- c('gene', 'nEdges', 'component')
  df <- df[order(df$nEdges, decreasing = TRUE), ]
  if (!is.null(dataset))
    df$dataset <- dataset
  return(df)
}

#' Compute the gene degrees in several filtered overlap data frame
#'
#' This function computes the gene degrees several filtered overlap data frames.
#'
#' @param overlapDFs List of overlap data frames
#' @param datasets Character vector containing the names used for the datasets
#' @param cutoffs Numeric vector allowing further filtering of the data frames
#'
#' @return A data frame showcasing the degrees of genes from multiple overlap
#' data frames
#'
#' @export
#'
geneDegrees <- function(overlapDFs, datasets, cutoffs = NULL){
  message('Finding gene degrees...')
  #Check if input is a single overlap data frame
  if (is(overlapDFs) == 'data.frame')
    return(geneDegreesCore(overlapDFs))
  dfList <- mapply(geneDegreesCore, overlapDFs, datasets, cutoffs, SIMPLIFY = FALSE)
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
#' @inheritParams geneDegrees
#'
#' @return A data frame containing the coordinates of the genes
#'
geneCoords <- function(overlapDFs, datasets, cutoffs = NULL){
  degreesDF <- geneDegrees(overlapDFs, datasets, cutoffs)
  distFreqDF <- distFreq(degreesDF)
  message('Finding gene coordinates...')
  df <- data.frame(matrix(nrow = 0, ncol = 3))
  for (i in 1:nrow(distFreqDF))
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
#' @param coordsDF Dataframe wih gene coordinates
#' @param extraCircles Number of circles drawn beyond those needed to include the
#' points representing the genes. Default is 0
#'
#' @return A data frame containing the coordinates of the genes
#'
circlesInfo <- function(coordsDF, extraCircles = 0){
  minDegree <- coordsDF$nEdges[nrow(coordsDF)] - extraCircles
  maxDegree <- coordsDF$nEdges[1]
  nCircles <- maxDegree - minDegree + 1 + extraCircles
  hasSharedMax <- coordsDF$nEdges[1] == coordsDF$nEdges[2]
  df <- data.frame(
    x = rep(0, nCircles),
    y = rep(0, nCircles),
    r = seq(nCircles + hasSharedMax - 0.5, hasSharedMax + 0.5, -1),
    nEdges = seq(minDegree, maxDegree)
  )
  return(df)
}


