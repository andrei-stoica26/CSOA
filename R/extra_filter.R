#' Compute the Jaccard similarity index of two sets
#'
#' This function computes the Jaccard similarity index of two sets
#'
#' @param a A vector
#' @param b A vector
#' @return A numeric value
#'
#'
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return (intersection/union)
}

#' Compute all the neighbors of all genes in an overlap data frame
#'
#' This function computes all the neighbors of all genes in an overlap data frame
#'
#' @param overlapDF An overlap data frame
#' @return A numeric value
#'
#'
geneNeighbors <- function(overlapDF){
  genes <- overlapGenes(overlapDF)
  return(sapply(genes, function(gene){
    geneLeft <- subset(overlapDF, gene1 == gene)
    geneRight <- subset(overlapDF, gene2 == gene)
    neighbors <- c(geneLeft$gene2, geneRight$gene1)
    return(neighbors)
  }))
}

#' Compute the Jaccard score of the neighbor sets of each gene in an overlap
#'
#' This function computes the Jaccard score of the neighbor sets of each gene
#' in an overlap for all overlap pairs in a data frame
#'
#' @param overlapDF An overlap data frame
#'
#' @return An overlap data frame with an added column of neighbor Jaccard scores
#'
#' @export
#'
#'
neighborJaccard <- function(overlapDF){
  neighbors <- geneNeighbors(overlapDF)
  overlapDF$neighborJac <- sapply(1:nrow(overlapDF), function(i) {
    neighbors1 <- neighbors[[overlapDF[i, ]$gene1]]
    neighbors2 <- neighbors[[overlapDF[i, ]$gene2]]
    return(jaccard(neighbors1, neighbors2))
    }
    , USE.NAMES = F)

  return(overlapDF)
}

#' Remove overlap pairs with low Jaccard scores
#'
#' This function iteratively removes all overlap pairs with Jaccard score below
#' a fixed cutoff until no overlap pairs can be removed
#'
#' @param overlapDF An overlap data frame
#' @param cutoff A cutoff used in the filtering of edges with low Jaccard scores
#' @param doConnComp Whether to calculate the connected components
#'
#' @return An overlap data frame in which edges with low Jaccard scores have
#' been removed
#'
#' @export
#'
#'
breakWeakTies <- function(overlapDF, cutoff = 1/3, doConnComp = F){
  prevNEdges <- -1
  nEdges <- nrow(overlapDF)
  while(prevNEdges != nEdges){
    overlapDF <- neighborJaccard(overlapDF)
    overlapDF <- subset(overlapDF, neighborJac > cutoff)
    prevNEdges <- nEdges
    nEdges <- nrow(overlapDF)
    message(paste0(prevNEdges - nEdges), ' edges with low Jaccard scores have been removed.')
  }
  overlapDF <- rankOverlaps(overlapDF)
  if (doConnComp)f
    overlapDF <- connectedComponents(overlapDF)
  return(overlapDF)
}
