#' @importFrom bayesbio jaccardSets
#'
NULL

#' Compute all the neighbors of all genes in an overlap data frame
#'
#' This function computes all the neighbors of all genes in an overlap data
#' frame
#'
#' @param overlapDF An overlap data frame
#' @return A numeric value
#'
#' @noRd
#'
geneNeighbors <- function(overlapDF){
  genes <- overlapGenes(overlapDF)
  neighbors <- lapply(genes, function(gene) {
    leftNeighbors <- overlapDF$gene2[overlapDF$gene1 == gene]
    rightNeighbors <- overlapDF$gene1[overlapDF$gene2 == gene]
    return(union(leftNeighbors, rightNeighbors))
  })
  names(neighbors) <- genes
  return(neighbors)
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
#' @noRd
#'
neighborJaccard <- function(overlapDF){
  neighbors <- geneNeighbors(overlapDF)
  overlapDF$neighborJac <- mapply(function(x, y)
    jaccardSets(neighbors[[x]], neighbors[[y]]),
    overlapDF$gene1, overlapDF$gene2)
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
breakWeakTies <- function(overlapDF, cutoff = 1/3, doConnComp = FALSE){
  prevNEdges <- -1
  nEdges <- nrow(overlapDF)
  message(nEdges, ' overlap', rep('s', nEdges != 1),
          ' have been selected for Jaccard-based filtering.')
  while(prevNEdges != nEdges){
    overlapDF <- neighborJaccard(overlapDF)
    overlapDF <- overlapDF[overlapDF$neighborJac > cutoff, ]
    prevNEdges <- nEdges
    nEdges <- nrow(overlapDF)
    message(prevNEdges - nEdges, ' edge', rep('s', nEdges != 1),
            ' with low neighbor Jaccard scores have been removed.')
  }
  overlapDF <- rankOverlaps(overlapDF)
  if (doConnComp)
    overlapDF <- connectedComponents(overlapDF)
  return(overlapDF)
}
