#' Find the connected components of the graph determined by the overlaps
#'
#' This function finds the connected components of the graph having the filtered
#' overlaps as edges
#'
#' @inheritParams computeCellScores
#' @param raiseWarning If the data frame contains more overlaps than this number,
#' users will be warned that they may have introduced the raw overlap data frame
#' as input
#'
#' @return An overlap data frame with a column indicated the number of the
#' connected component
#'
#' @export
#'
connectedComponents <- function(overlapDF, raiseWarning = 2000){
  if (nrow(overlapDF) > raiseWarning)
    warning(paste0('The number of overlaps in the data frame is very large (', nrow(overlapDF),
    '). Are you sure you filtered the overlap data frame?'))
  if(!nrow(overlapDF))
    stop('Error: The dataframe has no rows.')
  overlapDF$component <- -1
  rownames(overlapDF) <- 1:dim(overlapDF)[1]
  vertices <- union(overlapDF$gene1, overlapDF$gene2)
  seen <- c()
  nextComp = 1
  for (v in vertices){
    if (v %in% seen)
      next
    currVertices <- c(v)
    while (length(currVertices)){
      v <- currVertices[1]
      leftdf <- subset(overlapDF, gene1 == v)
      rightdf <- subset(overlapDF, gene2 == v)
      seen <- c(seen, v)
      newEdges <- as.integer(c(rownames(leftdf), rownames(rightdf)))
      overlapDF$component[newEdges] <- nextComp
      neighbors <- setdiff(c(leftdf$gene2, rightdf$gene1), c(currVertices, seen))
      currVertices <- c(currVertices, neighbors)
      currVertices <- currVertices[-1]
    }
    nextComp <- nextComp + 1
  }
  return(overlapDF)
}

#' Run CSOA separately on the connected components of the overlap graph
#'
#' This function runs CSOA on connected components of the graph having the
#' filtered overlaps as edges
#'
#' @inheritParams runCSOA
#' @param overlapDF An overlap data frame with connected component indices added
#' to the "component" column
#' @inheritParams computeCellScores
#'
#' @return An overlap data frame with a column indicated the number of the
#' connected component
#'
#' @export
#'
scoreModules <- function(scObj, overlapDF, colStr = 'Module'){
  expression <- expMat(scObj)
  scoredComponents <- lapply(seq_len(max(overlapDF$component)), function(i){
    message(str_c('Scoring ', colStr, i, '...'))
    overlapComp <- subset(overlapDF, component == i)
    overlapComp <- scoreOverlaps(overlapDF)
    genes <- unique(union(overlapComp[, 1], overlapComp[, 2]))
    message('Normalizing expression matrix by rows...')
    normExp <- kerntools::minmax(expression[genes, ], rows=T)
    scoreDF <- computeCellScores(overlapDF, normExp, colnames(scObj), colStr)
    return(scoreDF)
  })
  return(joinCellScores(scObj, scoredComponents))
}
