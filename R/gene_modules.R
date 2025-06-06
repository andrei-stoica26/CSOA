#' Find the connected components of the graph determined by the overlaps
#'
#' This function finds the connected components of the graph having the filtered
#' overlaps as edges
#'
#' @inheritParams warnUnfiltered
#' @return An overlap data frame with a column indicated the number of the
#' connected component
#'
#' @export
#'
connectedComponents <- function(overlapDF){
  warnUnfiltered(overlapDF)
  if(!nrow(overlapDF))
    stop('Error: The dataframe has no rows.')
  overlapDF$component <- -1
  rownames(overlapDF) <- 1:dim(overlapDF)[1]
  vertices <- overlapGenes(overlapDF)
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

#' Join scores for multiple gene sets
#'
#' This function joins multiple data frames with CSOA scores and returns the
#' scored object
#'
#' @inheritParams runCSOA
#' @param scoreDFList A list of data frames with CSOA scores
#'
#' @return An object (Seurat, SingleCellExpression or matrix, depending on the
#' input) containing all the scores
#'
joinCellScores <- function(scObj, scoreDFList){
  allScoresDF <- Reduce(cbind, scoreDFList)
  return(storeCellScores(scObj, allScoresDF))
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
    genes <- overlapGenes(overlapComp)
    message('Normalizing expression matrix by rows...')
    normExp <- kerntools::minmax(expression[genes, ], rows=T)
    scoreDF <- computeCellScores(overlapDF, normExp, colnames(scObj), colStr)
    return(scoreDF)
  })
  return(joinCellScores(scObj, scoredComponents))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for CSOA-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname edgeLists
#' @export
#'
edgeLists.default <- function(overlapObj, ...)
  stop('Unrecognized input type: overlapObj must be a data frame or a list of data frames.')

#' @rdname edgeLists
#' @export
#'
edgeLists.data.frame <- function(overlapObj, ...){
  if (!'component' %in% colnames(overlapObj))
    overlapObj <- connectedComponents(overlapObj)
  components <- unique(overlapObj$component)
  overlapObj <- lapply(components, function(i) {
    df <- subset(overlapObj, component == i)
    df <- df[, c('gene1', 'gene2', 'component')]
    colnames(df)[3] <- 'group'
    return(df)
    })
  names(overlapObj) <- components
  return(overlapObj)
}

#' @param groupNames Names of groups. If provided, must be a vector of the same
#' length as the list of overlap data frames
#' @param cutoff Number of retained edges from each overlap data frame after
#' refiltering. If NULL (as default), no refiltering will be performed
#'
#' @rdname edgeLists
#' @export
#'
edgeLists.list <- function(overlapObj, groupNames, cutoff = NULL, ...){
  overlapObj <- lapply(seq_along(groupNames), function(i) {
    df <- overlapObj[[i]]
    if (!is.null(cutoff))
      df <- df[seq_len(cutoff), ]
    df$group <- groupNames[[i]]
    df <- df[, c('gene1', 'gene2', 'group')]
    return(df)
  })
  names(overlapObj) <- groupNames
  return(overlapObj)
}
