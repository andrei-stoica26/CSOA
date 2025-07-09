#' Find the connected components of the graph determined by the overlaps
#'
#' This function finds the connected components of the graph having the filtered
#' overlaps as edges
#'
#' @param df A data frame with gene1 and gene2 columns.
#' @param colName Name of the connected components column to be added.
#' @return An overlap data frame with a column indicated the number of the
#' connected component.
#'
#' @export
#'
connectedComponents <- function(df, colName = 'component'){
  warnUnfiltered(df)
  if(!nrow(df))
    stop('The dataframe has no rows.')
  df[[colName]] <- -1
  rownames(df) <- 1:dim(df)[1]
  vertices <- overlapGenes(df)
  seen <- c()
  nextComp = 1
  for (v in vertices){
    if (v %in% seen)
      next
    currVertices <- c(v)
    while (length(currVertices)){
      v <- currVertices[1]
      leftdf <- subset(df, gene1 == v)
      rightdf <- subset(df, gene2 == v)
      seen <- c(seen, v)
      newEdges <- as.integer(c(rownames(leftdf), rownames(rightdf)))
      df[newEdges, colName] <- nextComp
      neighbors <- setdiff(c(leftdf$gene2, rightdf$gene1), c(currVertices, seen))
      currVertices <- c(currVertices, neighbors)
      currVertices <- currVertices[-1]
    }
    nextComp <- nextComp + 1
  }
  return(df)
}

#' Run CSOA separately on the connected components of the overlap graph
#'
#' This function runs CSOA on the connected components of the graph having the
#' filtered overlaps as edges
#'
#' @inheritParams runCSOA
#' @param df A data frame with gene1, gene2 and component columns
#' @param components Vector of connected components that will be scored
#' @param colStrTemplate Character used in the naming of the component gene sets
#' @param ... Additional parameters to other functions
#'
#' @return An object of the same class as scObj with CSOA scores corresponding to
#' the genes defining each connected components assinged for each cell
#'
#' @export
#'
scoreModules <- function(scObj, df, components, colStrTemplate = 'CSOA_component', ...){
  geneSets <- lapply(components, function(i) overlapGenes(subset(df, component == i)))
  geneSetNames <- paste0(colStrTemplate, components)
  return(runCSOAMultiple(scObj, geneSets, geneSetNames, ...))
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
    overlapObj <- connectedComponents(overlapObj, 'group')
  overlapObj <- overlapObj[, c('gene1', 'gene2', 'group')]
  components <- split(overlapObj, overlapObj$group)
  names(components) <- unique(overlapObj$group)
  return(components)
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
