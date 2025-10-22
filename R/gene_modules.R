#' Run CSOA separately on the connected components of the overlap graph
#'
#' This function runs CSOA on the connected components of the graph having the
#' filtered overlaps as edges.
#'
#' @inheritParams runCSOA
#' @param networkDF A data frame with \code{gene1}, \code{gene2} and component
#' columns.
#' @inheritParams overlapGenes
#' @param colStrTemplate Character used in the naming of the component
#' gene sets.
#' @param ... Additional parameters passed to \code{runCSOAMultiple}.
#'
#' @return An object of the same class as scObj with CSOA scores
#' corresponding to the genes defining each connected components
#' assigned for each cell.
#'
#' @examples
#' mat <- matrix(0, 500, 300)
#' rownames(mat) <- paste0('G', seq(500))
#' colnames(mat) <- paste0('C', seq(300))
#' mat[sample(8000)] <- runif(8000, max=13)
#' genes1 <- paste0('G', seq(100))
#' mat[genes1, 20:50] <- matrix(runif(100 * 31, min = 14, max = 15),
#' nrow = 100, ncol = 31)
#' genes2 <- paste0('G', seq(101, 200))
#' mat[genes2, 70:100] <- matrix(runif(100 * 31, min = 14, max = 15),
#' nrow = 100, ncol = 31)
#' genes <- union(genes1, genes2)
#' mat <- mat[genes, ]
#' overlapDF <- generateOverlaps(mat)
#' overlapDF <- processOverlaps(overlapDF)
#' overlapDF <- henna::connectedComponents(overlapDF)
#' df <- scoreModules(mat, overlapDF, unique(overlapDF$component))[[2]]
#' head(df)
#'
#' @export
#'
scoreModules <- function(scObj, networkDF, components,
                         colStrTemplate = 'CSOA_component', ...){
    geneSets <- overlapGenes(networkDF, components)
    names(geneSets) <- paste0(colStrTemplate, components)
    return(runCSOA(scObj, geneSets, ...))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for CSOA-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname edgeLists
#'
#' @keywords internal
#'
edgeLists.default <- function(overlapObj, ...)
    stop('Unrecognized input type: `overlapObj` must be',
         ' a data frame or a list of data frames.')

#' @rdname edgeLists
#'
#' @keywords internal
#'
edgeLists.data.frame <- function(overlapObj, ...){
    if (!'component' %in% colnames(overlapObj))
        overlapObj <- connectedComponents(overlapObj, 'group')
    overlapObj <- overlapObj[, c('gene1', 'gene2', 'group')]
    components <- split(overlapObj, overlapObj$group)
    names(components) <- unique(overlapObj$group)
    return(components)
}

#' @param groupNames Names of groups. If provided, must be a vector
#' of the same length as the list of overlap data frames.
#' @param cutoff Number of retained edges from each overlap data frame after
#' refiltering. If \code{NULL} (as default), no refiltering will be performed.
#'
#' @rdname edgeLists
#'
#' @keywords internal
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
