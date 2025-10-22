#' Extracts the data expression matrix from object
#'
#' This function extracts the data expression matrix from object as a
#' non-sparse matrix. Selected genes can be specified as input.
#'
#' @param scObj A Seurat object, SingleCellExperiment object, or
#' expression matrix.
#' @param ... Additional arguments.
#'
#' @return An expression matrix.
#'
#' @examples
#' library(Seurat)
#' mat <- matrix(0, 6, 4)
#' mat[sample(length(mat), 7)] <- sample(3, 7, TRUE)
#' seuratObj <- CreateSeuratObject(counts = mat)
#' seuratObj <- NormalizeData(seuratObj)
#' expMat(seuratObj)
#'
#' @export
#'
expMat <- function(scObj, ...)
    UseMethod(generic='expMat', object=scObj)

#' Extract the edge list from overlap data frame or list of overlap
#' data frames
#'
#' This function creates a list of data frames with three columns:
#' gene1, gene2 and group. If \code{overlapObj} is an overlap data frame,
#' the groups correspond to the connected components. If it is a
#' list of overlap data frames, the groups must be specified as
#' groupNames.
#'
#' @param overlapObj An overlap data frame or list of overlap data frames.
#' @param ... Additional arguments.
#'
#' @return A list of data frames.
#'
#' @keywords internal
#'
edgeLists <- function(overlapObj, ...)
    UseMethod(generic='edgeLists', object=overlapObj)

#' Attach CSOA scores to object
#'
#' This function attaches the data frame of CSOA scores to the input object.
#'
#' @details If the input object is of the Seurat or SingleCellExpression class,
#' it will be returned with added CSOA scores. Otherwise, a list containing the
#' expression matrix and the CSOA scores data frame will be returned.
#'
#' @inheritParams expMat
#' @param ... Additional arguments.
#'
#' @examples
#' library(Seurat)
#' mat <- matrix(0, 500, 300)
#' rownames(mat) <- paste0('G', seq(500))
#' colnames(mat) <- paste0('C', seq(300))
#' mat[sample(8000)] <- sample(20, 8000, TRUE)
#' seuratObj <- CreateSeuratObject(mat)
#' seuratObj <- NormalizeData(seuratObj)
#' scores <- data.frame(CSOA = runif(300))
#' seuratObj <- attachCellScores(seuratObj, scores)
#' head(seuratObj$CSOA)
#'
#' @export
#'
attachCellScores <- function(scObj, ...)
    UseMethod(generic='attachCellScores', object=scObj)

