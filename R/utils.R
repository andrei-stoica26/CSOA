#' @importFrom methods is
#' @importFrom kerntools minmax
#' @importFrom qs qread
#' @importFrom stats runif
#' @importFrom SeuratObject LayerData
#' @importFrom sgof BY
#'
NULL

#' Adjust a dataframe column of p-values with Benjamini-Yekutieli
#'
#' This function performs the Benjamini-Yekutieli correction for multiple
#' testing in a dataframe column of p-values. It also offers an option of
#' filtering the dataframe based on p-values.
#'
#' @param df A dataframe with a column of p-values
#' @param pvalThr p-value threshold
#' @param colStr Name of the column of p-values
#'
#' @return The data frame with Benjamini-Yekutieli-corrected p-values
#'
#' @export
#'
byCorrectDF <- function(df, pvalThr = 0.05, colStr = 'pval'){
  df <- df[order(df[[colStr]]), ]
  df$pvalAdj <- BY(df[[colStr]], pvalThr)$Adjusted.pvalues
  df <- subset(df, pvalAdj < pvalThr)
  return(df)
}

#' Show the distribution of cell sets among cells
#'
#' This function returns a matrix that shows the presence of cell sets among
#' cells.
#'
#' @param cellSets A list of character vectors
#' @param allCells A character vector. If not specified, the union of the cell
#' sets
#'
#' @return A logical matrix with genes as rows and cells as columns
#'
#' @export
#'
cellDistribution <- function(cellSets, allCells = Reduce(union, cellSets)){
  df <- data.table::transpose(data.frame(lapply(cellSets, function(x) allCells %in% x)))
  rownames(df) <- names(cellSets)
  colnames(df) <- allCells
  return(as.matrix(df))
}

#' Get all unordered pairs of two elements from a vector
#'
#' This function returns all unorderded pairs of two elements from a vector as
#' a list of vectors of length 2
#'
#' @param v A vector
#'
#' @return A list of vectors of length 2
#'
#' @export
#'
getPairs <- function(v)
  return(utils::combn(v, 2, simplify=FALSE))

#' Get all genes from an overlap data frame
#'
#' This function gets all genes from an overlap data frame
#'
#' @inheritParams rankOverlaps
#'
#' @return A character vector of genes
#'
#' @export
overlapGenes <- function(overlapDF)
  return(union(overlapDF$gene1, overlapDF$gene2))

#' Generate the coordinates of points on a circle centered at origin
#'
#' This function generates nPoints on a circle of radius r centered at origin
#'
#' @param r Radius
#' @param nPoints Number of points
#'
#' @return A data frame with the coordinates of the points
#'
#' @export
#'
pointsOnCircle <- function(r, nPoints){
  angleOffset <- runif(n=1, min=0, max=2 * pi)
  theta <- 2 * pi / nPoints
  points <- lapply(1:nPoints, function(k) c(r * cos(k * theta + angleOffset),
                                            r * sin(k * theta + angleOffset)))
  df <- data.table::transpose(data.frame(points))
  colnames(df) <- c('x', 'y')
  return(df)
}

#' Run LayerData from Seurat and return an error when the requested layer does
#' not exist
#'
#' This function calls LayerData from Seurat and returns an error when the
#' requested layer does not exist
#'
#' @param seuratObj A Seurat object
#' @param layer Layer
#'
#' @return The output of LayerData if layer exists
#'
safeLayerData <- function(seuratObj, layer){
  layerData <- suppressWarnings(LayerData(seuratObj, layer=layer))
  if (!dim(layerData)[1])
    stop(paste0('The Seurat object has no ', layer, ' layer.'))
  return(layerData)
}

#' Filter matrix using rows and convert the matrix to non-sparse
#'
#' This internal functions selects input rows from a matrix in sorted name order.
#' If the rows parameter is set to NULL.
#'
#' @param matObj Matrix object
#' @param rows Rows
#'
#' @return A non-sparse matrix
#'
matrixRowFilter <- function(matObj, rows = NULL){
  if(!is.null(rows)){
    if(length(setdiff(rows, rownames(matObj))))
      stop('Some input genes do not exist in the expression matrix')
    matObj <- matObj[sort(rows), ]
    return(as.matrix(matObj))
  }
  return(as.matrix(matObj)[sort(rownames(matObj)), ])
}

#' Read and delete a .qs file
#'
#' This functions reads a .qs file, deletes it, and returns it content
#'
#' @param qsFile Name of .qs file including its path
#'
#' @return The content of the .qs file
#'
#' @export
#'
qGrab <- function(qsFile){
  res <- qread(qsFile)
  file.remove(qsFile)
  return(res)
}

#' Helper function to ensure easy testing of different rank methods
#'
#' This function controls the choice of rank functions everywhere in the package.
#'
#' @param v Vector
#'
#' @return Ranked vector
#'
#'
rankFun <- function(v)
  return(rank(v, ties.method='min'))

#' Replace a column by its rank
#'
#' This functions orders a data frame by the values in a column and replaces
#' them by the resulting rank.
#'
#' @param df A data frame
#' @param colName The name of a numeric column
#' @param rankSign 1 to rank the column increasingly, -1 to rank it decreasingly
#'
#' @return The data frame ordered by colName (decreasingly by default), in which the original
#' values of colName have been replaced by ranks
#'
#' @export
#'
rankReplace <- function(df, colName, rankSign = 1){
  df <- df[order(df[, colName], decreasing=rankSign - 1), ]
  df[, colName] <- rankFun(rankSign * df[, colName])
  return(df)
}

#' Applies kerntools minmax-normalization on a vector using
#'
#' This functions applies kerntools::minmax on a vector
#'
#' @param v Numeric vector
#'
#' @return A minmax-normalized vector
#'
#' @export
#'
vMinmax <- function(v)
  return(as.numeric(kerntools::minmax(as.matrix(v))))

#' Raise a warning that the overlap data frame may have been not filtered
#'
#' This function raises a warning that the overlap data frame may have been not
#' filtered based on the number of overlaps
#'
#' @inheritParams computeCellScores
#' @param raiseWarning If the data frame contains more overlaps than this number,
#' users will be warned that they may have introduced the raw overlap data frame
#' as input
#'
warnUnfiltered <- function(overlapDF, raiseWarning = 1000)
  if (nrow(overlapDF) > raiseWarning)
    warning(paste0('The number of overlaps in the data frame is very large (', nrow(overlapDF),
                   '). Are you sure you filtered the overlap data frame?'))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for CSOA-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @param genes Genes retained in the expression matrix. If NULL, all genes will
#' be retained
#'
#' @rdname expMat
#' @export
#'
expMat.default <- function(scObj, genes = NULL, ...)
  stop('Unrecognized input type: scObj must be a Seurat object with a data assay, a SingleCellExperiment with a logcounts assay, a matrix or a dgCMatrix object')

#' @rdname expMat
#' @export
#'
expMat.Seurat <- function(scObj, ...)
  return(matrixRowFilter(safeLayerData(scObj, layer='data'), ...))

#' @rdname expMat
#' @export
#'
expMat.SingleCellExperiment <- function(scObj, ...)
  return(matrixRowFilter(assay(scObj, 'logcounts'), ...))

#' @rdname expMat
#' @export
#'
expMat.dgCMatrix <- function(scObj, ...)
  return(matrixRowFilter(scObj, ...))

#' @rdname expMat
#' @export
#'
expMat.matrix <- function(scObj, ...)
  return(matrixRowFilter(scObj, ...))

