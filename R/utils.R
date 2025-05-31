#' @importFrom methods is
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
#' @param colStr Name of the column of p-values
#' @param pvalThr p-value threshold
#'
#' @return Dataframe with Benjamini-Yekutieli-corrected p-values
#'
byCorrectDF <- function(df, colStr='pval', pvalThr=0.05){
  df <- df[order(df[[colStr]]), ]
  df$pval_adj <- BY(df[[colStr]], 0.05)$Adjusted.pvalues
  if(!is.null(pvalThr))
    df <- subset(df, pval_adj < pvalThr)
  return(df)
}

#' Get all unorderded pairs of two elements from a vector
#'
#' This function returns all unorderded pairs of two elements from a vector as
#' a list of vectors of length 2
#'
#' @param v A vector
#'
#' @return A list of vectors of length 2
#'
getPairs <- function(v)
  return(utils::combn(v, 2, simplify = F))

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
  angleOffset <- runif(n = 1, min = 0, max = 2 * pi)
  theta <- 2 * pi / nPoints
  points <- lapply(1:nPoints, function(k) c(r * cos(k * theta + angleOffset), r * sin(k * theta + angleOffset)))
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
  layerData <- suppressWarnings(LayerData(seuratObj, layer = layer))
  if (!dim(layerData)[1])
    stop(paste0('The Seurat object has no ', layer, ' layer.'))
  return(layerData)
}

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


#' @rdname expMat
#' @export
#'
expMat.default <- function(scObj)
  stop('Unrecognized input type: scObj must be a Seurat object with a data assay, a SingleCellExperiment with a logcounts assay or an expression matrix.')

#' @rdname expMat
#' @export
#'
expMat.Seurat <- function(scObj)
  return(as.matrix(safeLayerData(scObj, layer='data')))

#' @rdname expMat
#' @export
#'
expMat.SingleCellExperiment <- function(scObj)
  return(assay(scObj, 'logcounts'))

#' @rdname expMat
#' @export
#'
expMat.matrix <- function(scObj)
  return(scObj)

#' @rdname edgeLists
#' @export
#'
