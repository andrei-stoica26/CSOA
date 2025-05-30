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

