#' @importFrom SeuratObject LayerData
#' @importFrom sgof BY
#'
NULL

#' Extracts the data expression matrix from Seurat object
#'
#' This function extracts the "data" expression matrix from Seurat object and
#' converts it to a non-sparse matrix
#'
#' @param seuratObj A Seurat object
#'
#' @return An expression matrix of matrix class
#'
#' @export
#'
expMat <- function(seuratObj)
  return(as.matrix(LayerData(seuratObj, layer='data')))

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
byCorrectDF <- function(df, colStr='pval', pvalThr=0.05){
  df <- df[order(df[[colStr]]), ]
  df$pval_adj <- BY(df[[colStr]], 0.05)$Adjusted.pvalues
  if(!is.null(pvalThr))
    df <- subset(df, pval_adj < pvalThr)
  return(df)
}
