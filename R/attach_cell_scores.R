#' @param scoreDF Dataframe of CSOA scores
#'
#' @rdname attachCellScores
#' @export
#'
attachCellScores.default <- function(scObj, scoreDF, ...)
    stop('Unrecognized input type: scObj must be a Seurat object with a',
         ' data assay, a SingleCellExperiment with a logcounts assay',
         ' a matrix or a dgCMatrix.')

#' @rdname attachCellScores

#' @return A Seurat object with CSOA scores added to metadata.
#'
#' @export
#'
attachCellScores.Seurat <- function(scObj, scoreDF, ...){
    for (colName in colnames(scoreDF))
        if (colName %in% colnames(scObj[[]]))
            scObj[[]][[colName]] <- c()
    scObj[[]] <- cbind(scObj[[]], scoreDF)
    return(scObj)
}

#' @rdname attachCellScores
#'
#' @return A SingleCellExperiment object with CSOA scores added to
#' \code{colData}.
#'
#' @export
#'
attachCellScores.SingleCellExperiment <- function(scObj, scoreDF, ...){
    for (colName in colnames(scoreDF))
        if (colName %in% colnames(colData(scObj)))
            colData(scObj)[[colName]] <- c()
    colData(scObj) <- cbind(colData(scObj), scoreDF)
    return(scObj)
}

#' @rdname attachCellScores
#'
#' @return A list containing the expression matrix and the CSOA scores data
#' frame.
#'
#' @export
#'
attachCellScores.matrix <- function(scObj, scoreDF, ...)
    return(list(object = scObj, scores = scoreDF))

#' @rdname attachCellScores
#'
#' @return A list containing the expression matrix and the CSOA scores data
#' frame.
#'
#' @export
#'
attachCellScores.dgCMatrix <- function(scObj, scoreDF, ...)
    return(list(object = scObj, scores = scoreDF))
