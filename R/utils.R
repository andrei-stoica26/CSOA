#' @importFrom methods is
#' @importFrom kerntools minmax
#' @importFrom qs2 qs_read qs_save
#' @importFrom SeuratObject LayerData
#'
NULL

#' Get all unordered pairs of two elements from a vector
#'
#' This function returns all unorderded pairs of two elements
#' from a vector.
#'
#' @param v A vector.
#'
#' @return A list of vectors of length 2.
#'
#' @examples
#' v <- c('ASD', 'VBN', 'HJKL')
#' getPairs(v)
#'
#' @export
#'
getPairs <- function(v)
    return(utils::combn(v, 2, simplify=FALSE))

#' Run LayerData from Seurat and return an error
#' when the requested layer does not exist
#'
#' This function calls LayerData from Seurat and
#' returns an error when the requested layer
#' does not exist.
#'
#' @param seuratObj A Seurat object.
#' @param layer Layer.
#'
#' @return The output of LayerData if layer exists.
#'
#' @noRd
#'
safeLayerData <- function(seuratObj, layer){
    layerData <- LayerData(seuratObj, layer=layer)
    if (!dim(layerData)[1])
        stop('The Seurat object has no ', layer, ' layer.')
    return(layerData)
}

#' Filter matrix using rows and convert the matrix to non-sparse
#'
#' This internal functions selects input rows from a
#' matrix in sorted name order. If rows is set to NULL,
#' it selects all rows.
#'
#' @param matObj Matrix object.
#' @param rows Rows.
#'
#' @return A non-sparse matrix.
#'
#' @noRd
#'
matrixRowFilter <- function(matObj, rows = NULL){
    if(!is.null(rows)){
        if(length(setdiff(rows, rownames(matObj))))
        stop('Some input genes do not exist in the',
             ' expression matrix.')
        matObj <- matObj[sort(rows), ]
        return(as.matrix(matObj))
    }
    return(as.matrix(matObj)[sort(rownames(matObj)), ])
}

#' Get all genes from an overlap data frame
#'
#' This function gets all genes from an overlap data frame.
#'
#' @inheritParams processOverlaps
#' @param components A numeric vector representing the connected components
#' of the overlap data frame graph.
#'
#' @return A character vector of genes.
#'
#' @examples
#' overlapDF <- data.frame(gene1 = paste0('G', c(1, 2, 3)),
#' gene1 = paste0('G', c(2, 7, 8)))
#' overlapGenes(overlapDF)
#'
#' @export
#'
overlapGenes <- function(overlapDF, components=NULL){
    if(is.null(components))
        return(union(overlapDF$gene1, overlapDF$gene2))
    return(lapply(components,
                  function(i)
                      overlapGenes(overlapDF[overlapDF$component == i, ])))
}

#' Extract gene pairs from overlap data frame
#'
#' This function extracts the gene pairs from an overlap data frame.
#'
#' @inheritParams processOverlaps
#'
#' @return A list of gene pairs.
#'
#' @examples
#' overlapDF <- data.frame(gene1 = paste0('G', c(1, 2, 3)),
#' gene1 = paste0('G', c(2, 7, 8)))
#' overlapPairs(overlapDF)
#'
#' @export
#'
#'
overlapPairs <- function(overlapDF)
    return(apply(overlapDF, 1, function(x)
        as.character(x[c(1, 2)]), simplify=FALSE))

#' Extract subset defined using gene pairs from overlap matrix
#'
#' This function extracts the subset determined by input
#' gene pairs from an overlap matrix.
#'
#' @inheritParams rankOverlaps
#' @param pairs Gene pairs corresponding to the extracted overlaps.
#'
#' @return An overlap data frame corresponding to the
#' selected gene pairs.
#'
#' @noRd
#'
overlapSlice <- function(overlapDF, pairs)
    return(overlapDF[overlapPairs(overlapDF) %in% pairs,])

#' Read and delete a .qs2 file
#'
#' This functions reads a .qs2 file, deletes it, and returns its content.
#'
#' @param qs2File Name of .qs2 file with path.
#'
#' @return The content of the .qs2 file.
#'
#' @examples
#' library(qs2)
#' qs_save(c(1, 2, 3), 'temp.qs2')
#' qGrab('temp.qs2')
#'
#' @export
#'
qGrab <- function(qs2File){
    res <- qs_read(qs2File)
    file.remove(qs2File)
    return(res)
}

#' Helper function to ensure easy testing of different
#' rank methods
#'
#' This function controls the choice of rank functions
#' everywhere in the package.
#'
#' @param v Vector.
#'
#' @return Ranked vector.
#'
#' @noRd
#'
rankFun <- function(v)
    return(rank(v, ties.method='min'))

#' Replace a column by its rank
#'
#' This functions orders a data frame by the values
#' in a column and replaces them by the resulting rank.
#'
#' @param df A data frame.
#' @param colName The name of a numeric column.
#' @param rankSign 1 to rank the column increasingly, -1
#' to rank it decreasingly.
#'
#' @return The data frame ordered by colName (decreasingly by default),
#' in which the original values of colName have been replaced by ranks.
#'
#' @noRd
#'
rankReplace <- function(df, colName, rankSign = 1){
    df <- df[order(df[, colName], decreasing=rankSign - 1), ]
    df[, colName] <- rankFun(rankSign * df[, colName])
    return(df)
}

#' Applies kerntools minmax-normalization on a vector
#'
#' This functions applies kerntools::minmax on a vector.
#'
#' @param v Numeric vector.
#'
#' @return A minmax-normalized vector.
#'
#' @noRd
#'
vMinmax <- function(v)
    return(as.numeric(kerntools::minmax(as.matrix(v))))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for CSOA-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' @param genes Genes retained in the expression matrix. If NULL, all genes
#' will be retained
#'
#' @rdname expMat
#' @export
#'
expMat.default <- function(scObj, genes = NULL, ...)
    stop('Unrecognized input type: scObj must be a',
         'Seurat object with a data assay,',
         'a SingleCellExperiment with a logcounts assay,',
         'a matrix or a dgCMatrix object.')

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

