#' Calculate gene degrees from edges data frame
#'
#' This function calculates gene degrees from a data frame with columns
#' 'gene1', gene2' and 'group'.
#' @param edgesDF A data frame of edges generated with edgeLists.
#'
#' @return A gene degrees data frame.
#'
#' @noRd
#'
geneDegreesCore <- function(edgesDF){
    genes <- union(edgesDF$gene1, edgesDF$gene2)
    df <- as.data.frame(table(c(edgesDF$gene1, edgesDF$gene2)))
    colnames(df) <- c('gene', 'nEdges')
    df$group <- as.factor(edgesDF$group[1])
    df <- df[order(df$nEdges, decreasing=TRUE), ]
    return(df)
}

#' Calculate gene degrees from multiple data frames of edges
#'
#' This function calculates gene degrees from the list of data frames of edges
#' generated with edgeLists.
#' @param edgesDFs A list of data frames of edges generated with edgeLists.
#'
#' @return A gene degrees data frame.
#'
#' @noRd
#'
geneDegrees <- function(edgesDFs){
    dfList <- lapply(edgesDFs, geneDegreesCore)
    df <- do.call(rbind, dfList)
    df <- df[order(df$nEdges, decreasing=TRUE), ]
    return(df)
}
