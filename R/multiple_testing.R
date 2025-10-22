#' @importFrom sgof BH BY
#'
NULL


#' Perform multiple testing correction and filtering with Bonferroni
#'
#' This function performs the Bonferroni correction for multiple
#' testing in a dataframe column of p-values and filters the data-frame
#' based on p-values.
#'
#' @param df A dataframe with a column of p-values.
#' @param nTests Number of tests. Default to the number of rows in the data
#' frame.
#' @param pvalThr p-value threshold.
#' @param colStr Name of the column of p-values.
#' @param newColStr Name of the column of adjusted p-values that will be
#' created.
#'
#' @return The data frame with Benjamini-Yekutieli-corrected p-values.
#'
#' @keywords internal
#'
bfCorrectDF <- function(df,
                        nTests = nrow(df),
                        pvalThr = 0.05,
                        colStr = 'pval',
                        newColStr = 'pvalAdj'){
    df[[newColStr]] <- df[[colStr]] * nTests
    df <- df[df[, newColStr] < pvalThr, ]
    return(df)
}

#' Perform multiple testing correction by controlling the false discovery rate
#'
#' This function perform multiple testing correction by controlling the false
#' discovery rate.
#'
#' @inheritParams bfCorrectDF
#'
#' @return The data frame with p-values corrected using the method of choice
#' (Benjamini-Hochberg or Benjamini-Yekutieli)
#'
#' @keywords internal
#'
fdrCorrectDF <- function(df,
                         fdrControlFun,
                         pvalThr = 0.05,
                         colStr = 'pval',
                         newColStr = 'pvalAdj'){
    df <- df[order(df[[colStr]]), ]
    df[[newColStr]] <- fdrControlFun(df[[colStr]], pvalThr)$Adjusted.pvalues
    df <- df[df[, newColStr] < pvalThr, ]
    return(df)
}

#' Perform multiple testing correction and filtering with Benjamini-Hochberg
#'
#' This function performs the Benjamini-Hochberg correction for multiple
#' testing in a dataframe column of p-values and filters the data-frame
#' based on p-values.
#'
#' @inheritParams bfCorrectDF
#' @param ... Additional arguments passed to \code{fdrCorrectDF}.
#'
#' @return The data frame with Benjamini-Hochberg-corrected p-values.
#'
#' @keywords internal
#'
bhCorrectDF <- function(df, ...)
    return(fdrCorrectDF(df, BH, ...))

#' Perform multiple testing correction and filtering with Benjamini-Yekutieli
#'
#' This function performs the Benjamini-Yekutieli correction for multiple
#' testing in a dataframe column of p-values and filters the data-frame
#' based on p-values.
#'
#' @inheritParams bfCorrectDF
#' @param ... Additional arguments passed to \code{fdrCorrectDF}.
#'
#' @return The data frame with Benjamini-Yekutieli-corrected p-values.
#'
#' @keywords internal
#'
byCorrectDF <- function(df, ...)
    return(fdrCorrectDF(df, BY, ...))

#' Perform multiple testing correction
#'
#' This function performs correction for multiple testing in a dataframe column
#' of p-values and filters the data-frame based on p-values.
#'
#' @details This function calls \code{bfCorrectDF}, \code{bhCorrectDF} or
#' \code{byCorrectDF}, depending on the selected option. See their
#' documentation for additional parameters.
#'
#' @inheritParams bfCorrectDF
#' @param mtMethod Multiple testing correction method. Options are
#' Bonferroni ('bf'), Benjamini-Hochberg('bh'), and Benjamini-Yekutieli ('by').
#' @param ... Additional arguments passed to the multiple testing correction
#' method.
#'
#' @keywords internal
#'
mtCorrectDF <- function(df, mtMethod = c('bf', 'bh', 'by'), ...){
    method <- match.arg(mtMethod, c('bf', 'bh', 'by'))
    fun <- eval(as.name(paste0(mtMethod, 'CorrectDF')))
    return(fun(df, ...))
}

