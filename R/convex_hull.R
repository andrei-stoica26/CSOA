
#' Find the indices of points on the left sector of the upper convex hull
#'
#' This function finds the indices of points on left sector of the upper convex
#' hull.
#'
#' @param df A data frame.
#' @param valIndex Index of column storing values used for computing the
#' convex hull.
#'
#' @return A numeric vector of indices.
#'
#'
upperConvexSemihullIndices <- function(df, valIndex = 2){
  lineIndices <- which.min(df[, valIndex])
  maxVal <- df[lineIndices, valIndex]
  if (nrow(df) > 1){
    for (i in seq_len(nrow(df)))
      if (df[i, valIndex] > maxVal){
        lineIndices <- c(lineIndices, i)
        maxVal <- df[i, valIndex]
      }
  }
  return(lineIndices)
}

#' Find the upper convex hull of a set of points
#'
#' This function extracts the upper convex hull of a set of points from a data
#' frame.
#'
#' @inheritParams upperConvexSemihullIndices
#' @param xIndex Index of the column representing the x axis.
#'
#' @return A data frame comprising the points on the upper complex hull.
#'
#' @export
#'
upperConvexHull <- function(df, xIndex = 1, valIndex = 2){
  if (length(colnames(df)) < xIndex)
    stop('xIndex too high; df does not have enough columns')
  if (length(colnames(df)) < valIndex)
    stop('valIndex too high; df does not have enough columns')
  df <- df[order(df[, xIndex]), ]
  leftIndices <- upperConvexSemihullIndices(df, valIndex)
  df <- df[order(df[, xIndex], decreasing=TRUE), ]
  rightIndices <- nrow(df) + 1 - rev(upperConvexSemihullIndices(df, valIndex))
  if (rightIndices[1] == leftIndices[length(leftIndices)])
    leftIndices <- leftIndices[seq_len(length(leftIndices) - 1)]
  df <- df[order(df[, xIndex]), ]
  return(df[c(leftIndices, rightIndices), ])
}

#' Construct a data frame of segments from a data frame of points
#'
#' This function constructs a data frame of segments from a data frame of points.
#'
#' @param pointsDF A data frame with the x and y coordinates of the points.
#' @inheritParams upperConvexHull
#' @param yIndex Index of the column representing the y axis.
#'
#' @return A data frame of segments.
#'
#' @export
#'
pointsToSegments <- function(pointsDF, xIndex = 1, yIndex = 2){
  df <- data.frame(x = pointsDF[seq_len(nrow(pointsDF) - 1), xIndex],
                   y = pointsDF[seq_len(nrow(pointsDF) - 1), yIndex],
                   xEnd = pointsDF[seq(2, nrow(pointsDF)), xIndex],
                   yEnd = pointsDF[seq(2, nrow(pointsDF)), yIndex])
  return(df)
}

#' Construct a data frame of segments from a data frame of points
#'
#' This function constructs a data frame of segments from a data frame of points.
#'
#' @param hull A data frame representing an upper convex hull of a set of points
#' with columns rank and freq.
#' @param rankCutoff Rank cutoff
#' @param type Whether to compute the polygon corresponding to the retained
#' overlaps ('in', default value) or the one corresponding to the discarded
#' overlaps (any other value).
#'
#' @return A data frame of polygon vertices.
#'
#' @export
#'
hullToPolygon <- function(hull, rankCutoff, type = 'in'){
  maxFreq <- max(hull$freq)
  if (type == 'in'){
    df <- subset(hull, rank <= rankCutoff)
    if(df$rank[nrow(df)] != rankCutoff)
      df <- rbind(df, c(rankCutoff, maxFreq))
    df <- rbind(df, c(rankCutoff, min(df$freq)))
  } else{
    df <- subset(hull, rank > rankCutoff)
    if(df$rank[nrow(df)] != rankCutoff)
      df <- rbind(c(rankCutoff, maxFreq), df)
    df <- rbind(c(rankCutoff, min(df$freq)), df)
  }
  return(df)
}
