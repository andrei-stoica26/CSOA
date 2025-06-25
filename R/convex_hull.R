
#' Find the indices of points on the left sector of the upper convex hull
#'
#' This function finds the indices of points on left sector of the upper convex
#' hull.
#'
#' @param df A data frame.
#' @param yIndex Index of column storing the y coordinates of the points.
#'
#' @return A numeric vector of indices representing the points in the upper
#' convex semihull.
#'
#'
upperConvexSemihullIndices <- function(df, yIndex = 2){
  lineIndices <- which.min(df[, yIndex])
  maxVal <- df[lineIndices, yIndex]
  if (nrow(df) > 1){
    for (i in seq_len(nrow(df)))
      if (df[i, yIndex] > maxVal){
        lineIndices <- c(lineIndices, i)
        maxVal <- df[i, yIndex]
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
#' @param xIndex Index of the column storing the x coordinates of the points.
#'
#' @return A data frame comprising the points on the upper complex hull.
#'
#' @export
#'
upperConvexHull <- function(df, xIndex = 1, yIndex = 2){
  if (length(colnames(df)) < xIndex)
    stop('xIndex too high; df does not have enough columns')
  if (length(colnames(df)) < yIndex)
    stop('yIndex too high; df does not have enough columns')
  df <- df[, c(xIndex, yIndex)]
  colnames(df) <- c('x', 'y')
  df <- df[order(df$x), ]
  leftIndices <- upperConvexSemihullIndices(df, 2)
  df <- df[order(df$x, decreasing=TRUE), ]
  rightIndices <- nrow(df) + 1 - rev(upperConvexSemihullIndices(df, 2))
  if (rightIndices[1] == leftIndices[length(leftIndices)])
    leftIndices <- leftIndices[seq_len(length(leftIndices) - 1)]
  df <- df[order(df$x), ]
  return(df[c(leftIndices, rightIndices), ])
}

#' Construct a data frame of segments from a data frame of points
#'
#' This function constructs a data frame of segments from a data frame of points.
#'
#' @param pointsDF A data frame with the x and y coordinates of the points.
#' @inheritParams upperConvexHull
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
#' @param xInt The x coordinate of the point where a vertical line will be drawn
#' @param type Whether to compute the polygon corresponding to the retained
#' overlaps ('in', default value) or the one corresponding to the discarded
#' overlaps (any other value).
#'
#' @return A data frame of polygon vertices.
#'
#' @export
#'
hullToPolygon <- function(hull, xInt, type = 'in'){
  yMax <- max(hull$y)
  if (type == 'in'){
    df <- subset(hull, x <= xInt)
    if(df$x[nrow(df)] != xInt)
      df <- rbind(df, c(xInt, yMax))
    df <- rbind(df, c(xInt, min(df$y)))
  } else{
    df <- subset(hull, x > xInt)
    if(df$x[nrow(df)] != xInt)
      df <- rbind(c(xInt, yMax), df)
    df <- rbind(c(xInt, min(df$y)), df)
  }
  return(df)
}
