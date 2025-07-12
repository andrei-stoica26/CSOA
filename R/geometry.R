#' Generate the coordinates of points on a circle centered at origin
#'
#' This function generates nPoints on a circle of radius r
#' centered at origin.
#'
#' @param r Radius.
#' @param nPoints Number of points.
#'
#' @return A data frame with the coordinates of the points.
#'
#' @noRd
#'
pointsOnCircle <- function(r, nPoints){
    angleOffset <- runif(n=1, min=0, max=2 * pi)
    theta <- 2 * pi / nPoints
    points <- lapply(seq(nPoints),
                     function(k) c(r * cos(k * theta + angleOffset),
                                   r * sin(k * theta + angleOffset)))
    res <- do.call(rbind, points)
    colnames(res) <- c('x', 'y')
    return(res)
}

#' Construct a data frame of segments from a data frame of points
#'
#' This function constructs a data frame of segments from a data frame of
#' points.
#'
#' @param pointsDF A data frame with the x and y coordinates of the points.
#' @inheritParams upperConvexHull
#'
#' @return A data frame of segments.
#'
#' @noRd
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
#' This function constructs a data frame of segments from a data frame of
#' points.
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
#' @noRd
#'
hullToPolygon <- function(hull, xInt, type = 'in'){
    yMax <- max(hull$y)
    if (type == 'in'){
        df <- subset(hull, x <= xInt)
        if(df$x[nrow(df)] != xInt)
            df <- rbind(df, c(xInt, yMax))
        yMin <- min(df$y)
        if (yMax != yMin)
        df <- rbind(df, c(xInt, yMin))
    } else{
        df <- subset(hull, x >= xInt)
        if(df$x[nrow(df)] != xInt)
            df <- rbind(c(xInt, yMax), df)
        yMin <- min(df$y)
        if (yMax != yMin)
            df <- rbind(c(xInt, yMin), df)
    }
    return(df)
}
