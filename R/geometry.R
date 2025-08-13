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
#' Each point must appear only once.
#' @param joinEnds Whether to join the last point with the first one.
#'
#' @return A data frame of segments.
#'
#' @noRd
#'
pointsToSegments <- function(pointsDF,
                             joinEnds = TRUE){
    df <- data.frame(x = pointsDF[seq_len(nrow(pointsDF) - 1), 1],
                     y = pointsDF[seq_len(nrow(pointsDF) - 1), 2],
                     xEnd = pointsDF[seq(2, nrow(pointsDF)), 1],
                     yEnd = pointsDF[seq(2, nrow(pointsDF)), 2])
    if(joinEnds)
        df <- rbind(df, c(df$xEnd[nrow(df)],
                          df$yEnd[nrow(df)],
                          df$x[1],
                          df$y[1]))
    return(df)
}

#' Construct the convex hull of a set of points
#'
#' This function constructs the convex hull of a set points.
#'
#' @details The points must be provided as a data frame with two columns.
#'
#' @param pointsDF A data frame with the x and y coordinates of the points.
#' @param hullIndices Precalculated hull indices. Default is \code{NULL}: hull
#' indices are not provided, but they are calculated by \code{convexHull}.
#'
#' @return The points on the convex hull of the original set of points.
#'
#' @noRd
#'
convexHull <- function(pointsDF, hullIndices=NULL){
    if(!is.null(hullIndices))
        hull <- pointsDF[hullIndices, c(1, 2)] else
            hull <- pointsDF[chull(pointsDF[, 1], pointsDF[, 2]), c(1, 2)]
        colnames(hull) <- c('x', 'y')
        return(hull)
}

#' Find the coordinates where vertical or horizontal line intersects the hull
#'
#' This function finds the coordinates where vertical or horizontal line
#' intersects the hull.
#'
#'
#' @param df A four-column data frame representing segments.
#' @param axis An integer representing the axis intersected by the vertical or
#' horizontal line, x (1) or y (2).
#' @param axisIntersect The coordinate where the vertical or horizontal line
#' intersects the relevant axis.
#'
#' @return A vector of size two representing the coordinates of the two
#' intersection points between the vertical or horizontal line and the convex
#' hull on the axis different from the input axis.
#'
#' @keywords internal
#'
borderCoords <- function(df, axis, axisIntersect){

    if(is.null(axisIntersect))
        return(c('None', 'None'))

    otherAxis <- axis %% 2 + 1
    axisEnd <- axis + 2
    otherAxisEnd <- otherAxis + 2
    axisVals <- df[, axis]
    axisEndVals <- df[, axisEnd]

    coords <- c(df[axisVals == axisIntersect, ][, otherAxis])

    if (length(coords) < 2){
        df <- df[axisVals < axisIntersect & axisEndVals > axisIntersect |
                     axisVals > axisIntersect & axisEndVals < axisIntersect, ]

        axisVals <- df[, axis]
        axisEndVals <- df[, axisEnd]
        otherAxisVals <- df[, otherAxis]
        otherAxisEndVals <- df[, otherAxisEnd]

        df$diffRatio <- (otherAxisEndVals - otherAxisVals) /
            (axisEndVals  - axisVals)
        df$newCoord <- df$diffRatio * (axisIntersect - axisVals) +
            otherAxisVals

        coords <- c(coords, df$newCoord)
    }

    return(sort(coords))
}
