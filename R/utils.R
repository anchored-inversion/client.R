# Bilinear interpolation.
#
bilinear <- function(
    grid,   # 2D grid definition.
    image,  # Data matrix (or image) defined on 'grid'.
    xyout   # Coord matrix (with 2 col) of locations to interpolate.
    )
{
    ij <- grid_xyz2ijk(grid, xyout, rounded = FALSE)

    ix1 <- floor(ij[, 1])
    iy1 <- floor(ij[, 2])
        # Value range is 0, 1,..., nx and 0, 1,..., ny.
        # This indicates the grid point (cell center)
        # _before_ the receiver location.
    rx1 <- ix1 + 1 - ij[, 1]
    ry1 <- iy1 + 1 - ij[, 2]

    # Deal with locations outside of the range of grid points.

    idx <- which(ix1 < 1)
    if (length(idx))
    {
        ix1[idx] <- 1
        a <- 1 - ij[idx, 1]
        rx1[idx] <- (1 + a) / (1 + a + a)
    }
    idx <- which(ix1 >= grid$len[1])
    if (length(idx))
    {
        ix1[idx] <- grid$len[1] - 1
        a <- ij[idx, 1] - grid$len[1]
        rx1[idx] <- a / (1 + a + a)
    }

    idx <- which(iy1 < 1)
    if (length(idx))
    {
        iy1[idx] <- 1
        a <- 1 - ij[idx, 2]
        ry1[idx] <- (1 + a) / (1 + a + a)
    }
    idx <- which(iy1 >= grid$len[2])
    if (length(idx))
    {
        iy1[idx] <- grid$len[2] - 1
        a <- ij[idx, 2] - grid$len[2]
        ry1[idx] <- a / (1 + a + a)
    }

    return(
        image[cbind(ix1, iy1)] * rx1 * ry1 +
        image[cbind(ix1, iy1+1)] * rx1 * (1 - ry1) +
        image[cbind(ix1+1, iy1)] * (1 - rx1) * ry1 +
        image[cbind(ix1+1, iy1+1)] * (1 - rx1) * (1 - ry1)
        )
}



#' Moving average with non-overlapping windows of size \code{Q} along a vector (line).
#' Resultant length of data vector is roughly \code{1/Q} of the original.
#'
#' @export
average_line <- function(obj, Q = 2)
{
    if (is.matrix(obj)) {
        if (ncol(obj) == 2)
            obj <- list(x = obj[, 1], y = obj[, 2])
        else if (ncol(obj) == 1)
            obj <- list(x = seq_along(obj), y <- c(obj))
        else
            stop('wrong input format')
    } else if (is.atomic(obj))
        obj <- list(x = seq_along(obj), y = obj)
    else
        stopifnot(is.list(obj))
            # With two vector members 'x' and 'y'
            # of the same length. Not checked.

    M <- length(obj$x)
    MI <- trunc(M/Q)
    indQ <- 1 : Q

    x <- rep(NA, MI)
    y <- rep(NA, MI)

    for (j in 1:MI) {
        idx <- indQ + (j-1) * Q
        x[j] <- mean(obj$x[idx])
        y[j] <- mean(obj$y[idx], na.rm = TRUE)
    }

    list(x = x, y = y, Q = Q)
}




# This function is taken from the 'fields' package version 6.6.2.
average_image <- function(obj, Q = 2) {
    # fast method to sum over a QXQ block in image.
    # Q is the number of elements to average over in each dimension
    # e.g.  Q=5 --  blocks of 25 values are averaged to one grid cell.
    if (is.matrix(obj)) {
        obj <- list(x = 1:nrow(obj), y = 1:ncol(obj), z = obj)
    }
    M <- length(obj$x)
    N <- length(obj$y)
    Mi <- trunc(M/Q)
    Ni <- trunc(N/Q)
    # space to hold results
    z <- matrix(NA, nrow = Mi, ncol = N)
    x2 <- rep(NA, Mi)
    y2 <- rep(NA, Ni)
    indQ <- 1:Q
    # sum over block of rows and handle x grid values
    for (j in 1:Mi) {
        x2[j] <- mean(obj$x[indQ + (j - 1) * Q])
        z[j, ] <- colMeans(obj$z[indQ + (j - 1) * Q, ], na.rm = TRUE)
    }
    # sum over blocks of columns  and average y grid values
    for (k in 1:Ni) {
        y2[k] <- mean(obj$y[indQ + (k - 1) * Q])
        z[, k] <- rowMeans(z[, indQ + (k - 1) * Q], na.rm = TRUE)
    }
    return(list(x = x2, y = y2, z = z[1:Mi, 1:Ni], Q = Q))
}




#' Print the summary that is returned by `Session$summarize_project`.
#'
#' @param obj The output of `Session$summarize_project`.
#' @param digits passed on to \code{print}.
#' @param \dots additional arguments to \code{print}.
#'
#' @export
print_summary <- function(obj, digits = 3, ...)
{
    cat('Dimensions of model grid:', paste(obj$grid.dim, sep = ' x '), '\n')

    cat('\n')
    cat('  dimension of parameter space:', obj$x.dim, '\n')
    cat('    of which,\n')
    cat('      geostat parameters: ', obj$geost.dim,
        ' (', obj$geost.labels, ')\n')
    cat('      anchors: ', obj$anchor.dim, '\n')
    cat('  number of mixture components retained:', obj$n.mixtures, '\n')

    if (obj$n.models < 2)
    {
        cat('\n')
        cat('Only initial approx is available.\n')
        return(invisible())
    }

    cat('\n')
    cat('Started at     ', obj$date.creation, '\n')
    cat('Last updated at', obj$date.last.update, '\n')
    cat('Number of updates:', obj$n.models - 1, '\n')
    cat('In each update:\n')
    cat('  dimension of conditioning data:', obj$n.forward, '\n')
    if (any(obj$n.forward.reduced != obj$n.forward)) {
        cat('    reduced to:', obj$n.forward.reduced, '\n')
    }

    cat('  reproduction of conditioning data:\n')
    z <- cbind(n = obj$diag$n.samples,
            n.anchors = obj$diag$n.anchors,
            intLL = obj$diag$intLL,
            mad.50 = obj$diag$mad.50,
            mad.75 = obj$diag$mad.75,
            mad.95 = obj$diag$mad.95,
            mad.max = obj$diag$mad.max)
    rownames(z) <- paste('      iter', 1:nrow(z), '  ')
    print(z, digits = digits, ...)
    invisible()
}

