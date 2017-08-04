#' Bilinear interpolation.
#'
#' @export
bilinear <- function(
    grid,   # 2D grid definition.
    image,  # Data matrix (or image) defined on 'grid'.
    xyout   # Coord matrix (with 2 col) of locations to interpolate.
    )
{
    ij <- grid.xyz2ijk(grid, xyout, rounded = FALSE)

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
