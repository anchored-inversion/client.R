#' Seismic traveltime calculation by finite difference.
#' An R interface for the C implementation.
#'
#' This is an R interface to a C implementation.
#'
#' @param grid
#'      grid definition.
#' @param slowness
#'      matrix of constant slowness in each grid box.
#'      As row number increases, \code{X} coordinate increases;
#'      as column number increases, \code{Y} coordinate increases.
#' @param xy.source
#'      2-col matrix. Location of wave source.
#'      Must be within the domain of the grid.
#' @param xy.detectors
#'      2-col matrix. Locations of detectors.
#'      If missing, travel time to all grid points are returned
#'      as a matrix; if provided, travel times to these detectors
#'      are returned as a vector.
#' @param inner.radius
#'      size of inner radius of search.
#'
#' @return
#'      A list with members \code{grid}, \code{traveltime} and
#'      \code{flag}.
#'
#' @export
#' @useDynLib AnchoredInversionExamples c_travel_time_2d
traveltime.2d <- function(
    grid, slowness, xy.source, xy.detectors = NULL,
    inner.radius = 3
    )
{
    stopifnot(
        inner.radius > 0,
        length(grid$from) == 2,
        isTRUE(all.equal(grid$by[1], grid$by[2]))
            # Algorithm requires square grid cells.
        )

    dx <- grid$by[1]

    xy.source <- xy.source - rep(grid$from, each = length(xy.source)/2)
        # Algorithm assumes lower-left grid point (or cell
        # center) is (0, 0).

    # FIXME
    # Use .Call to avoid copying 'slowness', and pass
    # 'xy.detectors' to the C side in the meantime to reduce
    # size of the returned array.

    z <- .C(c_travel_time_2d,
        as.integer(nrow(slowness)),
        as.integer(ncol(slowness)),
        as.double(dx),
        as.double(xy.source),
        as.double(slowness),
        as.integer(inner.radius),
        tt = double(length(slowness)),
        nsing = integer(1),
        nneg = integer(1)
        )

    if (is.null(xy.detectors))
        tt <- array(z$tt, dim = dim(slowness))
    else
        tt <- bilinear(grid,
                array(z$tt, dim = dim(slowness)),
                xy.detectors)

    list(grid = grid,
         traveltime = tt,
         flag = list(n.singular = z$nsing, n.negative = z$nneg))
}

