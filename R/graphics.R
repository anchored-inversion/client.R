#' Plot a 1D field curve
#'
#' Plot a 1D field curve given the field value vector and model grid definition.
#'
#' One scenario where this function is useful happens in synthetic studies,
#' where the true field is known (constructed by the user).
#' In this case one would like to compare the model simulated field realizations,
#' requested from the server, against the truth.
#' Note that the simualtions returned by the server are in model scale.
#' The plots of simulations shown on the website \url{www.anchored-inversion.com}
#' are in model scale.
#' To be compared with these plots, the parameter \code{field} to this function
#' should be in model scale.
#'
#' @param field A vector of the field values.
#' @param grid Grid definition. See `grid.R`.
#' @param ... Additional paramters passed to `lattice::xyplot`.
#'
#' @export
plot_field_1d <- function(field, grid, xlab = 'X', ylab = 'Y', type = 'l', ...)
{
    lattice::xyplot(field ~ grid_fullxyz(grid),
           type = type,
           xlab = xlab,
           ylab = ylab,
           ...)

    # To display the returned object, say `z`, do something like
    #    trellis.device(color = FALSE, width = 5, height = 4)
    #    print(z)
    #
    # To overlay some markers on a dashed curve, do something like this:
    #
    #    plot.field.1d(field, grid,
    #                  panel = function(...) {
    #                      panel.xyplot(...)
    #                      panel.points(marker.x, marker.y, pch=20)
    #                      },
    #                  type = 'l',
    #                  lty = 'dashed'
    #                 )
}



plot_field_2d <- function(
    field, grid, xlab = 'X', ylab = 'Y',
    col.cuts = 99,
    col.palette = terrain.colors,
    col.at = NULL,
    ...)
{
    xy <- grid_fullxyz(grid)
    if (is.null(col.at)) {
        dd <- diff(range(field)) * 2
        col.at <- seq(min(field) - dd, max(field) + dd, length = col.cuts)
    }
    lattice::levelplot(c(field) ~ xy[, 1] * xy[, 2],
              col.regions = col.palette(col.cuts + 1),
              cuts = col.cuts,
              at = col.at,
              xlab = xlab,
              ylab = ylab,
              ...
             )

    # To overlay some markers in the field, do something like this:
    #    plot.field.2d(field, grid,
    #                  panel = function(...) {
    #                      panel.levelplot(...)
    #                      panel.points(
    #                          marker.x, marker.y,
    #                          pch = 19, col = 'black')
    #                      }
    #                        # May use more 'panel.points' for other sets of markers.
    #                 )
}

