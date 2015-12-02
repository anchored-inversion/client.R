#' Steady state 1-D groundwater flow equation, using Darcy's law.
#'
#' Solve the 1-D ODE: \code{d(K dh/dx) / dx = s}.
#'
#' @param K hydraulic conductivity (m/s).
#'      The possible value range is huge (say 13 orders of magnitude).
#'      Typical value range for good groundwater conditions:
#'      1e-5---1e-1 m/s.
#' @param dx step size (m) of numerical grid.
#' @param source volumetric source rate per unit volume
#'      (m^3/s / m^3, ie 1/s).
#' @param b0.type type of upstream boundary condition.
#' @param b0 value of upstream boundary condition.
#'      If \code{b0.type} is \code{'dirichlet'},
#'      \code{b0} provides a head value (m).
#'      If \code{b0.type} is \code{'neumann'},
#'      \code{b0} provides water flux at boundary.
#'      Flux is volumetric rate per volume; same unit as
#'      \code{source}.
#' @param b1.type type of downstream boundary condition.
#' @param b1 value of downstream boundary condition.
#' @param return.idx indices of grid points to return
#'
#' The model domain is discretized uniformly.
#' \code{K}, \code{s}, \code{h} are all evaluated at the model grids.
#'
#' @return Statedy-state hydraulic head (m) in each grid cell.
#'      Same length as \code{K}.
#'
#' @references
#'  Schwartz and Zhang, Fundamentals of Ground Water, Wiley, 2003.
#' @export
darcy <- function(
    K,
    dx,
    source = rep(0, length(K)),
    b0.type = c('dirichlet', 'neumann'),
    b0,
    b1.type = c('dirichlet', 'neumann'),
    b1,
    return.idx = NULL
    )
{
    b0.type <- match.arg(b0.type)
    b1.type <- match.arg(b1.type)

    n <- length(K)
    A <- matrix(0, n, n)
    # diags <- list(rep(0, n), rep(0, n), rep(0, n))
    b <- matrix(source * dx * dx, ncol = 1)

    for (i in 2 : (n-1))
    {
        A[i, (i-1) : (i+1)] <- c(K[i-1], -(K[i-1] + K[i]), K[i])
        # diags[[1]][i - 1] <- K[i-1]
        # diags[[2]][i] <- -K[i-1] - K[i]
        # diags[[3]][i] <- K[i]
    }

    if (b0.type == 'dirichlet')
    {
        A[1,1] <- 1
        # diags[[2]][1] <- 1
        b[1] <- b0
    } else
    {
        A[1, 1:2] <- c(1, -1)
        # diags[[2]][1] <- 1
        # diags[[3]][1] <- -1
        b[1] <- - b0 * dx
    }

    if (b1.type == 'dirichlet')
    {
        A[n, n] <- 1
        # diags[[2]][n] <- 1
        b[n] <- b1
    } else
    {
        A[n, (n-1) : n] <- c(-1, 1)
        # diags[[1]][n-1] <- -1
        # diags[[2]][n] <- 1
        b[n] <- b1 * dx
    }

#   A <- bandSparse(n, k = c(-1,0,1), diag = diags, symm = FALSE)

        # It turned out using the sparse matrix package (Matrix) makes
        # things much slower. It probably is due to the overhead of the
        # package and the small size of the test problem.

    h <- solve(A, b)
        # This could fail, if the input conditions are problematic.

    # as.vector(h)
    if (is.null(return.idx))
        c(h)
    else
        c(h)[return.idx]
}


