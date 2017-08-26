#' 1-D surface rainfall-runoff finite difference simulator.
#'
#' See Fortran procedure in `../src'.
#'
#' @param dx spatial step size (meters) of the numerical grid.
#' @param nx number of cells.
#' @param dt time step size (seconds).
#' @param nt number of time steps to simulate.
#' @param h0
#'      initial flowdepth.  Scalar, or vector of length \code{nx + 1}.
#' @param q0
#'      initial discharge. Scalar, or vector of length \code{nx + 1}.
#' @param qt.ub
#'      upstream discharge (that is, inflow) at all times.
#'      Scalar, or vector of length \code{nt}.
#' @param rain rainfall input (mm/hr).
#'      Scalar, or \code{nx} by \code{nt} matrix,
#'      one column for each time step.
#' @param rough roughness. Scalar, or vector of length \code{nx}.
#' @param s0 scalar, or vector of length \code{nx}.
#' @param return.config
#'
#' @return
#'      A list with members \code{discharge} and \code{flowdepth},
#'      each being a matrix.
#'      Each row is the requested times at one requested location;
#'      each col is a temporal snapshot (at one time step)
#'      for the requested locations.
#'      Flow velocity can be obtained as \code{discharge / flowdepth}.
#'
#' @export
#' @useDynLib AnchoredInversionClient f_runoff_1d
runoff <- function(
    dx, nx, dt, nt,
    h0 = 0, q0 = 0, qt.ub = 0,
    rain, rough, s0,
    return.config
        # List containing members 'discharge' and 'flowdepth'.
        # If either is missing, don't return that data.
        # Each member contains elements 'space.idx' and 'time.idx'.
    )
{
    if (length(h0) == 1) h0 <- rep(h0, nx + 1)
    if (length(q0) == 1) q0 <- rep(q0, nx + 1)
    if (length(qt.ub) == 1) qt.ub <- rep(qt.ub, nt)
    if (length(rough) == 1) rough <- rep(rough, nx)
    if (length(s0) == 1) s0 <- rep(s0, nx)

    rain <- rain / 1000 / 3600
        # mm/hr --> m/sec
    if (length(rain) == 1)
        rain <- matrix(rain, nx, nt)
    else if (length(rain) == nt)
        rain <- matrix(rain, nrow = nx, ncol = nt, byrow = TRUE)
    else
        stopifnot(is.matrix(rain), dim(rain) == c(nx, nt))

    ret.discharge.space <- return.config$discharge$space.idx
        # May be NULL.
    ret.discharge.time <- return.config$discharge$time.idx
        # May be NULL.

    ret.flowdepth.space <- return.config$flowdepth$space.idx
        # May be NULL.
    ret.flowdepth.time <- return.config$flowdepth$time.idx
        # May be NULL.

    z <- .Fortran(f_runoff_1d,
        dx = as.double(dx),
        nx = as.integer(nx),
        dt = as.double(dt),
        nt = as.integer(nt),
        h0 = as.double(h0),
        q0 = as.double(q0),
        qt_ub = as.double(qt.ub),
        rain = as.double(rain),
        rough = as.double(rough),
        s0 = as.double(s0),
        discharge_space_idx = as.integer(ret.discharge.space),
        discharge_space_n = as.integer(length(ret.discharge.space)),
        discharge_time_idx = as.integer(ret.discharge.time),
        discharge_time_n = as.integer(length(ret.discharge.time)),
        flowdepth_space_idx = as.integer(ret.flowdepth.space),
        flowdepth_space_n = as.integer(length(ret.flowdepth.space)),
        flowdepth_time_idx = as.integer(ret.flowdepth.time),
        flowdepth_time_n = as.integer(length(ret.flowdepth.time)),
        discharge = double(length(ret.discharge.space) *
            length(ret.discharge.time)),
        flowdepth = double(length(ret.flowdepth.space) *
            length(ret.flowdepth.time)),
        flag = integer(1)
        )

    if (z$flag != 0)
        stop('Numerical stability criterion failed')
    else
    {
        z <- z[c('discharge', 'flowdepth')]
        if (length(ret.discharge.space) && length(ret.discharge.time))
        {
            if (length(ret.discharge.space) > 1L &&
                length(ret.discharge.time) > 1L)
                dim(z$discharge) <- c(length(ret.discharge.space),
                    length(ret.discharge.time))
        } else
            z$discharge <- NULL

        if (length(ret.flowdepth.space) && length(ret.flowdepth.time))
        {
            if (length(ret.flowdepth.space) > 1L &&
                length(ret.flowdepth.time) > 1L)
                dim(z$flowdepth) <- c(length(ret.flowdepth.space),
                    length(ret.flowdepth.time))
        } else
            z$flowdepth <- NULL

        z
    }
}
