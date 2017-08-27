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



runoff.1d <- function(
    n.x = 100,
    seed = NULL
  )
{
    if (is.null(seed)) {
        seed <- sample(1000, 1)
    }
    set.seed(seed)


    #----------------------
    #     read in data
    #----------------------

    # data(volcano)
    # myfield <- volcano[, sample(ncol(volcano), 1)]
    #     # Vector of length 87.
    # myfield <- average.line(myfield, 2) $y

    data(Denali)
    i <- sample(length(Denali$data) - n.x, 1)
    myfield <- Denali$data[i : (i+n.x)]
    myfield <- average.line(myfield, 2) $y

    mygrid <- list(
        from = 2,
        by = 4,
        len = length(myfield)
        )

    # For the runoff model to run stably,
    # space step should be large compared to time step.

    myfield <- myfield - (5 * min(myfield) - max(myfield)) / (5 - 1)
        # Scale the data to make the value range span 5 times.
    myfield <- myfield / mean(myfield) * .02
        # Make the average close to .02,
        # so that this data are reasonable as synthetic
        # Manning's roughness coef.


    f.field.transform <- function(x, reverse = FALSE)
        log.transform(x, lower = 1e-10, reverse = reverse)

    myfield <- f.field.transform(myfield)



    #---------------------
    #     linear data
    #---------------------

    n.linear <- 0
    linear.data <- NULL

    #---------------------------------
    #    forward function and data
    #---------------------------------

    rain.config <- list(
            # Can have multiple elements of the same structure.
            # Each element represent one experiment, ie rain event.
        list(
            dt = 1,
                # Time step in seconds.
                # An integer.
                # Make 60 divisible by 'dt'.
            t.total = 180,
                # Total experiment time in minutes.
                # An integer.
            rain.stop = 60,
                # Rain stops at the end of this minute.
                # (Rain starts at the beginning of the whole period.
                # It is useless to consider otherwise.)
                # An integer.
            rain.int = runif(1, 2, 10)
                # Rain intensity in mm/hr, constant during this period.
            ),
        list(
            dt = 1,
            t.total = 240,
            rain.stop = 30,
            rain.int = runif(1, 10, 30)
            )
        )

    runoff.config <- lapply(rain.config,
        function(x) {
            list(
                dt = x$dt,
                nt = ceiling(x$t.total * 60 / x$dt),
                    # Total number of numerical time steps.
                rain = rep(c(x$rain.int, 0),
                        c(x$rain.stop, x$t.total - x$rain.stop) * 60 / x$dt)
                    # Rain intensity vector correponding to each
                    # numerical time step.
                )
        } )

    # What model output to take as forward data?
    forward.config <- list(
            # Have as many elements as 'runoff.config' does.
        list(
            discharge = list(
                space.idx = prod(mygrid$len),
                time.idx = seq(from = 1,
                    by = round(600 / rain.config[[1]]$dt),
                        # One measurement per 10 minutes.
                    to = runoff.config[[1]]$nt)
                    # Time index at which observations are extracted.
                )
            ),
        list(
            flowdepth = list(
                space.idx = prod(mygrid$len),
                time.idx = seq(from = 30,
                    by = round(300 / rain.config[[2]]$dt),
                        # One measurement per 5 minutes.
                    to = runoff.config[[2]]$nt)
                )
            )
        )


    f.runoff <- function(
        x,       # roughness field in its natural unit.
        grid = mygrid,
        runoff.config,
        return.config
        )
    {
        z <- vector('list', length(runoff.config))
        for (i.config in seq_along(z))
        {
            zz <- tryCatch(
                do.call(AnchoredInversionClient::runoff.1d,
                    c(list(
                        dx = grid$by,
                        nx = grid$len,
                        h0 = 0, q0 = 0, qt.ub = 0,
                        rough = x,
                        s0 = .01,
                        return.config = return.config[[i.config]]),
                      runoff.config[[i.config]]
                      )
                    ),
                error = function(...) NULL
                )
            if (is.null(zz))
                z[[i.config]] <- lapply(return.config[[i.config]],
                    function(x) rep(NA, prod(sapply(x, length))))
            else
                z[[i.config]] <- zz
        }

        z
    }


    f.forward.transform <- function(x, reverse = FALSE)
    {
        if (reverse)
            log.transform(x, reverse = TRUE) + 1e-100
        else
            log.transform(x + 1e-100)
                # Guard against '0' values in 'x'.
    }

    forward.fun <- function(
        x, # 'x' is a list of fields in the _transformed_ unit.
        grid,
        forward.config,
        runoff.config
        )
    {
        x.is.list <- is.list(x)
        if (!x.is.list)
            x <- list(x)

        z.runoff <- lapply(
            lapply(x, f.field.transform, rev = TRUE),
            f.runoff,
            grid = grid,
            runoff.config = runoff.config,
            return.config = forward.config
            )

        result <- lapply(z.runoff, function(x)
                    f.forward.transform(unname(unlist(x))))

        if (x.is.list)
            result
        else
            result[[1]]
    }

    forward.args <- list(
        grid = mygrid,
        forward.config = forward.config,
        runoff.config = runoff.config)

    f.forward <- function(x) { do.call(forward.fun, c(list(x), forward.args)) }

    forward.data <- f.forward(myfield)

    # forward.data.groups is to help plotting.
    forward.data.groups <- lapply(
        seq_along(forward.config),
        function(i) {
            paste(
                rep(names(forward.config[[i]]),
                    sapply(forward.config[[i]],
                        function(x) prod(sapply(x, length)))
                    ),
                i)
            }
        )
    forward.data.groups <- as.factor(do.call(c, forward.data.groups))

    list(
        mygrid = mygrid,
        myfield = myfield,
        f.field.transform = f.field.transform,
        f.forward = f.forward,
        f.forward.transform = f.forward.transform,
        forward.data = forward.data,
        forward.data.groups = forward.data.groups,
        linear.data = linear.data
    )
}
