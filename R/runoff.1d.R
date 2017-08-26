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
