require(AnchoredInversion)
require(AnchoredInversionUtils)
require(AnchoredInversionExamples)
require(lattice)

options(error = utils::recover)

my.seed <- sample(1000, 1)
#my.seed <- 919
    # Tough cases:
set.seed(my.seed)


#----------------------
#     read in data
#----------------------

# data(volcano)
# myfield <- volcano[, sample(ncol(volcano), 1)]
#     # Vector of length 87.
# myfield <- average.line(myfield, 2) $y

data(Denali)
i <- sample(length(Denali$data) - 300, 1)
myfield <- Denali$data[i : (i+300)]
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
    AnchoredInversionUtils::log.transform(x, lower = 1e-10, reverse = reverse)

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
            do.call(runoff.1d,
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
        AnchoredInversionUtils::log.transform(x, reverse = TRUE) + 1e-100
    else
        AnchoredInversionUtils::log.transform(x + 1e-100)
            # Guard against '0' values in 'x'.
}

f.forward <- function(
    x, # 'x' is a list of fields in the _transformed_ unit.
    grid,
    forward.config,
    runoff.config,
    parallel = AI.getOptions('parallel')
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

forward.fun <- function(x) { do.call(f.forward, c(x, forward.args)) }

forward.data <- forward.fun(list(myfield))
    # Forward data.

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


#--- plot the actual field ---

trellis.device(color = FALSE, width = 5, height = 4)
print(xyplot(
    f.field.transform(myfield, rev = TRUE) ~ AnchoredInversionUtils::grid.fullxyz(mygrid),
    type = 'l',
    xlab = 'X',
    ylab = 'Roughness')
    )

trellis.device(color = FALSE, width = 5, height = 4)
print(xyplot(
    myfield ~ grid.fullxyz(mygrid),
    type = 'l',
    xlab = 'X',
    ylab = 'Log roughness')
    )

#--- plot the forward process and measurements ---

tmp.forward.config <- forward.config
for (i.config in seq_along(tmp.forward.config))
    tmp.forward.config[[i.config]][[1]]$time.idx <-
        1 : runoff.config[[i.config]]$nt

z <- f.runoff(
    f.field.transform(myfield, rev = TRUE),
    grid = mygrid,
    runoff.config = runoff.config,
    return.config = tmp.forward.config)
for (i.runoff in seq_along(z))
{
    for (i.forward in seq_along(tmp.forward.config[[i.runoff]]))
    {
        forward.type <- names(tmp.forward.config[[i.runoff]])[i.forward]

        y <- z[[i.runoff]][[forward.type]]
            # if 'flowdepth, 'm'; # if 'discharge', m^2 / s
        x <- seq_along(y) * runoff.config[[i.runoff]]$dt / 60
            # minutes
        if (forward.type == 'discharge')
        {
            y <- y * 60  # m^2 / min
            ylab <- expression(paste('Discharge (', m^2, '/min)', sep = ''))
        } else
        {
            ylab <- 'Flow depth (m)'
        }

        x11(width = 4, height = 4)
        print(xyplot(y ~ x,
            type = 'l',
            xlab = 'Time (min)',
            ylab = ylab,
            main = forward.type
            ))
    }
}


#------------------------------------
# Model inference, summary, plotting.
#------------------------------------

save.file.prefix <- 'runoff'

source('anchor.common.R')

forward.data.x <- do.call(c,
    lapply(forward.config,
        lapply,
        function(x) x$time.idx / 60)
            # Time in 'minutes'
    )

z <- plot.anchorit.forward.sample(ai.out,
    data.forward.groups = forward.data.groups,
    x.at = forward.data.x,
    xlab = 'Time (minutes)')
for (zz in unpack.lattice.plots(z)) {x11(); print(zz)}

rm(z, zz)
