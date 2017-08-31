require(AnchoredInversionClient)
require(lattice)
options(error = utils::recover)

my.seed <- sample(1000, 1)
#my.seed <- 636
    # tough cases: 636, 702, 356
set.seed(my.seed)


#--- synthetic data ---

NX <- 180
NY <- 120
NXY.smoothing <- 2
NXY.upscaling <- 3

data(HalfDome)
myfield <- HalfDome$data
dd <- dim(myfield)
d1 <- sample(dd[1] - NX, 1)
d2 <- sample(dd[2] - NY, 1)
myfield <- myfield[d1 : (d1+NX), d2 : (d2+NY)]

#data(Appalachian)
#myfield <- Appalachian$data

if (NXY.smoothing > 1)
    myfield <- AnchoredInversionClient::average.image(myfield, NXY.smoothing) $ z

mygrid <- list(
    from = c(0, 0),
    by = c(1, 1),
    len = dim(myfield)
    )

myfield <- myfield / mean(myfield) * (1/5000 + 1/8000)/2
    # Use these values as synthetic 'slowness',
    # which is 1/speed.
    # Realistic values of seismic P wave velocity are
    # 5-8 km/s, that is, 5000-8000 m/s,
    # then a realistic value range of 'slowness' is
    #  1/8000 - 1/5000 s/m.
    # Scale the data to make the average within
    # this range.

f.field.transform <- function(x, reverse = FALSE)
    AnchoredInversionClient::log.transform(x, lower = 1e-10, reverse = reverse)

myfield <- f.field.transform(myfield)


#--- tomography config ---

# Time consumption of forward model is proportional to
# the number of sources.
n.sources.left <- 6
n.sources.right <- 6

# Time consumption of forward model is not affected by
# the number of receivers.
# So, can use a large number of receivers.
n.receivers.left <- 10
n.receivers.right <- 10
n.receivers.top <- 15

# Locate sensors.
# Place the sensors somewhat away from the very borders so that
# if the grid is made coarser, they are not too far out of
# the range of grid points.
yspan <- mygrid$by[2] * (mygrid$len[2] - 6)
ystep <- yspan / n.sources.left
sources.left <- cbind(
    mygrid$from[1] + mygrid$by[1]*2,
    seq(from = mygrid$from[2] + mygrid$by[2]*2.5 + ystep/2,
        by = ystep,
        len = n.sources.left)
    )
ystep <- yspan / n.receivers.left
receivers.left <- cbind(
        sources.left[1,1],
        seq(from = mygrid$from[2] + mygrid$by[2]*2.5 + ystep / 2,
            by = ystep,
            len = n.receivers.left)
        )

yspan <- mygrid$by[2] * (mygrid$len[2] - 6)
ystep <- yspan / n.sources.right
sources.right <- cbind(
    mygrid$from[1] + mygrid$by[1] * (mygrid$len[1] - 3),
    seq(from = mygrid$from[2] + mygrid$by[2]*2.5 + ystep/2,
        by = ystep,
        len = n.sources.right)
    )
ystep <- yspan / n.receivers.right
receivers.right <- cbind(
        sources.right[1,1],
        seq(from = mygrid$from[2] + mygrid$by[2]*2.5 + ystep/2,
            by = ystep,
            len = n.receivers.right)
        )

xspan <- mygrid$by[1] * (mygrid$len[1] - 6)
xstep <- xspan / n.receivers.top
receivers.top <- cbind(
        seq(from = mygrid$from[1] + mygrid$by[1]*2.5 + xstep/2,
            by = xstep,
            len = n.receivers.top),
        mygrid$from[2] + mygrid$by[2] * (mygrid$len[2] - 3)
        )


#--- linear data---
n.linear <- 0
linear.data <- NULL


#--- forward function, nonlinear data ---

f.forward <- function(x,
    grid,
    sources.left, sources.right,
    receivers.left, receivers.right, receivers.top,
    f.field.transform, f.forward.transform,
    parallel = AI.getOptions('parallel'))
{
    f.traveltime <- function(x)
    {
        ttime <- c()
        for (isource in seq_len(nrow(sources.left)))
            ttime <- c(
                ttime,
                AnchoredInversionExamples::traveltime.2d(
                    grid = grid, slowness = x,
                    xy.source = sources.left[isource, ],
                    xy.detectors = rbind(receivers.right,
                        receivers.top)) $traveltime
                )
        for (isource in seq_len(nrow(sources.right)))
            ttime <- c(
                ttime,
                AnchoredInversionExamples::traveltime.2d(
                    grid = grid, slowness = x,
                    xy.source = sources.right[isource, ],
                    xy.detectors = rbind(receivers.left,
                        receivers.top)) $traveltime
                )
        ttime
    }


    x.is.list <- is.list(x)
    if (!x.is.list)
        x <- list(x)
    x <- lapply(x, f.field.transform, rev = TRUE)

    y <- lapply(
            x = x,
            FUN = f.traveltime
            )

    y <- lapply(y,
            function(x) {
                if (any(is.na(x) | is.infinite(x)))
                    rep(NA, length(x))
                else
                    unname(f.forward.transform(x))
                }
            )

    if (x.is.list) y else y[[1]]
}

# Transform forward data to real line.
f.forward.transform <- function(x, reverse = FALSE)
    AnchoredInversionClient::log.transform(x, lower = 1e-10, reverse = reverse)

myfield.full <- myfield
mygrid.full <- mygrid

forward.args.full <- list(
    grid = mygrid.full,
    sources.left = sources.left,
    sources.right = sources.right,
    receivers.left = receivers.left,
    receivers.right = receivers.right,
    receivers.top = receivers.top,
    f.field.transform = f.field.transform,
    f.forward.transform = f.forward.transform)

forward.data <- do.call(f.forward,
    c(list(myfield.full), forward.args.full))
    # 'forward.data' is obtained by 'forward.args.full',
    # i.e. in a fine grid, whereas the forward model
    # used for inversion uses 'forward.args' which
    # corresponds to a coarse grid.

forward.data.groups <- NULL

z <- AnchoredInversionClient::average.image(list(
        x = do.call(seq, lapply(mygrid.full, '[[', 1)),
        y = do.call(seq, lapply(mygrid.full, '[[', 2)),
        z = myfield.full),
        NXY.upscaling)
mygrid <- list(
    from = c(z$x[1], z$y[1]),
    by = c(z$x[2] - z$x[1], z$y[2] - z$y[1]),
    len = c(length(z$x), length(z$y))
    )
myfield <- z$z


forward.args <- modifyList(forward.args.full, list(grid = mygrid))

forward.fun <- function(x) { do.call(f.forward, c(x, forward.args)) }
f.forward <- forward.fun


field_value_range = range(myfield) +
    c(-1, 1) * runif(2, 2, 10) * diff(range(myfield))
    #   c(-1, 1) * diff(range(myfield)),
    # A guessed range of the field values.
    # Use a wide range to make the problem more difficult.
# However, the field is defined on the entire real line.

source('anchor.common.R')

