require(AnchoredInversion)
require(AnchoredInversionExamples)
require(lattice)
require(AnchoredInversionUtils)

options(error = utils::recover)

my.seed <- sample(1000, 1)
#my.seed <- 115
    #Tough cases: 721, 115
set.seed(my.seed)


#=================================================================
# Prepare synthetic field and parameterize the field
# by geostatistical formulations.
# Synthetic field is used to check the the result;
# in applications we do not have such 'synthetic' (ie real) field.
# Geostat parameterization is reflected in 'corr.args'
# defined in 'util.R', and in 'mygrid'.
# Another important element is to define how the actual field
# should be transformed to make it take value in (-infty, infty).
# A common transformation is log.
#---------------

K.xmin <- 1e-10
    # Lower boundary of possible value.
    # Specify this lower boundary to avoid generating essentially
    # '0' conductivity, which will fail the Darcy equation.
K.my.mean <- 1e-5
    # Mean value of synthetic field.
    # Choose a value that makes the synthetic data
    # reasonable representation of hydraulic conductivity.
    # Unit is m/s.

# data(volcano)
# myfield <- volcano[, sample(ncol(volcano), 1)]
#     # Vector of length 87.
# myfield <- average.line(myfield, 2) $y

data(Denali)
i <- sample(length(Denali$data) - 400, 1)
myfield <- Denali$data[i : (i+400)]
myfield <- average.line(myfield, 4) $y

lb <- K.xmin + runif(1) * (K.my.mean - K.xmin)/2
    # Min value of synthetic data.
    # Transform 'myfield' to
    #   y = lb + (myfield - min(myfield) * b)
    # and require
    #   mean(y) = K.my.mean
    # then
    #   b = (K.my.mean - lb) / (mean(myfield) - min(myfield))
b <- (K.my.mean - lb) / (mean(myfield) - min(myfield))
myfield <- lb + (myfield - min(myfield)) * b


f <- function(lower)
{
    force(lower)
    function(x, reverse = FALSE)
        AnchoredInversionUtils::log.transform(x, reverse = reverse, lower = lower)
}

f.field.transform <- f(lower = K.xmin)


myfield <- f.field.transform(myfield)


mygrid <- list(
    from = .5,
    by = 1,
    len = length(myfield)
    )


#====================================================
# Prepare linear grid data, that is, linear functions
# of the field on the numerical grid, such as direct
# measurements. Such data may or may not be available.
#----------------

n.linear <- 1
linear.data <- lapply(sample(length(myfield), n.linear),
    function(idx) list(points = idx, value = myfield[idx]))


#=================================================================
# Define forward function and prepare forward data.
# The forward function takes a random field, or a list thereof,
# return corresponding outputs.
# An important element is how to transform the physical outcome
# onto (-infty, infty). This transformation is regarded as part
# of the forward function.
# This transformation should be performed on the actual measurement
# as well.
#
# In applications, this step typically will require calling a non-R
# numerical code. It may be necessary to use intermediate files
# for data storage and transfer. Fixed initial and boundary conditions
# also need to be handled.
#--------------------------

f.darcy <- function(x,
    dx,
    b0.type = 'dirichlet', b0 = 1,
    b1.type = 'dirichlet', b1 = 0,
    ...)
{
    tryCatch(darcy(K = x, dx = dx,
        b0.type = b0.type, b0 = b0,
        b1.type = b1.type, b1 = b1, ...),
        error = function(...) rep(NA, length(x)))
}

darcy.args <- list(dx = mygrid$by)


f.forward.transform <- function(x, ...)
    AnchoredInversionUtils::logit.transform(x, lower = -.01, upper = 1.01, ...)


f.forward <- function(
    x,
    f.darcy,
    darcy.args,
    f.field.transform,
    f.forward.transform,
    forward.idx,
    parallel = AI.getOptions('parallel'))
{
    x.is.list <- is.list(x)
    if (!x.is.list)
        x <- list(x)

    x <- lapply(x, f.field.transform, rev = TRUE)

    z.result <- lapply(x,
        FUN = function(x, darcy.args) {
                do.call(f.darcy, c(list(x), darcy.args))
            },
        darcy.args = c(darcy.args, list(return.idx = forward.idx))
        )

    y <- lapply(z.result, function(z) {
        if (any(is.na(z))) rep(NA, length(z))
        else {
            z <- unname(f.forward.transform(z))
            if (any(!is.finite(z)))
                # If this happens, the input field
                # must be crazy and the darcy output
                # is unreasonable.
                rep(NA, length(z))
            else
                z
            }
        })

    if (x.is.list) y else y[[1]]
}

n.forward <- 12
forward.data.idx <- seq(from = 1, to = mygrid$len, len = n.forward + 2)
forward.data.idx <- round(forward.data.idx[2 : (n.forward + 1)])
    # Uniform sampling
#forward.data.idx <- sample(3 : (mygrid$len - 2), n.forward, replace = FALSE)
    # Random sampling

forward.args <- list(
    f.darcy = f.darcy,
    darcy.args = darcy.args,
    f.field.transform = f.field.transform,
    f.forward.transform = f.forward.transform,
    forward.idx = forward.data.idx)


forward.fun <- function(x) { do.call(f.forward, c(x, forward.args)) }


forward.data <- forward.fun(list(myfield))
    # Forward data.
forward.data.groups <- NULL


#--- plot the actual field ---

trellis.device(color = FALSE, width = 5, height = 4)
print(xyplot(
    f.field.transform(myfield, rev = TRUE)
        ~ AnchoredInversionUtils::grid.fullxyz(mygrid),
    type = 'l',
    xlab = 'X',
    ylab = 'Conductivity')
    )

trellis.device(color = FALSE, width = 5, height = 4)
print(xyplot(
    myfield ~ AnchoredInversionUtils::grid.fullxyz(mygrid),
    type = 'l',
    xlab = 'X',
    ylab = 'Log conductivity')
    )

#--- plot the forward process and measurements ---

forward.curve <- list(
    x = AnchoredInversionUtils::grid.fullxyz(mygrid),
    y = do.call(f.darcy,
        c(list(f.field.transform(myfield, rev = TRUE)),
          darcy.args))
    )
forward.data.x <- grid.ijk2xyz(mygrid, forward.data.idx)
forward.data.y <- f.forward.transform(forward.data, rev = TRUE)

trellis.device(color = FALSE, width = 4, height = 4)
print(xyplot(
    forward.curve$y ~ forward.curve$x,
    panel = function(...) {
            panel.xyplot(...)
            lpoints(forward.data.x, forward.data.y, pch = 20)
        },
    type = 'l',
    lty = 'dashed',
    xlab = 'X (spatial coord)',
    ylab = 'Z (head)'
    ))


#------------------------------------
# Model inference, summary, plotting.
#------------------------------------

save.file.prefix <- 'darcy'

source('anchor.common.R')

trellis.device(color = TRUE, width = 6, height = 7)
z <- plot.anchorit.forward.sample(ai.out, x.at = forward.data.x,
    xlab = 'X coord')
print(z)

