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
    b1
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

    c(h)
}



#' Set up environments for a 1-D Darcy example application.
#'
#' @export
darcy.1d <- function(
    n.linear.data = 1,
    n.forward = 12,
    seed = NULL
    )
{
    if (is.null(seed))
        seed <- sample(1000, 1)
        #Tough seed cases: 721, 115
    set.seed(seed)
    flog.info('set.seed called with parameter seed=%s', seed)

    #=================================================================
    # Prepare synthetic field and parameterize the field
    # by geostatistical formulations.
    # Synthetic field is used to check the result;
    # in applications we do not have such 'synthetic' (ie real) field.
    # Geostat parameterization is reflected in 'corr.args'
    # and in 'mygrid'.
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
        # a reasonable representation of hydraulic conductivity.
        # Unit is m/s.

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

    f.field.transform <- function(x, reverse = FALSE) {
        log.transform(x, reverse = reverse, lower = K.xmin)
    }

    myfield <- f.field.transform(myfield)

    mygrid <- list(
        from = .5,
        by = 1,
        len = length(myfield)
        )

    #====================================================
    # Prepare linear grid data, that is, linear functions
    # of the field on the numerical grid, such as direct
    # measurements. Such data may or may not be available
    # in actual applications.
    #----------------

    if (n.linear.data < 1)
    {
        linear.data <- NULL
    } else
    {
        stopifnot(n.linear.data < mygrid$len / 10)
        linear.data <- lapply(
            sample(length(myfield), n.linear.data),
            function(idx) list(points = idx, value = myfield[idx]))
    }

    #=================================================================
    # Define forward function and prepare forward data.
    # The forward function takes a random field
    # and returns corresponding output.
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

    f.darcy <- function(
        x,
        b0.type = 'dirichlet',
        b0 = 1,
        b1.type = 'dirichlet',
        b1 = 0,
        ...)
    {
        tryCatch(
            darcy(
                K = x,
                dx = mygrid$by,
                b0.type = b0.type,
                b0 = b0,
                b1.type = b1.type,
                b1 = b1,
                ...),
        error = function(...) rep(NA, length(x))
        )
    }

    f.forward.transform <- function(x, ...)
        logit.transform(x, lower = -.01, upper = 1.01, ...)

    stopifnot(n.forward < mygrid$len / 3)
    forward.data.idx <- seq(from = 1, to = mygrid$len, len = n.forward + 2)
    forward.data.idx <- round(forward.data.idx[2 : (n.forward + 1)])
        # Uniform sampling
    #forward.data.idx <- sample(3 : (mygrid$len - 2), n.forward, replace = FALSE)
        # Random sampling

    f.forward <- function(x)
    {
        x <- f.field.transform(x, rev = TRUE)

        z <- f.darcy(x)
        if (!is.null(forward.data.idx)) {
            z <- z[forward.data.idx]
        }

        if (any(is.na(z)))
        {
            rep(NA, length(z))
        } else
        {
            z <- unname(f.forward.transform(z))
            if (any(!is.finite(z)))
                # If this happens, the input field
                # must be crazy and the darcy output
                # is unreasonable.
                rep(NA, length(z))
            else
                z
        }
    }

    forward.data <- f.forward(myfield)
    forward.data.groups <- NULL

    forward.data.x <- grid.ijk2xyz(mygrid, forward.data.idx)
    forward.data.y <- f.forward.transform(forward.data, rev = TRUE)

    list(
        mygrid = mygrid,
        myfield = myfield,
        f.field.transform = f.field.transform,
        f.darcy = f.darcy,
        f.forward = f.forward,
        f.forward.transform = f.forward.transform,
        forward.data.x = forward.data.x,
        forward.data.y = forward.data.y,
        forward.data = forward.data,
        forward.data.groups = forward.data.groups,
        linear.data = linear.data
        )
}
