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
#' @param b0_type type of upstream boundary condition.
#' @param b0 value of upstream boundary condition.
#'      If \code{b0_type} is \code{'dirichlet'},
#'      \code{b0} provides a head value (m).
#'      If \code{b0_type} is \code{'neumann'},
#'      \code{b0} provides water flux at boundary.
#'      Flux is volumetric rate per volume; same unit as
#'      \code{source}.
#' @param b1_type type of downstream boundary condition.
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
#'
darcy <- function(
    K,
    dx,
    source = rep(0, length(K)),
    b0_type = c('dirichlet', 'neumann'),
    b0,
    b1_type = c('dirichlet', 'neumann'),
    b1
    )
{
    b0_type <- match.arg(b0_type)
    b1_type <- match.arg(b1_type)

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

    if (b0_type == 'dirichlet')
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

    if (b1_type == 'dirichlet')
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


# Lower boundary of possible value.
# Specify this lower boundary to avoid generating essentially
# '0' conductivity, which will fail the Darcy equation.
# This is not only a concern in synthetic studies.
# In real applications, this is a config value used in
# the field transform function.
K_xmin <- 1e-10


darcy_forward_transform <- function(x, ...) {
    logit_transform(x, lower = -.01, upper = 1.01, ...)
}


# Another important element is to define how the actual field
# should be transformed to make it take value in (-infty, infty).
# A common transformation is log.
#
#' @export
darcy_field_transform <- function(x, reverse = FALSE) {
    log_transform(x, reverse = reverse, lower = K_xmin)
}


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
#
# The following function creates the forward function.
#
#' @export
make_darcy_forward_function <- function(mygrid, forward_data_idx = NULL)
{
    f_darcy <- function(
        x,
        b0_type = 'dirichlet',
        b0 = 1,
        b1_type = 'dirichlet',
        b1 = 0,
        ...)
    {
        tryCatch(
            darcy(
                K = x,
                dx = mygrid$by,
                b0_type = b0_type,
                b0 = b0,
                b1_type = b1_type,
                b1 = b1,
                ...),
        error = function(...) rep(NA, length(x))
        )
    }

    # This is the forward function with some config
    # stored in closure.
    function(x) {
        x <- darcy_field_transform(x, rev = TRUE)

        z <- f_darcy(x)
        if (!is.null(forward_data_idx)) {
            z <- z[forward_data_idx]
        }

        if (any(is.na(z)))
        {
            rep(NA, length(z))
        } else
        {
            z <- unname(darcy_forward_transform(z))
            if (any(!is.finite(z)))
                # If this happens, the input field
                # must be crazy and the darcy output
                # is unreasonable.
                rep(NA, length(z))
            else
                z
        }
    }
}


# Prepare synthetic field and parameterize the field
# by geostatistical formulations.
# Synthetic field is used to check the result;
# in applications we do not have such 'synthetic' (ie real) field.
# Geostat parameterization is reflected in 'corr_args'
# and in 'mygrid'.
#
#' @export
make_darcy_field <- function()
{
    K_my_mean <- 1e-5
        # Mean value of synthetic field.
        # Choose a value that makes the synthetic data
        # a reasonable representation of hydraulic conductivity.
        # Unit is m/s.

    data(Denali)
    i <- sample(length(Denali$data) - 400, 1)
    myfield <- Denali$data[i : (i+400)]
    myfield <- average_line(myfield, 4) $y

    lb <- K_xmin + runif(1) * (K_my_mean - K_xmin)/2
        # Min value of synthetic data.
        # Transform 'myfield' to
        #   y = lb + (myfield - min(myfield) * b)
        # and require
        #   mean(y) = K_my_mean
        # then
        #   b = (K_my_mean - lb) / (mean(myfield) - min(myfield))
    b <- (K_my_mean - lb) / (mean(myfield) - min(myfield))
    myfield <- lb + (myfield - min(myfield)) * b

    darcy_field_transform(myfield)
}


# Prepare linear grid data, that is, linear functions
# of the field on the numerical grid, such as direct
# measurements. Such data may or may not be available
# in actual applications.
#
#' @export
make_darcy_linear_data <- function(mygrid, myfield, n_linear)
{
    if (n_linear < 1)
    {
        NULL
    } else
    {
        stopifnot(n_linear < mygrid$len / 10)
        lapply(
            sample(length(myfield), n_linear),
            function(idx) list(points = idx, value = myfield[idx]))
    }
}

