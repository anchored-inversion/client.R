#' Steady state 1-D groundwater flow equation, using Darcy's law
#'
#' Solve the 1-D ODE: \code{d(K dh/dx) / dx = s}.
#'
#' @param K hydraulic conductivity (m/s); 1-D vector.
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
#' \code{K} (conductivity), \code{s}, \code{h} (head)
#' are all evaluated at the model grids.
#'
#' @return Steady-state hydraulic head (m) in each grid cell.
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


#' Minimum possible value of hydraulic conductivity
#'
#' While the hydraulic conductivity (K) is positive,
#' we do not assume it takes value in \eqn{(0, \infty)},
#' but rather in \eqn{(\text{K_xmin}, \infty)}.
#'
#' It is important to specify this lower boundary because
#' it helps avoid generating conductivity values that are
#' essentially 0, which would fail the Darcy equation.
#' This is necessary in both synthetic studies and real applications.
#'
#' This lower bound is used in function \code{\link{darcy_field_transform}},
#' which transforms the hydraulic conductivity values between natual and model scales.
#'
#' @seealso \code{\link{darcy_field_transform}}
K_xmin <- 1e-10


#' Transform the forward data vector between natural and model scales
#'
#' Because Anchored Inversion requires the forward data vector to be
#' defined on the real line, i.e. \eqn{(-\infty, \infty)} (the "model scale"),
#' whereas the forward data in the application may reside in a different range
#' (the "natural scale", e.g. certain data elements may be positive by definition),
#' a transformation is needed between these two scales.
#'
#' In the Darcy example application, the forward data are hydraulic heads
#' whose range are known because the value must lie between the two boundary
#' values, which are set by us (we set them to be 0 and 1).
#' We use a logit transformation between the natual and model scales.
#'
#' In a real world application, all the forward data items may not be of the same nature.
#' Consequently, different elements of the forward data vector may need different
#' transformations. That would be an internal detail of the transformation function.
#'
#' This transformation is applied on the forward data (either from actual measurement
#' or from running the forward model with a simulated field), which is on the natural scale,
#' before they are submitted to the Anchored Inversion server for model inference.
#' The inverse transformation can be useful for diagnostics.
#' On possible scenario is this: the forward model
#' may output the transformed forward data directly (ready to be submitted to the server),
#' but we want to reverse-transform them into the natual scales for comparison with
#' the field-measured forward data.
#'
#' @param x Forward data vector.
#' @param reverse If \code{TRUE}, transform from model to natural scale; if \code{FALSE},
#'     transform from natural to model scale.
#'
#' @return Transformed forward data vector.
#' @seealso \code{\link{make_darcy_forward_function}}
darcy_forward_transform <- function(x, reverse = FALSE, ...) {
    logit_transform(x, lower = -.01, upper = 1.01, reverse = reverse, ...)
}


#' Transform the field value vector between natural and model scales
#'
#' Because Anchored Inversion works on a spatial variable that is defined by
#' the entire real line, i.e. \eqn{(-\infty, \infty)},
#' this transformation is needed if the actual field
#' is not defined on the real line. A typical case is that the actual field value
#' is positive. In this case, the log transformation may be suitable.
#' If the actual field value is already define on the real line,
#' then this transformation is not needed.
#'
#' For the Darcy example application, we assume the field value is
#' on \eqn{[\text{K_xmin}, \infty)} (referred to as the "natural scale"),
#' and use a log transformation to map it onto the real line (referred to as the "model scale").
#'
#' Field simulations requested from the Anchored Inversion server are always in the model scale.
#' Before running the forward model with these fields, the fields need to be reverse-transformed
#' into the natual scale. If the field simulations are to be used for other purposes,
#' they need to be reverse-transformed similarly.
#'
#' @param x The field value vector.
#' @param reverse If \code{FALSE} (the default), transform from natural scale to model scale;
#'     If \code{TRUE} (that is, "reverse transformation"), transform from model scale to natural scale.
#' @return Transformed field value vector.
#'
#' @export
darcy_field_transform <- function(x, reverse = FALSE) {
    log_transform(x, reverse = reverse, lower = K_xmin)
}


#' Construct a forward model for the Darcy example application
#'
#' This is a "higher-order function". Calling it with required parameters
#' will get the forward function object in the return.
#'
#' In this example application, the forward
#' model implements the "Darcy's Law" on the input dydraulic conductivity field.
#' The result is the field of "hydraulic head" on the same model grid.
#' We assume that, in the imagined field study,
#' we measure the head values at select locations as our (forward) data.
#' In this synthetic example application, instead of "measure",
#' we simply "pick" the head values at the specified locations from the
#' vector of head values that are provided by running the forward model
#' (i.e. \code{\link{darcy}}) with the synthetic conductivity field.
#'
#' The forward function that is generated takes a single parameter:
#' the conductivity field vector \emph{on the model scale}.
#' This is expected to be a field simulation obtained from the
#' Anchored Inversion server.
#' In this function, the field vector is first inverse-transformed to the natural scale
#' (using \code{\link{darcy_field_transform}}).
#' Then function \code{\link{darcy}} is used to compute the hydraulic head values
#' based on Darcy's Law.
#' Then the specified head values are picked (using \code{\link{forward_data_idx}}) as the forward data.
#' At that point the forward data are on ther \emph{natural scale}.
#' They are then transformed to the \emph{model scale} using \code{\link{darcy_forward_transform}},
#' and returned.
#' If anything goes wrong in this whole process, a vector of \code{NA} is returned.
#' Otherwise, all elements of the vector should be valid, realistic values.
#' A vector of \code{NA} in the return usually suggests that the input field is unrealistic
#' so that Darcy's Law fails.
#' This happens routinely in the early iterations of the model.
#'
#' To recap, the generated forward function takes a field vector on the model scale
#' (straight from simulations by the Anchored Inversion server)
#' and returns a forward data vector on the model scale (ready to be submitted back
#' to the server).
#'
#' Note that the reverse-transformation on the natural-scale forward data
#' This transformation should be performed on the actual measurement
#' as well.
#'
#' In real applications, this forward function typically will require calling a non-R
#' numerical code. It may be necessary to use intermediate files
#' for data storage and transfer. It may use distributed and parallel computing
#' for performance.
#' Fixed initial and boundary conditions also need to be handled.
#'
#' In real applications, the field-measured forward data, which are on the natural scale,
#' need to be reverse-transformed to model scale before submitting to the server,
#' because the Anchored Inversion method expects the forward data to be on the model scale.
#' This is analogous to the last step in the generated forward function.
#'
#' @param grid Model grid definition.
#' @param forward_data_idx Model grid indices where the hydraulic head values
#'     shold be used as forward data. Default is \code{NULL}, meaning
#'     all the model grid points.
#' @param Generated forward function object.
#'
#' @seealso \code{\link{darcy_field_transform}}, \code{\link{darcy_forward_transform}}
#' @export
make_darcy_forward_function <- function(grid, forward_data_idx = NULL)
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
                dx = grid$by,
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
        # The input `x` is on the model scale.
        # Transform it to the natual scale before
        # applying the Darcy's Law.
        x <- darcy_field_transform(x, reverse = TRUE)

        z <- f_darcy(x)
        if (!is.null(forward_data_idx)) {
            z <- z[forward_data_idx]
        }

        if (any(is.na(z)))
        {
            rep(NA, length(z))
        } else
        {
            # The forward data returned from `f_darcy`
            # are in their natural scale.
            # Transform them to the model scale
            # so that they are ready to be submitted
            # to the Anchored Inversion server.
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
# and in 'grid'.
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
make_darcy_linear_data <- function(grid, myfield, n_linear)
{
    if (n_linear < 1)
    {
        NULL
    } else
    {
        stopifnot(n_linear < grid$len / 10)
        lapply(
            sample(length(myfield), n_linear),
            function(idx) list(points = idx, value = myfield[idx]))
    }
}

