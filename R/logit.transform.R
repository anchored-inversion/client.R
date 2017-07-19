## Logit transformation:
## transform between (-Inf, Inf) and (a, b)
## and take care of extreme values.
##
## Forward:  (a, b) --> (-Inf, Inf)
## Reverse:  (-Inf, Inf) --> (a, b)

## A related curve is the S curve (actually the inverse of logit):
##  y = 1 / (1 + exp(-x))

#' Logit transformation.
#'
#' Transform a real number on a interval to the whole real line,
#' or reversely.
#'
#' @param x
#'      Numerical.
#' @param lower
#'      Vector.
#' @param upper
#'      Vector.
#'      \code{lower} and \code{upper} must be numericals and such that
#'      \code{lower < upper}.
#'      This is not checked within this function.
#' @param reverse
#'      Reverse transform?
#'
#' @return
#'      If \code{reverse} is \code{TRUE}, return
#'      \code{(lower + upper * exp(x)) / (1 + exp(x))}.
#'      If \code{reverse} is \code{FALSE} (the default), return
#'      \code{log((x - lower) / (upper - x))}.
#'      The return could contain finite, infinite, or \code{NA} values.
#'
#' @export
logit.transform <- function(
    x,
    lower = 0,
    upper = 1,
    reverse = FALSE)
{
    if (is.vector(x))
    {
        stopifnot(
            length(lower) == 1L || length(lower) == length(x),
            length(upper) == 1L || length(upper) == length(x))
    } else if (is.matrix(x))
    {
        # Each row of 'x' is a 'case' or 'observation'.
        if (length(lower) > 1L)
        {
            stopifnot(length(lower) == ncol(x))
            lower <- matrix(rep(lower, each = nrow(x)), nrow = nrow(x))
        }
        if (length(upper) > 1L)
        {
            stopifnot(length(upper) == ncol(x))
            upper <- matrix(rep(upper, each = nrow(x)), nrow = nrow(x))
        }
    } else
        stop('"x" must be a vector or scalar')


    if (reverse)
    {
        y <- exp(x)
        z <- (lower + upper * y) / (1 + y)
    } else
    {
        idx <- which(is.finite(x) & (x <= lower | x >= upper))
        if (length(idx))
        {
            warning('some inputs are outside of expected limits')
            x[idx] <- NA
        }

        y <- (x - lower) / (upper - x)
        z <- log(y)
    }

    z   # This could contain finite, infinite, or NA values.
}

