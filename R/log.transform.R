## Log transformation:
## transform betwen (-Inf, Inf) and (lower, Inf)
## and take care of extreme values.
##
## Forward:  (lower, Inf) --> (-Inf, Inf)
## Reverse:  (-Inf, Inf) --> (lower, Inf)

#' Log transformation.
#'
#' Transform a real number on \code{[lower, Inf)}
#' to the whole real line, or reversely.
#'
#' @param x
#'      Numerical.
#' @param lower
#'      Numerical.
#' @param reverse
#'      Reverse transform?
#'
#' @return
#'      If \code{reverse} is \code{TRUE}, return
#'      \code{log(x - lower)}.
#'      If \code{reverse} is \code{FALSE}, return
#'      \code{lower + exp(x)}.
#'
#' @export log.transform
log.transform <- function(x, lower = 0, reverse = FALSE)
{
    if (reverse)
    {
        z <- lower + exp(x)
    } else
    {
        idx <- which(is.finite(x) & (x <= lower))
        if (length(idx))
        {
            warning('some inputs are below the lower bound')
            x[idx] <- NA
        }
        z <- log(x - lower)
    }
    z
}


