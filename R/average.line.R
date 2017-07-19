#' Sum over a segment of length \code{Q} along a vector (line).
#' Resultant length of data vector is roughly \code{1/Q} of the original.
#'
#' @seealso \code{average.image} in package \pkg{fields}.
#'
#' @export
average.line <- function(obj, Q = 2)
{
    if (is.matrix(obj)) {
        if (ncol(obj) == 2)
            obj <- list(x = obj[, 1], y = obj[, 2])
        else if (ncol(obj) == 1)
            obj <- list(x = seq_along(obj), y <- c(obj))
        else
            stop('wrong input format')
    } else if (is.atomic(obj))
        obj <- list(x = seq_along(obj), y = obj)
    else
        stopifnot(is.list(obj))
            # With two vector members 'x' and 'y'
            # of the same length. Not checked.

    M <- length(obj$x)
    MI <- trunc(M/Q)
    indQ <- 1 : Q

    x <- rep(NA, MI)
    y <- rep(NA, MI)

    for (j in 1:MI) {
        idx <- indQ + (j-1) * Q
        x[j] <- mean(obj$x[idx])
        y[j] <- mean(obj$y[idx], na.rm = TRUE)
    }

    list(x = x, y = y, Q = Q)
}


