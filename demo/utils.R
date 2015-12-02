# This function is taken from the 'fields' package version 6.6.2.
average.image <- function(obj, Q = 2) {
    # fast method to sum over a QXQ block in image.
    # Q is the number of elements to average over in each dimension
    # e.g.  Q=5 --  blocks of 25 values are averaged to one grid cell.
    if (is.matrix(obj)) {
        obj <- list(x = 1:nrow(obj), y = 1:ncol(obj), z = obj)
    }
    M <- length(obj$x)
    N <- length(obj$y)
    Mi <- trunc(M/Q)
    Ni <- trunc(N/Q)
    # space to hold results
    z <- matrix(NA, nrow = Mi, ncol = N)
    x2 <- rep(NA, Mi)
    y2 <- rep(NA, Ni)
    indQ <- 1:Q
    # sum over block of rows and handle x grid values
    for (j in 1:Mi) {
        x2[j] <- mean(obj$x[indQ + (j - 1) * Q])
        z[j, ] <- colMeans(obj$z[indQ + (j - 1) * Q, ], na.rm = TRUE)
    }
    # sum over blocks of columns  and average y grid values
    for (k in 1:Ni) {
        y2[k] <- mean(obj$y[indQ + (k - 1) * Q])
        z[, k] <- rowMeans(z[, indQ + (k - 1) * Q], na.rm = TRUE)
    }
    return(list(x = x2, y = y2, z = z[1:Mi, 1:Ni], Q = Q))
}
