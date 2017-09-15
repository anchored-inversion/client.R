.onLoad <- function(libname, pkgname) {
    old <- options(error = utils::recover)
    on.exit(options(old), add = TRUE)
    invisible()
}

