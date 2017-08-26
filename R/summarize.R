#' Print a summary of an \code{anchorit} object.
#'
#' @param digits passed on to \code{print}.
#' @param \dots additional arguments to \code{print}.
#'
#' @export
print_summary <- function(obj, digits = 3, ...)
{
    cat('Dimensions of model grid:', paste(obj$grid.dim, sep = ' x '), '\n')

    cat('\n')
    cat('  dimension of parameter space:', obj$x.dim, '\n')
    cat('    of which,\n')
    cat('      geostat parameters: ', obj$geost.dim,
        ' (', obj$geost.labels, ')\n')
    cat('      anchors: ', obj$anchor.dim, '\n')
    cat('  number of mixture components retained:', obj$n.mixtures, '\n')

    if (obj$n.models < 2)
    {
        cat('\n')
        cat('Only initial approx is available.\n')
        return(invisible())
    }

    cat('\n')
    cat('Started at     ', obj$date.creation, '\n')
    cat('Last updated at', obj$date.last.update, '\n')
    cat('Number of updates:', obj$n.models - 1, '\n')
    cat('In each update:\n')
    cat('  dimension of conditioning data:', obj$n.forward, '\n')
    if (any(obj$n.forward.reduced != obj$n.forward)) {
        cat('    reduced to:', obj$n.forward.reduced, '\n')
    }

    cat('  reproduction of conditioning data:\n')
    z <- cbind(n = obj$diag$n.samples,
            n.anchors = obj$diag$n.anchors,
            intLL = obj$diag$intLL,
            mad.50 = obj$diag$mad.50,
            mad.75 = obj$diag$mad.75,
            mad.95 = obj$diag$mad.95,
            mad.max = obj$diag$mad.max)
    rownames(z) <- paste('      iter', 1:nrow(z), '  ')
    print(z, digits = digits, ...)
    invisible()
}

