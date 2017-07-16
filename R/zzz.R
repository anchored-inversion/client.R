#' @importFrom futile.logger flog.debug flog.info flog.warn flog.error
.onLoad <- function(libname, pkgname)
{
    options(error = utils::recover)
}



.onUnload <- function(libpath)
{
    library.dynam.unload('AnchoredInversionClient', libpath)
}
