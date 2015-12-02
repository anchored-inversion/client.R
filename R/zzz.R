.onLoad <- function(libname, pkgname)
{
    invisible()
}



.onUnload <- function(libpath)
    library.dynam.unload('AnchoredInversionExamples', libpath)

