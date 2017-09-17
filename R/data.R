#' Denali dataset
#'
#' 1D elevation proile near Mt Denali.
#'
#' This is a 1D, west-east elevation profile south of Mt McKinley
#' (or Denali). It can be used as synthetic 1D data.
#'
#' The profile is south of Mt McKinley.
#' I avoided Mt McKinley because it would create a dominant peak,
#' making the profile less interesting.
#' The elevation data were downloaded from the National Map Seamless
#' Server, \url{http://seamless.usgs.gov},
#' National Elevation Dataset Alaska (NED) 2 Arc Second,
#' in `Grid Float' format, and read in using function \code{gridfloat}.
#'
#' @docType data
#' @name Denali
#' @format A list consisting of elements \code{data} (1D vector of length 2278)
#'     and \code{grid} (the "grid" definition).
#' @source \url{http://seamless.usgs.gov}
Denali <- function() {
    data(Denali)
}


#' HalfDome dataset
#'
#' 2D elevation map near Half Dome, Yosemite National Park.
#'
#' This is a 2D elevation map near Half Dome in Yosemite National Park.
#' It can be used as synthetic 2D data.
#'
#' The elevation data were downloaded from the National Map Seamless
#' Server, \url{http://seamless.usgs.gov}, in `Grid Float' format,
#' National Elevation Dataset (NED) 1 Arc Second,
#' and read in using function \code{gridfloat}.
#'
#' @docType data
#' @name HalfDome
#' @source \url{http://seamless.usgs.gov}
HalfDome <- function() {
    data(HalfDome)
}

