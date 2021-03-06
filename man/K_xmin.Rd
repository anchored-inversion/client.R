% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/darcy.R
\docType{data}
\name{K_xmin}
\alias{K_xmin}
\title{Minimum possible value of hydraulic conductivity}
\format{An object of class \code{numeric} of length 1.}
\usage{
K_xmin
}
\description{
While the hydraulic conductivity (K) is positive,
we do not assume it takes value in \eqn{(0, \infty)},
but rather in \eqn{(\text{K_xmin}, \infty)}.
}
\details{
It is important to specify this lower boundary because
it helps avoid generating conductivity values that are
essentially 0, which would fail the Darcy equation.
This is necessary in both synthetic studies and real applications.

This lower bound is used in function \code{\link{darcy_field_transform}},
which transforms the hydraulic conductivity values between natual and model scales.
}
\seealso{
\code{\link{darcy_field_transform}}
}
\keyword{datasets}
