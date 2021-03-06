% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphics.R
\name{plot_field_1d}
\alias{plot_field_1d}
\title{Plot a 1D field curve}
\usage{
plot_field_1d(field, grid, xlab = "X", ylab = "Y", type = "l", ...)
}
\arguments{
\item{field}{A vector of the field values.}

\item{grid}{Grid definition. See `grid.R`.}

\item{...}{Additional paramters passed to `lattice::xyplot`.}
}
\description{
Plot a 1D field curve given the field value vector and model grid definition.
}
\details{
One scenario where this function is useful happens in synthetic studies,
where the true field is known (constructed by the user).
In this case one would like to compare the model simulated field realizations,
requested from the server, against the truth.
Note that the simualtions returned by the server are in model scale.
The plots of simulations shown on the website \url{www.anchored-inversion.com}
are in model scale.
To be compared with these plots, the parameter \code{field} to this function
should be in model scale.
}
