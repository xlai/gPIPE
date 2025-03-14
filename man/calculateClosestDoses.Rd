% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/strategies-admissible-Closest.R
\name{calculateClosestDoses}
\alias{calculateClosestDoses}
\title{Calculate Doses Closest to the Boundary}
\usage{
calculateClosestDoses(vec, nrows, ncols)
}
\arguments{
\item{vec}{Binary vector representing the current dose configuration}

\item{nrows}{Number of rows in the dose matrix}

\item{ncols}{Number of columns in the dose matrix}
}
\value{
Binary matrix where TRUE indicates doses closest to the boundary
}
\description{
This function identifies dose combinations that lie on the boundary between
acceptable and unacceptable doses.
}
\keyword{internal}
