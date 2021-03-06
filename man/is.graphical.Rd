% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is.graphical.R
\name{is.graphical}
\alias{is.graphical}
\title{Determines if Network is Graphical}
\usage{
is.graphical(A)
}
\arguments{
\item{A}{A partial correlation network (adjacency matrix)}
}
\value{
Returns a TRUE/FALSE for whether network is graphical
}
\description{
Tests for whether the network is graphical.
Input must be a partial correlation network.
Function assumes that partial correlations were computed from a multivariate normal distribution
}
\examples{
\dontrun{
A <- LoGo(neoOpen, normal = TRUE, partial = TRUE)

is.graphical(A)
}

}
\author{
Alexander Christensen <alexpaulchristensen@gmail.com>
}
