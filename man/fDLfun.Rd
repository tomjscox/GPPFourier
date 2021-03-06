% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SunsetSunrise.R
\name{fDLfun}
\alias{fDLfun}
\title{Relative amount of daylight hours at a specified date and location.}
\usage{
fDLfun(date = "2016-07-01", phi = 51.176, lambda = 4.326, H = 0,
  method = c("insol", "almanac"))
}
\arguments{
\item{date}{POSIXct object or date specificied in unambiguous format. See \code{\link{as.POSIXct}}}

\item{phi}{Latitude}

\item{lambda}{Longitude}

\item{H}{Height of location where fDL is to be calculated}

\item{method}{Method for daylength calculation. "insol" (default) makes use of the function \code{daylength()} from the package \code{insol}; "almanac" makes use of \code{\link{SunRiseSet}} implementing the algorithm from the Astronomical Almanac}
}
\value{
The fraction of daylight hours at the specified date and location. Sunrise and Sunset are calculated with \code{\link{SunRiseSet}}.
}
\description{
Relative amount of daylight hours at a specified date and location.
}
\examples{
fDLfun("2016-06-21")
}
\seealso{
\code{\link{SunRiseSet}}
}
\author{
Tom Cox <tom.cox@uantwerp.be>
}
