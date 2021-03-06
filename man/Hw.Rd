% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FourierFunctions.R
\name{Hw}
\alias{Hw}
\title{Frequency response of moving average filter}
\usage{
Hw(w, L)
}
\arguments{
\item{w}{The frequency}

\item{L}{Length of the moving average filter}
}
\value{
Hw(w,L) = 1/L*(1-exp(-1i*L*w))/(1-exp(-1i*w))
}
\description{
Frequency response of moving average filter
}
\author{
Tom Cox <tom.cox@uantwerp.be>
}
\keyword{utilities}
