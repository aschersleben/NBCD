% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{getDataRotHyp2}
\alias{getDataRotHyp2}
\title{Generates data of "rotating hyperplane" example.}
\usage{
getDataRotHyp2(angle, n.obs = 1, time)
}
\arguments{
\item{angle}{\code{[numeric(1)]}\cr
Rotation angle of the hyperplane in relation to the x axis (clockwise).}

\item{n.obs}{\code{[numeric(1)]}\cr
Number of observations to draw from the "concept" (distribution).}

\item{time}{\code{[numeric]}\cr
The time stamp to attach to the data.
If missing, it's the "angle" argument.}
}
\value{
\code{[list]}\cr
}
\description{
Generates data of "rotating hyperplane" example.
}
\examples{
set.seed(1909)
dat <- getDataRotHyp2(45, 100)
plot(dat$x, col = dat$class)

}

