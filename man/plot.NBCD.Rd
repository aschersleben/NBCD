% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NBCD.R
\name{plot.NBCD}
\alias{plot.NBCD}
\title{Plot predicted class values of NBCD model.}
\usage{
\method{plot}{NBCD}(x, time, use.lm = FALSE, ..., n.models = Inf,
  .verb = TRUE)
}
\arguments{
\item{x}{\code{[NBCD]}\cr
NBCD model.}

\item{time}{\code{[numeric(1)]}\cr
Time for which to predict (passed to predict.NBCD(...)).}

\item{use.lm}{\code{[logical(1)]}\cr
Use lm models to forecast the movements of mean and sd? (Default TRUE)
If FALSE, use mean and sd of old.model.
(passed to predict.NBCD(...))}

\item{...}{\code{[any]}\cr
Arguments passed to plot.nb2(...).}

\item{n.models}{\code{[numeric(1)]}\cr
Number of past model values for the linear prediction.}

\item{.verb}{\code{[logical(1)]}\cr
Show messages? Default TRUE.}
}
\description{
Plot predicted class values of NBCD model.
}

