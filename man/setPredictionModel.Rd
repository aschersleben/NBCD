% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NBCD.R
\name{setPredictionModel}
\alias{setPredictionModel}
\title{Create a prediction model from parts of an NBCD model.}
\usage{
setPredictionModel(old.model, new.model, var.name, n.models = Inf)
}
\arguments{
\item{old.model}{\code{[list]}\cr
"old" part of NBCD model.}

\item{new.model}{\code{[list]}\cr
"current" part of NBCD model.}

\item{var.name}{\code{[character]}\cr
Name of variable for which to create prediction model.}

\item{n.models}{\code{[numeric(1)]}\cr
Number of past model values for the linear prediction.}
}
\value{
\code{[list]}\cr
  To be used as "pred.mod" part of "old" part of NBCD model.
  For functions makeNBCDmodel(...) and predict.NBCD(...).
}
\description{
Create a prediction model from parts of an NBCD model.
}

