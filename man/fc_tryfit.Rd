% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fc_tryfit.R
\name{fc_tryfit}
\alias{fc_tryfit}
\title{Error handling for fitting failCompare models}
\usage{
fc_tryfit(
  y = y,
  y_sfrac = NULL,
  fit_call,
  model = "weibull3",
  non_cen = NULL,
  Hess = NULL,
  inits = NULL,
  ...
)
}
\arguments{
\item{y}{numeric time argument of failure times carried through}

\item{y_sfrac}{survival fraction}

\item{fit_call}{call to dependent model fitting functions.}

\item{model}{model argument passed from fc_fit()}

\item{non_cen}{logical indicating censored variables for use by flexsurv and vitality models}

\item{Hess}{logical argument to fc_fit() carried through}

\item{inits}{initial parameter values for optimization of a single model}

\item{...}{additional arguments passed to dependent functions}
}
\value{
model fitting output for internal use by fc_fit
}
\description{
Error handling for fitting failCompare models
}
\details{
Prevents errors from interrupting single- and multi-model runs using fc_fit
}
