% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fc_pred.R
\name{fc_pred}
\alias{fc_pred}
\title{Failure time predictions}
\usage{
fc_pred(mod_obj = NULL, times, pars = NULL, model = NULL)
}
\arguments{
\item{mod_obj}{model object (class = \code{fc_obj})}

\item{times}{time vector}

\item{pars}{parameter estimates, if mod_obj absent}

\item{model}{survival model name, if mod_obj absent}
}
\value{
numeric vector failure/survival probability
}
\description{
generates predictions from failure time model objects.
}
