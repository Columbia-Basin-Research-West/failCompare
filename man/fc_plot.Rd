% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fc_plot.R
\name{fc_plot}
\alias{fc_plot}
\title{Plotting failure time and sample survival function}
\usage{
fc_plot(
  time,
  surv,
  censorID,
  group = NULL,
  hist = T,
  surv_curv = T,
  main,
  ylim,
  xlim,
  ylab,
  xlab,
  ...
)
}
\arguments{
\item{time}{failure time (x axis)}

\item{surv}{survival function (y axis)}

\item{censorID}{binary or logical variable the same length as \code{time} indicating censored observations, with zeros or FALSE indicating a censored observation}

\item{group}{grouping variable, limit of 3}

\item{hist}{logical show histogram of failure times}

\item{surv_curv}{logical show histogram of failure times}

\item{main}{title for scatterplot}

\item{ylim}{y axis limits for survival plot, used to override default of c(0,1)}

\item{xlim}{x axis limits for survival plot}

\item{ylab}{y axis limits for survival plot, used to override default of c(0,1)}

\item{xlab}{x axis limits for survival plot}

\item{...}{additional argument passed to plot()}
}
\value{
histogram of failure times and/or scatter plot of sample survival function.
}
\description{
Plotting failure time and sample survival function
}
\details{
plot of type "data" shown by default. For "residual" type plot showing (Kaplan-Meier estimates - parametric model fit), plot a singular model of class = fc_obj.
}
