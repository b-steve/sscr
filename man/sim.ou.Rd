% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.r
\name{sim.ou}
\alias{sim.ou}
\title{Simulating an Ornstein-Uhlenbeck process}
\usage{
sim.ou(mu, tau, sigma, n.steps, start = NULL)
}
\arguments{
\item{mu}{Point of attraction (conceptually the home-range centre).}

\item{tau}{Temporal autocorrelation parameter.}

\item{sigma}{Spatial autocorrelation parameter.}

\item{n.steps}{Total number of time steps.}

\item{start}{Coordinates of starting location. Default is to
generate a location from the long-term stationary distribution.}
}
\description{
Simulates animal movement from an Ornstein-Uhlenbeck process.
}
\author{
This is a modified version of a function provided by Theo
    Michelot.
}
