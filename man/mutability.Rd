% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mutability.R
\name{mutability}
\alias{mutability}
\title{Mutability characteristic. Searches interval with the maximum growth rate in the \strong{time.window}
and backtraces hypothetical resistant population.}
\usage{
mutability(
  tables,
  tables.gr,
  tables.col,
  tables.cum,
  time.window = c(25, 50),
  tube = "Tube1",
  MIC = 0.004,
  Dtime = 8,
  PPower = 100,
  lt = 0.15,
  V = 20
)
}
\arguments{
\item{tables}{- List of tables generated by \emph{calaculate_tables} function.}

\item{tables.gr}{- List of tables generated by \emph{growth.rates} function.}

\item{tables.cum}{- List of tables generated by \emph{cumulative}}

\item{time.window}{- A vector that specifies interval where a point with maximum growth rate will initiate
backtracing of a hypothetical population.}

\item{MIC}{- Minimum Inhibitory Concentration.}

\item{Dtime}{- Dilution time.}

\item{PPower}{- A power of peristaltic pumps. Default = 100.}

\item{V}{- Volume of of media in the reactors. Default = 20.}
}
\description{
Mutability characteristic. Searches interval with the maximum growth rate in the \strong{time.window}
and backtraces hypothetical resistant population.
}
