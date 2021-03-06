% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/confidence_seq.R
\name{cross_prob_cs_search}
\alias{cross_prob_cs_search}
\title{Grid-searching the upper bound on the boundary crossing probability for confidence sequences.}
\usage{
cross_prob_cs_search(g, nmax, nmin = 1L, m_upper = 1000L)
}
\arguments{
\item{g}{Constant threshold (positive numeric).}

\item{nmax}{Upper bound of the target time interval.}

\item{nmin}{Lower bound of the target time interval. (default = 1L)}

\item{m_upper}{Upper bound on the grid for searching m value. (default = 1e+3L)}
}
\value{
Grid-search result of the upper bound on the boundary crossing probability for confidence sequences.
}
\description{
\code{cross_prob_cs_search} is used to grid-search the constant threshold for the GLR-like confidence sequences with finite target time interval.
}
