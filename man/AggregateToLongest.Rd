% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregation_and_disaggregation.R
\name{AggregateToLongest}
\alias{AggregateToLongest}
\title{Function to Aggregate to the longest seasonal period buckets}
\usage{
AggregateToLongest(
  history,
  periods.agg = c(7),
  agg.fun = "sum",
  cols.agg = c("actual", "actual.lower", "actual.upper"),
  dt.format = ".AsPOSIXlt"
)
}
\arguments{
\item{history:}{Data frame for the unaggregated history.
The required fields are dt and actual.}

\item{periods.agg:}{the shorter seasonal periods to aggregate over}

\item{agg.fun:}{the aggregation function}

\item{cols.agg:}{vector of columns to aggregate}
}
\value{
: weekly.history: Data frame for weekly history
}
\description{
Function to Aggregate to the longest seasonal period buckets
}
