% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ensemble_forecast.R
\name{GetForecastRange}
\alias{GetForecastRange}
\title{Function to get a forecast sensitivty range from an ensemble of forecasts}
\usage{
GetForecastRange(
  forecast.list,
  lower.method = "min",
  upper.method = "max",
  pred.level = 0.8,
  transform = "none",
  box.cox.lambda = 1
)
}
\arguments{
\item{forecast.list:}{list of forecast dataframes (need to contain fields
for dt, forecast, forecast.lower, foreast.upper)}

\item{lower.method:}{function to apply to ensemble of forecasts
(pointwise) to get the lower end of the forecast sensitivity range}

\item{upper.method:}{function to apply to ensemble of forecasts
(pointwise) to get the upper end of the forecast sensitivity range}
}
\value{
: a dataframe containing the range forecast and uncertainty
}
\description{
Function to get a forecast sensitivty range from an ensemble of forecasts
}
