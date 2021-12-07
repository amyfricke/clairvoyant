#' Function to find a Box Cox transformation on a series with trend and
#'  (possibly) seasonality
#'
#' @param history dataframe of dts and actuals to find a Box Cox transform for
#' @param seasonality boolean indicating if seasonality should be applied
#' @param periods integers giving the length of the seasonal periods
#' @return The paramater controlling the Box Cox transformation found
#' @noRd
GetBoxCoxTransform <- function(history, seasonality=TRUE,
                               periods=c(7, 364)) {
  # Find the Box Cox transformation for a time series
  history$dt.idx <- c(1:nrow(history))
  if (seasonality && any(periods > 1)) {
    for (period in periods) {
      history[[paste0('season.idx', period)]] <-
        as.factor(history$dt.idx %% period)
    }
    bxcx.formula <- paste0('actual ~ dt.idx +',
                           paste0('season.idx', periods,
                                  collapse = ' + '))

    bxcx <- MASS::boxcox(formula(bxcx.formula),
                         lambda=seq(0, 1, by=.1), plotit=F, data=history)
  } else {
    bxcx <- MASS::boxcox(actual ~ dt.idx, lambda=seq(0, 1, by=.1),
                         plotit=F, data=history)
  }
  box.cox.lambda <- bxcx$x[which.max(bxcx$y)]
  return(box.cox.lambda)
}

#' Function to aggregate and transform a smoothed weekly time series
#'
#' @param history: smoothed history to transform acutals columns
#' @param cols.transform: columns to transform
#' @param transform: transformation to apply
#' @param period: seasonality period
#' @return transformed + aggregated data frame and aggregated dataframe
#' @noRd
.AggregateAndTransform <- function(history,
                                   periods.agg=NULL,
                                   agg.fun='sum',
                                   periods=c(7, 364),
                                   cols.transform=c('actual',
                                                    'actual.lower',
                                                    'actual.upper'),
                                   dt.format='.AsPOSIXlt',
                                   transform='box_cox') {

  if (length(periods.agg) > 0 &&
      all(periods.agg > 1)) {
    aggregated.list <- AggregateToLongest(history,
                                          periods.agg=periods.agg,
                                          agg.fun=agg.fun,
                                          cols.agg=c('actual', 'actual.lower',
                                                     'actual.upper'),
                                          dt.format=dt.format)
    pre.transform <- aggregated.list[[paste(max(periods.agg))]]
    periods <- periods[!(periods %in% periods.agg)]
    periods <- periods / (max(periods.agg))

  } else {
    aggregated.list <- NULL
    pre.transform <- history
    periods <- periods
  }
  history.list <- list()
  history.list$aggregated <- aggregated.list
  transformed <- pre.transform
  history.list$periods <- periods

  if (transform == 'box_cox') {
    box.cox.lambda <- GetBoxCoxTransform(pre.transform,
                                         seasonality=any(periods > 1),
                                         periods=periods)
  } else {
    box.cox.lambda <- 1
  }
  history.list$box.cox.lambda <- box.cox.lambda
  for (i in 1:length(cols.transform)) {
    transformed[cols.transform[i]] <-
      .TransformTimeSeries(pre.transform[cols.transform[i]],
                           transform=transform,
                           box.cox.lambda=box.cox.lambda)
  }
  history.list$transformed <- transformed
  return(history.list)
}

#' Function to back transform a smoothed weekly forecast
#'
#' @param history: smoothed history to transform acutals columns
#' @param cols.transform: columns to transform
#' @param transform: transformation to apply
#' @param box.cox.lambda: box cox parameter for transformation
#' @return transformed data frame
#' @noRd
.BackTransformForecast <- function(forecast.transformed,
                                   cols.transform=c('forecast',
                                                    'forecast.lower',
                                                    'forecast.upper'),
                                   transform='box_cox',
                                   box.cox.lambda=1) {
  forecast.backtransformed <- forecast.transformed
  for (i in 1:length(cols.transform)) {
    forecast.backtransformed[cols.transform[i]] <-
      .BackTransformTimeSeries(forecast.transformed[cols.transform[i]],
                               transform=transform,
                               box.cox.lambda=box.cox.lambda)
  }
  return(list(forecast=forecast.backtransformed))
}

#' Function to obtain holiday features as dummy variables
#'  to use as external regressors
#'
#' @param history dataframe of dts to get holiday features for
#' @param holidays.df dataframe of dts and holidays to get features for
#' @return matrix of dummy variables for each holiday to be used as external
#'   regressors
#' @export
FeaturizeHolidays <- function(history, holidays.df) {

  holidays.matrix <- subset(history, select=dt)
  unique.holidays <- unique(holidays.df$holiday)
  for (hday in unique.holidays) {

    holiday.dts <- subset(holidays.df, holiday==hday)
    holidays.matrix[hday] <-
      (holidays.matrix$dt %in% holiday.dts$dt)

  }
  return(holidays.matrix)
}

#' Function to obtain holiday features as dummy variables
#'  to use as external regressors (and aggregate as necessary)
#'
#' @param history.start history start dt to get holiday features for
#' @param fcst.dts list of the start and end dts of the forecast
#' @param periods.agg vector of periods that the history is aggregated over
#' @param holidays.df dataframe of dts and holidays to get features for
#' @return matrix of dummy variables for each holiday to be used as external
#'   regressors for both the past and future
#' @export
GetHolidayFeatures <- function(history.start, fcst.dts,
                               dt.units='days', dt.format='.AsPOSIXlt',
                               periods.agg, holidays.df) {
  dt.unit <- substr(dt.units, 1, (nchar(dt.units) - 1))
  total.dts <- data.frame(dt=seq(history.start,
                                 fcst.dts$end.dt, by=paste(1, dt.units)))
  holidays.matrix <- FeaturizeHolidays(total.dts, holidays.df)
  holidays.all <- unique(holidays.df$holiday)
  if (length(periods.agg) > 0 && max(periods.agg) > 1) {
    holidays.matrix <- AggregateToLongest(holidays.matrix,
                                          periods.agg,
                                          agg.fun='max',
                                          cols.agg=holidays.all,
                                          dt.format=dt.format)
    holidays.matrix <- holidays.matrix[[paste(max(periods.agg))]]
  }

  holidays.past <-
    subset(holidays.matrix, dt < fcst.dts$begin.dt, select=-c(dt))
  holidays.future <-
    subset(holidays.matrix, dt >= fcst.dts$begin.dt, select=-c(dt))

  return(list(holidays.past=holidays.past, holidays.future=holidays.future))

}

#' Function to obtain trigonometric seasonality features to use as
#'  external regressors (in place of built in dummy variable seasonality)
#'
#' @param history.len: length of the history that the trig seasonality features
#'   are for
#' @param fcst.len: length of the forecast that the trig seasonality features
#'  are for
#' @param periods: the vector of seasonal periods lengths to create
#'  trigonometric seasonality features of
#' @return matrix of trigonometric seasonality features for each
#'  seasonal period to be used as external regressors
#' @export
GetTrigSeasonalityFeatures <- function(history.len, fcst.len,
                                       periods) {

  overall.index <- c(1:(history.len + fcst.len))
  trig.seasonality.matrix <- data.frame(index=overall.index)
  for (period in periods) {
    period.index <- overall.index %% period
    trig.seasonality.matrix[paste0('sin.period', period)] <-
      sin(2 * pi / period * period.index)
    trig.seasonality.matrix[paste0('cos.period', period)] <-
      cos(2 * pi / period * period.index)
  }

  trig.seasonality.matrix.past <-
    subset(trig.seasonality.matrix, index <= history.len, select=-c(index))
  trig.seasonality.matrix.future <-
    subset(trig.seasonality.matrix, index > history.len, select=-c(index))

  return(list(trig.seasonality.past=trig.seasonality.matrix.past,
              trig.seasonality.future=trig.seasonality.matrix.future))

}
