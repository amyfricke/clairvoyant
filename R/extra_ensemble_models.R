#' Option for forecasting (aggregated or disaggregated, long term or short
#' term). Finds the Holt Winters forecast. Cannot incorporate external
#' regressors
#'
#' @param history: dataframe containing fields for dt, actual, actual.lower,
#'   actual.upper for forecasting
#' @param fcst.dts: named vector containing start and end dt of forecast
#' @param dt.units: units for the datetimes (dts) in the history, string
#' @param dt.format: string giving a formatting function for the dts
#' @param periods: length of seasonal period (in aggregated buckets of
#' periods.agg)
#' @param periods.agg: periods that training data (and forecast will be) are
#'   aggregated over
#' @param periods.trig: seasonal periods to fit trigonometric curves to as
#'   external regressors
#' @param pred.level: level for the prediction intervals
#' @param transform: the transformation applied and thus the inverse to apply
#'   to get forecast on the original scale
#' @param box.cox.lambda: the value of the box cox lambda in the event the
#'   transform is 'box_cox'
#' @param x.reg: an optional dataframe of external regressors
#' @param x.future an optional (but required if x.reg is not null) of future
#'    external regressors
#' @param holidays.df a dataframe containing holidays to model as features
#' @return: dataframe containing the dt, forecast, forecast lower prediction
#'   limit and upper prediction limit
HWinters <- function(history, fcst.dts,
                     dt.units='days', dt.format='.AsPOSIXlt',
                     periods=c(7), periods.agg=NULL, periods.trig=c(364),
                     pred.level=0.8, transform='box_cox', box.cox.lambda=1,
                     x.reg=NULL, x.future=NULL, holidays.df=NULL) {
  
  stopifnot(is.null(x.reg) && is.null(x.future) && length(period.trig) == 0 &&
              is.null(holidays.df))
  
  fcst.interval <- max(c(1, periods.agg))
  dt_units <- get(dt.units)
  fcst.dts.seq <- seq((fcst.dts$begin.dt + dt_units(fcst.interval - 1)),
                      fcst.dts$end.dt, by=paste(fcst.interval, dt.units))
  len.fcst <- length(fcst.dts.seq)
  
  ts.training <- ts(history$actual, frequency=max(periods))
  
  seasonality <- (period > 1)
  seasonal <- 'add'
  fcst <- predict(HoltWinters(ts.training,
                              seasonal=seasonal), len.fcst, seasonality,
                  level=pred.level)[, 1]
  
  ts.training.lower <- ts(history$actual.lower, frequency=period)
  hw.fit.lower <- HoltWinters(ts.training.lower, gamma=seasonality,
                              seasonal=seasonal)
  fcst.lower <- predict(hw.fit.lower, len.fcst, T, level=pred.level)[, 3]
  
  ts.training.upper <- ts(history$actual.upper, frequency=period)
  hw.fit.upper <- HoltWinters(ts.training.upper, gamma=seasonality,
                              seasonal=seasonal)
  fcst.upper <- predict(hw.fit.upper, len.fcst, T, level=pred.level)[, 2]
  
  fcst.df <- data.frame(dt=fcst.dts.seq, forecast=fcst,
                        forecast.lower=fcst.lower,
                        forecast.upper=fcst.upper)
  fcst.df <- .BackTransformForecast(fcst.df, transform=transform,
                                    box.cox.lambda=box.cox.lambda)
  return(fcst.df)
  
}

#' Option for forecasting (aggregated or disaggregated, long term or
#' short term). Gets an STL decomposition and then forecasts using ets
#'
#' @param history: dataframe containing fields for dt, actual, actual.lower,
#'   actual.upper for forecasting
#' @param fcst.dts: named vector containing start and end dt of forecast
#' @param dt.units: units for the datetimes (dts) in the history, string
#' @param dt.format: string giving a formatting function for the dts
#' @param periods: length of seasonal period (in aggregated buckets of
#' periods.agg)
#' @param periods.agg: periods that training data (and forecast will be) are
#'   aggregated over
#' @param periods.trig: seasonal periods to fit trigonometric curves to as
#'   external regressors
#' @param pred.level: level for the prediction intervals
#' @param transform: the transformation applied and thus the inverse to apply
#'   to get forecast on the original scale
#' @param box.cox.lambda: the value of the box cox lambda in the event the
#'   transform is 'box_cox'
#' @param x.reg: an optional dataframe of external regressors
#' @param x.future an optional (but required if x.reg is not null) of future
#'    external regressors
#' @param holidays.df a dataframe containing holidays to model as features
#' @return: dataframe containing the dt, forecast, forecast lower prediction
#'   limit and upper prediction limit
StlEts <- function(history, fcst.dts, dt.units='days', dt.format='.AsPOSIXlt',
                   periods=c(7), periods.agg=NULL, periods.trig=c(364),
                   pred.level=0.8, transform='box_cox', box.cox.lambda=1,
                   x.reg=NULL, x.future=NULL, holidays.df=NULL) {
  
  stopifnot(is.null(x.reg) && is.null(x.future) && length(period.trig) == 0 &&
              is.null(holidays.df))
  
  fcst.interval <- max(c(1, periods.agg))
  dt_units <- get(dt.units)
  fcst.dts.seq <- seq((fcst.dts$begin.dt + dt_units(fcst.interval - 1)),
                      fcst.dts$end.dt, by=paste(fcst.interval, dt.units))
  len.fcst <- length(fcst.dts.seq)
  
  ts.training <- ts(history$actual, frequency=max(periods))
  
  stlets.fit <- suppressWarnings(
    forecast::stlm(ts.training, biasadj=T, method='ets'))
  
  fcst <- forecast::forecast(stlets.fit, h=len.fcst, level=pred.level,
                             allow.multiplicative.trend=TRUE)
  
  fcst.df <- data.frame(dt=fcst.dts.seq, forecast=fcst$mean,
                        forecast.lower=as.vector(fcst$lower),
                        forecast.upper=as.vector(fcst$upper))
  
  fcst.df <- .BackTransformForecast(fcst.df, transform=transform,
                                    box.cox.lambda=box.cox.lambda)
  return(fcst.df)
  
}

#' Option for forecasting (aggregated or disaggregated, long term or short
#' term). Gets an STL decomposition and then forecasts using theta
#'
#' @param history: dataframe containing fields for dt, actual, actual.lower,
#'   actual.upper for forecasting
#' @param fcst.dts: named vector containing start and end dt of forecast
#' @param dt.units: units for the datetimes (dts) in the history, string
#' @param dt.format: string giving a formatting function for the dts
#' @param periods: length of seasonal period (in aggregated buckets of
#' periods.agg)
#' @param periods.agg: periods that training data (and forecast will be) are
#'   aggregated over
#' @param periods.trig: seasonal periods to fit trigonometric curves to as
#'   external regressors
#' @param pred.level: level for the prediction intervals
#' @param transform: the transformation applied and thus the inverse to apply
#'   to get forecast on the original scale
#' @param box.cox.lambda: the value of the box cox lambda in the event the
#'   transform is 'box_cox'
#' @param x.reg: an optional dataframe of external regressors
#' @param x.future an optional (but required if x.reg is not null) of future
#'    external regressors
#' @param holidays.df a dataframe containing holidays to model as features
#' @return: dataframe containing the dt, forecast, forecast lower prediction
#'   limit and upper prediction limit
StlTheta <- function(history, fcst.dts,
                     dt.units='days', dt.format='.AsPOSIXlt',
                     periods=c(7), periods.agg=NULL, periods.trig=c(364),
                     pred.level=0.8, transform='box_cox', box.cox.lambda=1,
                     x.reg=NULL, x.future=NULL, holidays.df=NULL) {
  
  stopifnot(is.null(x.reg) && is.null(x.future) && length(period.trig) == 0 &&
              is.null(holidays.df))
  
  fcst.interval <- max(c(1, periods.agg))
  dt_units <- get(dt.units)
  fcst.dts.seq <- seq((fcst.dts$begin.dt + dt_units(fcst.interval - 1)),
                      fcst.dts$end.dt, by=paste(fcst.interval, dt.units))
  len.fcst <- length(fcst.dts.seq)
  
  ts.training <- ts(history$actual, frequency=max(periods))
  stltheta.fit <- suppressWarnings(stl(ts.training, s.window= 7 +  4 * seq(6)))
  
  fcst <- forecast::forecast(stltheta.fit, h=len.fcst, forecastfunction=thetaf,
                             level=pred.level)
  fcst.df <- data.frame(dt=fcst.dts.seq, forecast=fcst$mean,
                        forecast.lower=as.vector(fcst$lower),
                        forecast.upper=as.vector(fcst$upper))
  
  fcst.df <- .BackTransformForecast(fcst.df, transform=transform,
                                    box.cox.lambda=box.cox.lambda)
  return(fcst.df)
  
}

#' Option for forecasting (aggregated or disaggregated, long term or short
#' term). Gets an STL decomposition and then forecasts using thief mlp
#'
#' @param history: dataframe containing fields for dt, actual, actual.lower,
#'   actual.upper for forecasting
#' @param fcst.dts: named vector containing start and end dt of forecast
#' @param dt.units: units for the datetimes (dts) in the history, string
#' @param dt.format: string giving a formatting function for the dts
#' @param periods: length of seasonal period (in aggregated buckets of
#' periods.agg)
#' @param periods.agg: periods that training data (and forecast will be) are
#'   aggregated over
#' @param periods.trig: seasonal periods to fit trigonometric curves to as
#'   external regressors
#' @param pred.level: level for the prediction intervals
#' @param transform: the transformation applied and thus the inverse to apply
#'   to get forecast on the original scale
#' @param box.cox.lambda: the value of the box cox lambda in the event the
#'   transform is 'box_cox'
#' @param x.reg: an optional dataframe of external regressors
#' @param x.future an optional (but required if x.reg is not null) of future
#'    external regressors
#' @param holidays.df a dataframe containing holidays to model as features
#' @return: dataframe containing the dt, forecast, forecast lower prediction
#'   limit and upper prediction limit
StlMlp <- function(history, fcst.dts, dt.units='days', dt.format='.AsPOSIXlt',
                   periods=c(7), periods.agg=NULL, periods.trig=c(364),
                   pred.level=0.8, transform='box_cox', box.cox.lambda=1,
                   x.reg=NULL, x.future=NULL, holidays.df=NULL) {
  
  fcst.interval <- max(c(1, periods.agg))
  dt_units <- get(dt.units)
  fcst.dts.seq <- seq((fcst.dts$begin.dt + dt_units(fcst.interval - 1)),
                      fcst.dts$end.dt, by=paste(fcst.interval, dt.units))
  len.fcst <- length(fcst.dts.seq)
  
  ts.training <- ts(history$actual, frequency=max(periods))
  stl.fit <- suppressWarnings(stl(ts.training, s.window= 7 +  4 * seq(6)))
  
  if ((!is.null(x.reg) && !is.null(x.future)) ||
      (!is.null(holidays.df)) || (length(periods.trig) > 0)) {
    
    if (!is.null(holidays.df)) {
      
      holidays.matrix <- GetHolidayFeatures(history.start=
                                              (min(history$dt) -
                                                 dt_units(fcst.interval - 1)),
                                            fcst.dts=fcst.dts,
                                            dt.units=dt.units,
                                            dt.format=dt.format,
                                            periods.agg=periods.agg,
                                            holidays.df=holidays.df)
      x.reg <- cbind(x.reg, holidays.matrix$holidays.past)
      x.future <- cbind(x.future, holidays.matrix$holidays.future)
      
    }
    
    if (length(periods.trig) > 0) {
      trig.matrix <- GetTrigSeasonalityFeatures(nrow(history), len.fcst,
                                                periods.trig)
      
      x.reg <- cbind(x.reg, trig.matrix$trig.seasonality.past)
      x.future <- cbind(x.future, trig.matrix$trig.seasonality.future)
    }
    
    fcst <- forecast::forecast(stl.fit, h=len.fcst,
                               forecastfunction=nnfor::mlp.thief,
                               level=pred.level, xreg=x.reg, newxreg=x.future)
  } else {
    fcst <- forecast::forecast(stl.fit, h=len.fcst,
                               forecastfunction=nnfor::mlp.thief,
                               level=pred.level)
  }
  pred.level2 <- (1 - pred.level) / 2
  fcst.std <- apply(fcst$all.mean, 1, FUN=sd)
  fcst.lower <- fcst$mean + qnorm(pred.level2) * fcst.std
  fcst.upper <- fcst$mean - qnorm(pred.level2) * fcst.std
  
  fcst.df <- data.frame(dt=fcst.dts.seq, forecast=fcst$mean,
                        forecast.lower=fcst.lower,
                        forecast.upper=fcst.upper)
  fcst.df <- .BackTransformForecast(fcst.df, transform=transform,
                                    box.cox.lambda=box.cox.lambda)
  return(fcst.df)
  
}

#' Option for forecasting (aggregated or disaggregated, long term or short
#' term). Gets an STL decomposition and then forecasts using thief elm
#'
#' @param history: dataframe containing fields for dt, actual, actual.lower,
#'   actual.upper for forecasting
#' @param fcst.dts: named vector containing start and end dt of forecast
#' @param dt.units: units for the datetimes (dts) in the history, string
#' @param dt.format: string giving a formatting function for the dts
#' @param periods: length of seasonal period (in aggregated buckets of
#' periods.agg)
#' @param periods.agg: periods that training data (and forecast will be) are
#'   aggregated over
#' @param periods.trig: seasonal periods to fit trigonometric curves to as
#'   external regressors
#' @param pred.level: level for the prediction intervals
#' @param transform: the transformation applied and thus the inverse to apply
#'   to get forecast on the original scale
#' @param box.cox.lambda: the value of the box cox lambda in the event the
#'   transform is 'box_cox'
#' @param x.reg: an optional dataframe of external regressors
#' @param x.future an optional (but required if x.reg is not null) of future
#'    external regressors
#' @param holidays.df a dataframe containing holidays to model as features
#' @return: dataframe containing the dt, forecast, forecast lower prediction
#'   limit and upper prediction limit
StlElm <- function(history, fcst.dts, dt.units='days', dt.format='.AsPOSIXlt',
                   periods=c(7), periods.agg=NULL, periods.trig=c(364),
                   pred.level=0.8, transform='box_cox', box.cox.lambda=1,
                   x.reg=NULL, x.future=NULL, holidays.df=NULL) {
  
  fcst.interval <- max(c(1, periods.agg))
  dt_units <- get(dt.units)
  fcst.dts.seq <- seq((fcst.dts$begin.dt + dt_units(fcst.interval - 1)),
                      fcst.dts$end.dt, by=paste(fcst.interval, dt.units))
  len.fcst <- length(fcst.dts.seq)
  
  ts.training <- ts(history$actual, frequency=max(periods))
  stl.fit <- suppressWarnings(stl(ts.training, s.window= 7 +  4 * seq(6)))
  
  if ((!is.null(x.reg) && !is.null(x.future)) ||
      (!is.null(holidays.df)) || (length(periods.trig) > 0)) {
    
    if (!is.null(holidays.df)) {
      
      holidays.matrix <- GetHolidayFeatures(history.start=
                                              (min(history$dt) -
                                                 dt_units(fcst.interval - 1)),
                                            fcst.dts=fcst.dts,
                                            dt.units=dt.units,
                                            dt.format=dt.format,
                                            periods.agg=periods.agg,
                                            holidays.df=holidays.df)
      x.reg <- cbind(x.reg, holidays.matrix$holidays.past)
      x.future <- cbind(x.future, holidays.matrix$holidays.future)
      
    }
    
    if (length(periods.trig) > 0) {
      trig.matrix <- GetTrigSeasonalityFeatures(nrow(history), len.fcst,
                                                periods.trig)
      
      x.reg <- cbind(x.reg, trig.matrix$trig.seasonality.past)
      x.future <- cbind(x.future, trig.matrix$trig.seasonality.future)
    }
    
    fcst <- forecast::forecast(stl.fit, h=len.fcst,
                               forecastfunction=nnfor::elm.thief,
                               level=pred.level, xreg=x.reg, newxreg=x.future)
  } else {
    fcst <- forecast::forecast(stl.fit, h=len.fcst,
                               forecastfunction=nnfor::elm.thief,
                               level=pred.level)
  }
  
  pred.level2 <- (1 - pred.level) / 2
  fcst.std <- apply(fcst$all.mean, 1, FUN=sd)
  fcst.lower <- fcst$mean + qnorm(pred.level2) * fcst.std
  fcst.upper <- fcst$mean - qnorm(pred.level2) * fcst.std
  
  fcst.df <- data.frame(dt=fcst.dts.seq, forecast=fcst$mean,
                        forecast.lower=fcst.lower,
                        forecast.upper=fcst.upper)
  fcst.df <- .BackTransformForecast(fcst.df, transform=transform,
                                    box.cox.lambda=box.cox.lambda)
  return(fcst.df)
  
}

#' Option for forecasting (aggregated or disaggregated, long term or short
#' term) using thief mlp
#'
#' @param history: dataframe containing fields for dt, actual, actual.lower,
#'   actual.upper for forecasting
#' @param fcst.dts: named vector containing start and end dt of forecast
#' @param dt.units: units for the datetimes (dts) in the history, string
#' @param dt.format: string giving a formatting function for the dts
#' @param periods: length of seasonal period (in aggregated buckets of
#' periods.agg)
#' @param periods.agg: periods that training data (and forecast will be) are
#'   aggregated over
#' @param periods.trig: seasonal periods to fit trigonometric curves to as
#'   external regressors
#' @param pred.level: level for the prediction intervals
#' @param transform: the transformation applied and thus the inverse to apply
#'   to get forecast on the original scale
#' @param box.cox.lambda: the value of the box cox lambda in the event the
#'   transform is 'box_cox'
#' @param x.reg: an optional dataframe of external regressors
#' @param x.future an optional (but required if x.reg is not null) of future
#'    external regressors
#' @param holidays.df a dataframe containing holidays to model as features
#' @return: dataframe containing the dt, forecast, forecast lower prediction
#'   limit and upper prediction limit
Mlp <- function(history, fcst.dts, dt.units='days', dt.format='.AsPOSIXlt',
                periods=c(7), periods.agg=NULL, periods.trig=c(364),
                pred.level=0.8, transform='box_cox', box.cox.lambda=1,
                x.reg=NULL, x.future=NULL, holidays.df=NULL) {
  
  fcst.interval <- max(c(1, periods.agg))
  dt_units <- get(dt.units)
  fcst.dts.seq <- seq((fcst.dts$begin.dt + dt_units(fcst.interval - 1)),
                      fcst.dts$end.dt, by=paste(fcst.interval, dt.units))
  len.fcst <- length(fcst.dts.seq)
  
  ts.training <- ts(history$actual, frequency=max(c(1, periods)))
  if ((!is.null(x.reg) && !is.null(x.future)) ||
      (!is.null(holidays.df)) || (length(periods.trig) > 0)) {
    
    if (!is.null(holidays.df)) {
      
      holidays.matrix <- GetHolidayFeatures(history.start=
                                              (min(history$dt) -
                                                 dt_units(fcst.interval - 1)),
                                            fcst.dts=fcst.dts,
                                            dt.units=dt.units,
                                            dt.format=dt.format,
                                            periods.agg=periods.agg,
                                            holidays.df=holidays.df)
      x.reg <- cbind(x.reg, holidays.matrix$holidays.past)
      x.future <- cbind(x.future, holidays.matrix$holidays.future)
      
    }
    
    if (length(periods.trig) > 0) {
      trig.matrix <- GetTrigSeasonalityFeatures(nrow(history), len.fcst,
                                                periods.trig)
      
      x.reg <- cbind(x.reg, trig.matrix$trig.seasonality.past)
      x.future <- cbind(x.future, trig.matrix$trig.seasonality.future)
    }
    
    x.all <- rbind(x.reg, x.future)
    mlp.fit <- suppressWarnings(nnfor::mlp(ts.training, xreg=x.all,
                                           lags=0, reps=50))
    fcst <- forecast::forecast(mlp.fit, h=len.fcst,
                               level=pred.level)
  } else {
    mlp.fit <- suppressWarnings(nnfor::mlp(ts.training, reps=50))
    fcst <- forecast::forecast(mlp.fit, h=len.fcst,
                               level=pred.level, xreg=x.all)
  }
  
  pred.level2 <- (1 - pred.level) / 2
  fcst.std <- apply(fcst$all.mean, 1, FUN=sd)
  fcst.lower <- fcst$mean + qnorm(pred.level2) * fcst.std
  fcst.upper <- fcst$mean - qnorm(pred.level2) * fcst.std
  
  fcst.df <- data.frame(dt=fcst.dts.seq, forecast=fcst$mean,
                        forecast.lower=fcst.lower,
                        forecast.upper=fcst.upper)
  
  fcst.df <- .BackTransformForecast(fcst.df, transform=transform,
                                    box.cox.lambda=box.cox.lambda)
  return(fcst.df)
  
}

#' Option for forecasting (aggregated or disaggregated, long term or short
#' term) using thief elm
#'
#' @param history: dataframe containing fields for dt, actual, actual.lower,
#'   actual.upper for forecasting
#' @param fcst.dts: named vector containing start and end dt of forecast
#' @param dt.units: units for the datetimes (dts) in the history, string
#' @param dt.format: string giving a formatting function for the dts
#' @param periods: length of seasonal period (in aggregated buckets of
#' periods.agg)
#' @param periods.agg: periods that training data (and forecast will be) are
#'   aggregated over
#' @param periods.trig: seasonal periods to fit trigonometric curves to as
#'   external regressors
#' @param pred.level: level for the prediction intervals
#' @param transform: the transformation applied and thus the inverse to apply
#'   to get forecast on the original scale
#' @param box.cox.lambda: the value of the box cox lambda in the event the
#'   transform is 'box_cox'
#' @param x.reg: an optional dataframe of external regressors
#' @param x.future an optional (but required if x.reg is not null) of future
#'    external regressors
#' @param holidays.df a dataframe containing holidays to model as features
#' @return: dataframe containing the dt, forecast, forecast lower prediction
#'   limit and upper prediction limit
Elm <- function(history, fcst.dts, dt.units='days', dt.format='.AsPOSIXlt',
                periods=c(7), periods.agg=NULL, periods.trig=c(364),
                pred.level=0.8, transform='box_cox', box.cox.lambda=1,
                x.reg=NULL, x.future=NULL, holidays.df=NULL) {
  
  fcst.interval <- max(c(1, periods.agg))
  dt_units <- get(dt.units)
  fcst.dts.seq <- seq((fcst.dts$begin.dt + dt_units(fcst.interval - 1)),
                      fcst.dts$end.dt, by=paste(fcst.interval, dt.units))
  len.fcst <- length(fcst.dts.seq)
  
  ts.training <- ts(history$actual, frequency=max(periods))
  if ((!is.null(x.reg) && !is.null(x.future)) ||
      (!is.null(holidays.df)) || (length(periods.trig) > 0)) {
    
    if (!is.null(holidays.df)) {
      
      holidays.matrix <- GetHolidayFeatures(history.start=
                                              (min(history$dt) -
                                                 dt_units(fcst.interval - 1)),
                                            fcst.dts=fcst.dts,
                                            dt.units=dt.units,
                                            dt.format=dt.format,
                                            periods.agg=periods.agg,
                                            holidays.df=holidays.df)
      x.reg <- cbind(x.reg, holidays.matrix$holidays.past)
      x.future <- cbind(x.future, holidays.matrix$holidays.future)
      
    }
    
    if (length(periods.trig) > 0) {
      trig.matrix <- GetTrigSeasonalityFeatures(nrow(history), len.fcst,
                                                periods.trig)
      
      x.reg <- cbind(x.reg, trig.matrix$trig.seasonality.past)
      x.future <- cbind(x.future, trig.matrix$trig.seasonality.future)
    }
    
    x.all <- rbind(x.reg, x.future)
    elm.fit <- suppressWarnings(nnfor::elm(ts.training, type='lasso',
                                           comb='median', lags=0, reps=50,
                                           xreg=x.all))
    fcst <- forecast::forecast(elm.fit, h=len.fcst,
                               level=pred.level, xreg=x.all)
  } else {
    elm.fit <- suppressWarnings(nnfor::elm(ts.training, type='lasso',
                                           comb='median', reps=50))
    fcst <- forecast::forecast(elm.fit, h=len.fcst,
                               level=pred.level)
  }
  
  pred.level2 <- (1 - pred.level) / 2
  
  fcst.std <- apply(fcst$all.mean, 1, FUN=sd)
  fcst.lower <- fcst$mean + qnorm(pred.level2) * fcst.std
  fcst.upper <- fcst$mean - qnorm(pred.level2) * fcst.std
  
  fcst.df <- data.frame(dt=fcst.dts.seq, forecast=fcst$mean,
                        forecast.lower=fcst.lower,
                        forecast.upper=fcst.upper)
  
  fcst.df <- .BackTransformForecast(fcst.df, transform=transform,
                                    box.cox.lambda=box.cox.lambda)
  return(fcst.df)
  
}

