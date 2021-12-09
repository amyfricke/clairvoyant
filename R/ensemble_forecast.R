#' Function to create needed hyperparameters for ensemble forecasting
#'
#' @param transform: the transformation applied to the data to find ensemble
#'   forecast for
#' @param periods: length of seasonal periods that are not aggregated over
#'   or fit as trigonometric curves (in periods.trig)
#' @param periods.trig: seasonal periods to fit trigonometric curves to as
#'   external regressors
#' @param models: vector of the model names included in the ensemble
#' @param x.features: matrix of x.features to include as regressors in the
#'   ensemble
#' @param holidays.df a dataframe containing holidays to model as features
#' @param consensus.method: the function applied to get the consensus forecast
#' @param pred.level: confidence level for prediction/confidence interval
#' @return: A list containing the hyperparameters needed to smooth events and
#'   estimate their effects
#' @export
EnsembleParameters <- function(transform='none',
                               periods=c(7),
                               periods.trig=c(364),
                               models=c('Arima011', 'Arima111', 'Arima013',
                                        'Arima113',  'Arima112', 'Arima012',
                                        'AutoArima', 'Bsts', 'ProphetLinear'),
                               x.features=NULL,
                               holidays.df=NULL,
                               consensus.method='median',
                               range.methods=c('LowerQuartile', 'UpperQuartile'),
                               pred.level=0.8) {
  ensemble.parameters <- list()
  
  stopifnot(tolower(transform) %in% c('log', 'box_cox', 'none'))
  ensemble.parameters$transform <- tolower(transform)
  
  
  stopifnot(is.null(periods.trig) ||
              all(round(periods.trig) - periods.trig == 0))
  ensemble.parameters$periods.trig <- periods.trig
  
  stopifnot(all(round(periods) - periods == 0))
  # Remove any doubly specified periods (keep them in
  #  periods trig because those are easier to fit models for)
  periods <- periods[!(periods %in% periods.trig)]
  ensemble.parameters$periods <- periods
  
  stopifnot(class(models) == 'character')
  ensemble.parameters$models <- models
  
  #stopifnot(round(n.models.keep) - n.models.keep == 0)
  #ensemble.parameters$n.models.keep <- n.models.keep
  
  stopifnot(is.null(x.features) || any(class(x.features) %in%
                                         c('data.frame', 'matrix')))
  ensemble.parameters$x.features <- x.features
  
  stopifnot(is.null(holidays.df) || any(class(holidays.df) %in%
                                          c('data.frame', 'matrix')))
  ensemble.parameters$holidays.df <- holidays.df
  
  stopifnot(consensus.method %in% c('median', 'mean', 'LowerQuartile',
                                    'UpperQuartile'))
  ensemble.parameters$consensus.method <- consensus.method
  
  stopifnot(length(range.methods) == 2 &&
              all(range.methods %in% c('min', 'max',
                                       'LowerQuartile', 'UpperQuartile',
                                       'median', 'mean')))
  ensemble.parameters$range.methods <- range.methods
  
  stopifnot(pred.level >= 0 && pred.level <= 1)
  ensemble.parameters$pred.level <- pred.level
  
  
  return(ensemble.parameters)
}


#'  Function to generate arima forecast model with order (p, d, q)
#'  (and optionally, seasonal order 0,1,1)
#'
#' @param pdq.order vector containing the autoregressive order, degree of
#'  differencing and the moving average order (ordered in this manner)
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
.Arimapdq <- function(pdq.order=c(1, 1, 1), history, fcst.dts,
                      dt.units='days', dt.format='.AsPOSIXlt',
                      periods=c(7), periods.agg=NULL, periods.trig=c(364),
                      pred.level=0.8, transform='box_cox', box.cox.lambda=1,
                      x.reg=NULL, x.future=NULL, holidays.df=NULL) {
  
  fcst.interval <- max(c(1, periods.agg))
  dt_units <- get(dt.units)
  fcst.dts.seq <- seq((fcst.dts$begin.dt + dt_units(fcst.interval - 1)),
                      fcst.dts$end.dt, by=paste(fcst.interval, dt.units))
  len.fcst <- length(fcst.dts.seq)
  
  if (length(periods) > 0) {
    periods <- periods[periods > 1]
  }
  period.ts <- ifelse(length(periods) > 0, min(periods), 1)
  ts.training <- ts(history$actual, frequency=period.ts)
  ts.training.lower <- ts(history$actual.lower, frequency=period.ts)
  ts.training.upper <- ts(history$actual.upper, frequency=period.ts)
  
  
  if (any(periods > 1) && !is.null(periods)) {
    seasonal.order <- list(order=c(0, 1, 1), period=min(periods))
  } else {
    seasonal.order <- c(0, 0, 0)
  }
  
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
      if (!is.null(x.reg) && !is.null(x.future)) {
        x.reg <- cbind(x.reg, holidays.matrix$holidays.past)
        x.future <- cbind(x.future, holidays.matrix$holidays.future)
      } else {
        x.reg <- holidays.matrix$holidays.past
        x.future <- holidays.matrix$holidays.future
      }
      
    }
    
    if (length(periods.trig) > 0 && max(periods.trig) > 1) {
      trig.matrix <- GetTrigSeasonalityFeatures(nrow(history), len.fcst,
                                                periods.trig)
      if (!is.null(x.reg) && !is.null(x.future)) {
        x.reg <- cbind(x.reg, trig.matrix$trig.seasonality.past)
        x.future <- cbind(x.future, trig.matrix$trig.seasonality.future)
      } else {
        x.reg <- trig.matrix$trig.seasonality.past
        x.future <- trig.matrix$trig.seasonality.future
      }
    }
    
    x.reg <- as.matrix(x.reg)
    x.future <- as.matrix(x.future)
    arima.fit <- arima(ts.training, order=pdq.order, xreg=x.reg,
                       seasonal=seasonal.order)
    fcst <- predict(arima.fit, len.fcst, newxreg=x.future, se.fit=TRUE)
    
    arima.fit.lower <- arima(ts.training.lower, order=pdq.order,
                             seasonal=seasonal.order, xreg=x.reg)
    fcst.lower.pred <- predict(arima.fit.lower, len.fcst, newxreg=x.future)
    arima.fit.upper <- arima(ts.training.upper, order=pdq.order,
                             seasonal=seasonal.order, xreg=x.reg)
    fcst.upper.pred <- predict(arima.fit.upper, len.fcst, newxreg=x.future)
  } else {
    arima.fit <- arima(ts.training, order=pdq.order, xreg=x.reg,
                       seasonal=seasonal.order)
    fcst <- predict(arima.fit, len.fcst, newxreg=x.future, se.fit=TRUE)
    
    arima.fit.lower <- arima(ts.training.lower, order=pdq.order,
                             seasonal=seasonal.order, xreg=x.reg)
    fcst.lower.pred <- predict(arima.fit.lower, len.fcst, newxreg=x.future)
    arima.fit.upper <- arima(ts.training.upper, order=pdq.order,
                             seasonal=seasonal.order, xreg=x.reg)
    fcst.upper.pred <- predict(arima.fit.upper, len.fcst, newxreg=x.future)
  }
  fcst.lower <-
    fcst.lower.pred$pred -
    qnorm(1 - ((1 - pred.level) / 2)) * fcst.lower.pred$se
  
  fcst.upper <-
    fcst.upper.pred$pred +
    qnorm(1 - ((1 - pred.level) / 2)) * fcst.upper.pred$se
  
  fcst.df <- data.frame(dt=fcst.dts.seq, forecast=fcst$pred,
                        forecast.lower=fcst.lower,
                        forecast.upper=fcst.upper,
                        se=fcst$se)
  model.summary <- summary(arima.fit)
  fitted.df <- data.frame(dt=history$dt,
                          predicted=history$actual - arima.fit$residuals,
                          residuals=arima.fit$residuals)
  
  fcst.bcktrans <- .BackTransformForecast(fcst.df, transform=transform,
                                          box.cox.lambda=box.cox.lambda)
  fcst.bcktrans$forecast$se <- NULL
  
  return(list(forecast=fcst.bcktrans$forecast, trans.forecast=fcst.df,
              model.summary=model.summary, fitted=fitted.df))
  
}


#'  Option for forecasting (aggregated or disaggregated, long or short term)
#'  Finds the arima forecast with order (0, 1, 1) (and optionally, seasonal
#'  order 0,1,1)
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
#' @export
Arima011 <- function(history, fcst.dts,
                     dt.units='days', dt.format='.AsPOSIXlt',
                     periods=c(7), periods.agg=NULL, periods.trig=c(364),
                     pred.level=0.8, transform='box_cox', box.cox.lambda=1,
                     x.reg=NULL, x.future=NULL, holidays.df=NULL) {
  rslt <- .Arimapdq(c(0, 1, 1), history, fcst.dts, dt.units, dt.format,
                    periods, periods.agg, periods.trig, pred.level,
                    transform, box.cox.lambda, x.reg, x.future, holidays.df)
  
  
  return(list(forecast=rslt$forecast, trans.forecast=rslt$trans.forecast,
              model.summary=rslt$model.summary, fitted=rslt$fitted))
  
}


#' Option for forecasting (aggregated or disaggregated, long or short term)
#'  Finds the arima forecast with order (0, 1, 2) (and optionally, seasonal
#'  order 0,1,1)
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
#' @export
Arima012 <- function(history, fcst.dts,
                     dt.units='days', dt.format='.AsPOSIXlt',
                     periods=c(7), periods.agg=NULL, periods.trig=c(364),
                     pred.level=0.8, transform='box_cox', box.cox.lambda=1,
                     x.reg=NULL, x.future=NULL, holidays.df=NULL) {
  
  rslt <- .Arimapdq(c(0, 1, 2), history, fcst.dts, dt.units, dt.format,
                    periods, periods.agg, periods.trig, pred.level,
                    transform, box.cox.lambda, x.reg, x.future, holidays.df)
  
  
  return(list(forecast=rslt$forecast, trans.forecast=rslt$trans.forecast,
              model.summary=rslt$model.summary, fitted=rslt$fitted))
  
}


#' Option for forecasting (aggregated or disaggregated, long or short term)
#'  Finds the arima forecast with order (1, 1, 1) (and optionally, seasonal
#'  order 0,1,1)
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
#' @export
Arima111 <- function(history, fcst.dts,
                     dt.units='days', dt.format='.AsPOSIXlt',
                     periods=c(7), periods.agg=NULL, periods.trig=c(364),
                     pred.level=0.8, transform='box_cox', box.cox.lambda=1,
                     x.reg=NULL, x.future=NULL, holidays.df=NULL) {
  
  rslt <- .Arimapdq(c(1, 1, 1), history, fcst.dts, dt.units, dt.format,
                    periods, periods.agg, periods.trig, pred.level,
                    transform, box.cox.lambda, x.reg, x.future, holidays.df)
  
  
  return(list(forecast=rslt$forecast, trans.forecast=rslt$trans.forecast,
              model.summary=rslt$model.summary, fitted=rslt$fitted))
  
}

#' Option for forecasting (aggregated or disaggregated, long or short term)
#'  Finds the arima forecast with order (1, 1, 2) (and optionally, seasonal
#'  order 0,1,1)
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
#' @export
Arima112 <- function(history, fcst.dts,
                     dt.units='days', dt.format='.AsPOSIXlt',
                     periods=c(7), periods.agg=NULL, periods.trig=c(364),
                     pred.level=0.8, transform='box_cox', box.cox.lambda=1,
                     x.reg=NULL, x.future=NULL, holidays.df=NULL) {
  
  rslt <- .Arimapdq(c(1, 1, 2), history, fcst.dts, dt.units, dt.format,
                    periods, periods.agg, periods.trig, pred.level,
                    transform, box.cox.lambda, x.reg, x.future, holidays.df)
  
  
  return(list(forecast=rslt$forecast, trans.forecast=rslt$trans.forecast,
              model.summary=rslt$model.summary, fitted=rslt$fitted))
  
}

#' Option for forecasting (aggregated or disaggregated, long or short term)
#'  Finds the arima forecast with order (0, 1, 3) (and optionally, seasonal
#'  order 0,1,1)
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
#' @export
Arima013 <- function(history, fcst.dts,
                     dt.units='days', dt.format='.AsPOSIXlt',
                     periods=c(7), periods.agg=NULL, periods.trig=c(364),
                     pred.level=0.8, transform='box_cox', box.cox.lambda=1,
                     x.reg=NULL, x.future=NULL, holidays.df=NULL) {
  
  rslt <- .Arimapdq(c(0, 1, 3), history, fcst.dts, dt.units, dt.format,
                    periods, periods.agg, periods.trig, pred.level,
                    transform, box.cox.lambda, x.reg, x.future, holidays.df)
  
  
  return(list(forecast=rslt$forecast, trans.forecast=rslt$trans.forecast,
              model.summary=rslt$model.summary, fitted=rslt$fitted))
  
}

#' Option for forecasting (aggregated or disaggregated, long or short term)
#'  Finds the arima forecast with order (1, 1, 3) (and optionally, seasonal
#'  order 0,1,1)
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
#' @export
Arima113 <- function(history, fcst.dts,
                     dt.units='days', dt.format='.AsPOSIXlt',
                     periods=c(7), periods.agg=NULL, periods.trig=c(364),
                     pred.level=0.8, transform='box_cox', box.cox.lambda=1,
                     x.reg=NULL, x.future=NULL, holidays.df=NULL) {
  
  rslt <- .Arimapdq(c(1, 1, 3), history, fcst.dts, dt.units, dt.format,
                    periods, periods.agg, periods.trig, pred.level,
                    transform, box.cox.lambda, x.reg, x.future, holidays.df)
  
  
  return(list(forecast=rslt$forecast, trans.forecast=rslt$trans.forecast,
              model.summary=rslt$model.summary, fitted=rslt$fitted))
  
}

#' Option for forecasting (aggregated or disaggregated, long or short term)
#'  Finds the arima forecast with order (2, 1, 1) (and optionally, seasonal
#'  order 0,1,1)
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
#' @export
Arima211 <- function(history, fcst.dts,
                     dt.units='days', dt.format='.AsPOSIXlt',
                     periods=c(7), periods.agg=NULL, periods.trig=c(364),
                     pred.level=0.8, transform='box_cox', box.cox.lambda=1,
                     x.reg=NULL, x.future=NULL, holidays.df=NULL) {
  
  rslt <- .Arimapdq(c(2, 1, 1), history, fcst.dts, dt.units, dt.format,
                    periods, periods.agg, periods.trig, pred.level,
                    transform, box.cox.lambda, x.reg, x.future, holidays.df)
  
  
  return(list(forecast=rslt$forecast, trans.forecast=rslt$trans.forecast,
              model.summary=rslt$model.summary, fitted=rslt$fitted))
  
}

#' Option for forecasting (aggregated or disaggregated, long or short term)
#'  Finds the optimal arima forecast
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
#' @export
AutoArima <- function(history, fcst.dts,
                      dt.units='days', dt.format='.AsPOSIXlt',
                      periods=c(7), periods.agg=NULL, periods.trig=c(364),
                      pred.level=0.8, transform='box_cox', box.cox.lambda=1,
                      x.reg=NULL, x.future=NULL, holidays.df=NULL) {
  
  fcst.interval <- max(c(1, periods.agg))
  dt_units <- get(dt.units)
  fcst.dts.seq <- seq((fcst.dts$begin.dt + dt_units(fcst.interval - 1)),
                      fcst.dts$end.dt, by=paste(fcst.interval, dt.units))
  len.fcst <- length(fcst.dts.seq)
  
  if (length(periods) > 0) {
    periods <- periods[periods > 1]
  }
  period.ts <- ifelse(length(periods) > 0, min(periods), 1)
  ts.training <- ts(history$actual, frequency=period.ts)
  
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
      if (!is.null(x.reg) && !is.null(x.future)) {
        x.reg <- cbind(x.reg, holidays.matrix$holidays.past)
        x.future <- cbind(x.future, holidays.matrix$holidays.future)
      } else {
        x.reg <- holidays.matrix$holidays.past
        x.future <- holidays.matrix$holidays.future
      }
      
    }
    
    if (length(periods.trig) > 0 && max(periods.trig) > 1) {
      trig.matrix <- GetTrigSeasonalityFeatures(nrow(history), len.fcst,
                                                periods.trig)
      
      if (!is.null(x.reg) && !is.null(x.future)) {
        x.reg <- cbind(x.reg, trig.matrix$trig.seasonality.past)
        x.future <- cbind(x.future, trig.matrix$trig.seasonality.future)
      } else {
        x.reg <- trig.matrix$trig.seasonality.past
        x.future <-trig.matrix$trig.seasonality.future
      }
      
    }
    
    autoarima.fit <- suppressWarnings(
      forecast::auto.arima(ts.training, biasadj=T, max.D=1, max.P=1, max.Q=1,
                           max.p=3, max.q=3, parallel=TRUE,
                           xreg=as.matrix(x.reg)))
    
    fcst <- forecast::forecast(autoarima.fit, h=len.fcst, level=pred.level,
                               xreg=as.matrix(x.future))
  } else {
    autoarima.fit <- suppressWarnings(
      forecast::auto.arima(ts.training, biasadj=T, max.D=1, max.P=1, max.Q=1,
                           max.p=3, max.q=3, parallel=TRUE))
    
    fcst <- forecast::forecast(autoarima.fit, h=len.fcst, level=pred.level)
    
  }
  
  fcst.se <- (as.vector(fcst$upper) - as.vector(fcst$lower)) /
    qnorm(0.5 + pred.level / 2)
  fcst.df <- data.frame(dt=fcst.dts.seq, forecast=fcst$mean,
                        forecast.lower=as.vector(fcst$lower),
                        forecast.upper=as.vector(fcst$upper),
                        se=fcst.se)
  
  model.summary <- summary(fcst)$model
  fitted.df <- data.frame(dt=history$dt, predicted=fcst$fitted,
                          residuals=fcst$residuals)
  
  fcst.bcktrans <- .BackTransformForecast(fcst.df, transform=transform,
                                          box.cox.lambda=box.cox.lambda)
  
  fcst.bcktrans$forecast$se <- NULL
  
  return(list(forecast=fcst.bcktrans$forecast, trans.forecast=fcst.df,
              model.summary=model.summary, fitted=fitted.df))
  
}


#' Option for forecasting (aggregated or disaggregated, long term or short
#' term) using bsts forecast
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
#' @export
Bsts <- function(history, fcst.dts, dt.units='days', dt.format='.AsPOSIXlt',
                 periods=c(7), periods.agg=NULL, periods.trig=c(364),
                 pred.level=0.8, transform='box_cox', box.cox.lambda=1,
                 x.reg=NULL, x.future=NULL, holidays.df=NULL) {
  
  set.seed(54321)
  fcst.interval <- max(c(1, periods.agg))
  dt_units <- get(dt.units)
  fcst.dts.seq <- seq((fcst.dts$begin.dt + dt_units(fcst.interval - 1)),
                      fcst.dts$end.dt, by=paste(fcst.interval, dt.units))
  len.fcst <- length(fcst.dts.seq)
  
  ts.training <- ts(history$actual, frequency=max(periods))
  pred.level2 <- (1 - pred.level) / 2
  
  if (length(periods.trig) > 0 && max(periods.trig) > 1) {
    freqs <- max(periods.trig) / periods.trig
  }
  
  
  if ((!is.null(x.reg) && !is.null(x.future)) ||
      (!is.null(holidays.df))) {
    
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
    
    for (xcol in colnames(x.reg)) {
      assign(xcol, x.reg[[xcol]])
    }
    bsts.form <- paste('ts.training ~', paste(colnames(x.reg), collapse=' + '))
    bsts.form <- as.formula(bsts.form)
    ss <- bsts::AddDynamicRegression(list(), bsts.form)
    ss <- bsts::AddSemilocalLinearTrend(ss, y=ts.training)
    bsts.regressors <- c('intercept', colnames(x.reg), 'trend.ar')
    if (length(periods.trig) > 0 && max(periods.trig) > 1) {
      ss <- bsts::AddTrig(ss, y=ts.training, period=max(periods.trig),
                          frequencies=freqs)
      bsts.regressors <- c(bsts.regressors, paste0('trig.', max(periods.trig)))
    }
    
    
    if (length(periods) > 0 && max(periods) > 1) {
      ss <- bsts::AddSeasonal(ss, y=ts.training,
                              nseasons=min(periods[periods > 1]))
      bsts.regressors <- c(bsts.regressors,
                           paste0('seasonal.', min(periods[periods > 1])))
    }
    
    bsts.fit <- try(bsts::bsts(bsts.form, state.specification=ss, niter=1000),
                    silent=F)
    fcst <- predict(bsts.fit, horizon=len.fcst, burn=200,
                    quantiles=c(pred.level2, 1-pred.level2), newdata=x.future)
    bsts.coefs <- c(colMeans(bsts.fit$coefficients),
                    mean(bsts.fit$trend.slope.ar.coefficient))
    bsts.sd <- c(apply(bsts.fit$coefficients, 2, sd),
                 sd(bsts.fit$trend.slope.ar.coefficient))
    
  } else {
    ss <- bsts::AddSemilocalLinearTrend(list(), y=ts.training)
    bsts.regressors <- c('trend.ar')
    if (length(periods.trig) > 0 && max(periods.trig) > 1) {
      ss <- bsts::AddTrig(ss, y=ts.training, period=max(periods.trig),
                          frequencies=freqs)
      bsts.regressors <- c(bsts.regressors, paste0('trig.', max(periods.trig)))
    }
    
    if (length(periods) > 0 && max(periods) > 1) {
      ss <- bsts::AddSeasonal(ss, y=ts.training, nseasons=periods)
      bsts.regressors <- c(bsts.regressors,
                           paste0('seasonal.', min(periods[periods > 1])))
    }
    bsts.fit <- try(bsts::bsts(ts.training, state.specification=ss, niter=1000),
                    silent=T)
    fcst <- predict(bsts.fit, horizon=len.fcst, burn=200,
                    quantiles=c(pred.level2, 1-pred.level2))
    bsts.coefs <- c(mean(bsts.fit$trend.slope.ar.coefficient))
    bsts.sd <- c(sd(bsts.fit$trend.slope.ar.coefficient))
    
  }
  
  if (length(periods.trig) > 0 && max(periods.trig) > 1) {
    bsts.coefs <- c(bsts.coefs, NA)
    bsts.sd <- c(bsts.sd, mean(bsts.fit[[paste0('trig.coefficient.sd.',
                                                max(periods.trig))]]))
  }
  if (length(periods) > 0 && max(periods) > 1) {
    min.periods.g1 <- min(periods[periods > 1])
    bsts.coefs <- c(bsts.coefs, NA)
    bsts.sd <- c(bsts.sd, sqrt(mean(bsts.fit[[paste0('sigma.seasonal.',
                                                     min.periods.g1)]])))
  }
  
  
  fcst.lower <- fcst$interval[1,]
  fcst.upper <- fcst$interval[2,]
  fcst.df <- data.frame(dt=fcst.dts.seq, forecast=fcst$median,
                        forecast.lower=fcst.lower,
                        forecast.upper=fcst.upper,
                        se=apply(fcst$distribution, 2, sd) / sqrt(800))
  model.summary <- data.frame(regressor=bsts.regressors,
                              coef=bsts.coefs,
                              sd=bsts.sd)
  model.summary$log.likelihood <- bsts.fit$log.likelihood[1000]
  
  bsts.residuals <-
    colMeans(bsts::bsts.prediction.errors(bsts.fit, burn=200)$in.sample)
  fitted.df <- data.frame(dt=history$dt,
                          predicted=history$actual - bsts.residuals,
                          residuals=bsts.residuals)
  
  fcst.bcktrans <- .BackTransformForecast(fcst.df, transform=transform,
                                          box.cox.lambda=box.cox.lambda)
  
  fcst.bcktrans$forecast$se <- NULL
  
  return(list(forecast=fcst.bcktrans$forecast, trans.forecast=fcst.df,
              model.summary=model.summary, fitted=fitted.df))
  
}

#' Option for forecasting (aggregated or disaggregated, long or short term).
#'  Finds a forecast using feed forward neural net autoregressive model
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
#' @export
Nnetar <- function(history, fcst.dts, dt.units='days', dt.format='.AsPOSIXlt',
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
      
      if (!is.null(holidays.df)) {
        
        holidays.matrix <- GetHolidayFeatures(history.start=
                                                (min(history$dt) -
                                                   dt_units(fcst.interval - 1)),
                                              fcst.dts=fcst.dts,
                                              dt.units=dt.units,
                                              dt.format=dt.format,
                                              periods.agg=periods.agg,
                                              holidays.df=holidays.df)
        if (!is.null(x.reg) && !is.null(x.future)) {
          x.reg <- cbind(x.reg, holidays.matrix$holidays.past)
          x.future <- cbind(x.future, holidays.matrix$holidays.future)
        } else {
          x.reg <- holidays.matrix$holidays.past
          x.future <- holidays.matrix$holidays.future
        }
        
      }
      
      if (length(periods.trig) > 0 && max(periods.trig) > 1) {
        trig.matrix <- GetTrigSeasonalityFeatures(nrow(history), len.fcst,
                                                  periods.trig)
        
        if (!is.null(x.reg) && !is.null(x.future)) {
          x.reg <- cbind(x.reg, trig.matrix$trig.seasonality.past)
          x.future <- cbind(x.future, trig.matrix$trig.seasonality.future)
        } else {
          x.reg <- trig.matrix$trig.seasonality.past
          x.future <-trig.matrix$trig.seasonality.future
        }
      }
    }
    nnetar.fit <- suppressWarnings(forecast::nnetar(ts.training, xreg=x.reg))
    
    fcst <- forecast::forecast(nnetar.fit, h=len.fcst, PI=TRUE,
                               level=pred.level, xreg=x.future)
  } else {
    nnetar.fit <- suppressWarnings(forecast::nnetar(ts.training))
    
    fcst <- forecast::forecast(nnetar.fit, h=len.fcst, PI=TRUE,
                               level=pred.level)
  }
  fcst.df <- data.frame(dt=fcst.dts.seq, forecast=fcst$mean,
                        forecast.lower=as.vector(fcst$lower),
                        forecast.upper=as.vector(fcst$upper))
  fcst.df$se <- (fcst.df$forecast.upper - fcst.df$forecast.lower) /
    (qnorm(0.5 + pred.level / 2))
  
  model.summary <- NULL
  fitted.df <- data.frame(dt=history$dt,
                          predicted=nnetar.fit$fitted,
                          residuals=history$actual - nnetar.fit$fitted)
  
  fcst.bcktrans <- .BackTransformForecast(fcst.df, transform=transform,
                                          box.cox.lambda=box.cox.lambda)
  fcst.bcktrans$forecast$se <- NULL
  
  return(list(forecast=fcst.bcktrans$forecast, trans.forecast=fcst.df,
              model.summary=model.summary, fitted=fitted.df))
  
}


#' Option for forecasting (aggregated or disaggregated, long term or short
#' term).  Finds the prophet forecast with linear trend.
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
#' @export
ProphetLinear <- function(history, fcst.dts,
                          dt.units='days', dt.format='.AsPOSIXlt',
                          periods=c(7), periods.agg=NULL, periods.trig=c(364),
                          pred.level=0.8, transform='box_cox', box.cox.lambda=1,
                          x.reg=NULL, x.future=NULL, holidays.df=NULL) {
  
  fcst.interval <- max(c(1, periods.agg))
  dt_units <- get(dt.units)
  fcst.by <-
    ifelse(fcst.interval==1, substr(dt.units, 1, (nchar(dt.units) - 1)),
           paste(fcst.interval, dt.units))
  fcst.dts.seq <- seq((fcst.dts$begin.dt + dt_units(fcst.interval - 1)),
                      fcst.dts$end.dt, by=fcst.by)
  len.fcst <- length(fcst.dts.seq)
  
  all.periods <- fcst.interval * c(1, periods, periods.trig)
  
  colnames(history) <- c('ds', 'y', 'y.lower', 'y.upper')
  
  if (!is.null(holidays.df)) {
    holidays.df <- holidays.df[c('dt', 'holiday')]
    colnames(holidays.df) <- c('ds', 'holiday')
  }
  prophet.m <- prophet::prophet(history, fit=F,
                                interval.width=pred.level,
                                yearly.seasonality=(any(all.periods %in%
                                                          c(364, 365))),
                                holidays=holidays.df,
                                mcmc.samples=1000,
                                changepoint.prior.scale=0.01)
  
  p.fcst.df <- data.frame(ds=fcst.dts.seq)
  periods.st <- c(1, 7, 364, 365, 24)
  periods.g <- fcst.interval * c(periods, periods.trig)
  
  if (!(all(periods.g %in% periods.st))) {
    for (period.ns in periods.g[!(periods.g %in% periods.st)]) {
      prophet.m <-
        prophet::add_seasonality(prophet.m,
                                 paste0('season_', period.ns),
                                 period.ns,
                                 fourier.order=2)
    }
  }
  
  
  if (!is.null(x.reg) && !is.null(x.future)) {
    for (x.reg.name in names(x.reg)){
      history[paste0('x.reg.', x.reg.name)] <- x.reg[x.reg.name]
      prophet.m <-
        prophet::add_regressor(prophet.m, paste0('x.reg.', x.reg.name))
      p.fcst.df[paste0('x.reg.', x.reg.name)] <-
        x.future[x.reg.name]
    }
    
  }
  prophet.f <- prophet::fit.prophet(prophet.m, history)
  p.fcst.df <- predict(prophet.f, p.fcst.df)
  
  fcst.df <- data.frame(dt=fcst.dts.seq, forecast=p.fcst.df$yhat,
                        forecast.lower=p.fcst.df$yhat_lower,
                        forecast.upper=p.fcst.df$yhat_upper,
                        se=(p.fcst.df$yhat_upper-p.fcst.df$yhat_lower)/
                          (qnorm(0.5 + pred.level / 2)))
  
  if (!is.null(x.reg) && !is.null(x.future)) {
    model.summary <- regressor_coefficients(prophet.f)
  } else {
    model.summary <- NULL
  }
  
  model.summary <- rbind(model.summary,
                         data.frame(regressor=c('trend'),
                                    regressor_mode=c('additive'),
                                    center=c(0),
                                    coef_lower=p.fcst.df$trend_lower[2]-
                                      p.fcst.df$trend_lower[1],
                                    coef=p.fcst.df$trend[2]-
                                      p.fcst.df$trend[1],
                                    coef_upper=p.fcst.df$trend_upper[2]-
                                      p.fcst.df$trend_upper[1]))
  model.summary$se <- ((model.summary$coef_upper - model.summary$coef_lower) /
                         (qnorm(0.5 + pred.level / 2)))
  
  prophet.fitted <- predict(prophet.f, history)
  
  fitted.df <- data.frame(dt=history$ds,
                          predicted=prophet.fitted$yhat,
                          residuals=history$y - prophet.fitted$yhat)
  
  fcst.bcktrans <- .BackTransformForecast(fcst.df, transform=transform,
                                          box.cox.lambda=box.cox.lambda)
  
  fcst.bcktrans$forecast$se <- NULL
  
  return(list(forecast=fcst.bcktrans$forecast, trans.forecast=fcst.df,
              model.summary=model.summary, fitted=fitted.df))
  
}

#' Option for forecasting (aggregated or disaggregated, long term or short
#' term).  Finds the prophet logistic forecast.
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
#' @export
ProphetLogistic <- function(history, fcst.dts,
                            dt.units='days', dt.format='.AsPOSIXlt',
                            periods=c(7), periods.agg=NULL, periods.trig=c(364),
                            pred.level=0.8, transform='box_cox', box.cox.lambda=1,
                            x.reg=NULL, x.future=NULL, holidays.df=NULL) {
  
  fcst.interval <- max(c(1, periods.agg))
  dt_units <- get(dt.units)
  fcst.dts.seq <- seq((fcst.dts$begin.dt + dt_units(fcst.interval - 1)),
                      fcst.dts$end.dt, by=paste(fcst.interval, dt.units))
  len.fcst <- length(fcst.dts.seq)
  
  colnames(history) <- c('ds', 'y', 'y.lower', 'y.upper')
  # TODO: change this to a manipulable parameter
  history$cap <- 10 * history$y
  
  if (!is.null(holidays.df)) {
    holidays.df <- holidays.df[c('dt', 'holiday')]
    colnames(holidays.df) <- c('ds', 'holiday')
  }
  all.periods <- fcst.interval * c(1, periods, periods.trig)
  prophet.m <- prophet::prophet(history, fit=F,
                                interval.width=pred.level,
                                yearly.seasonality=(any(all.periods %in%
                                                          c(364, 365))),
                                holidays=holidays.df,
                                #mcmc.samples=1000,
                                changepoint.prior.scale=0.01)
  
  p.fcst.df <- data.frame(ds=fcst.dts.seq, cap=max(history$cap))
  
  periods.st <- c(1, 7, 364, 365, 24)
  periods.g <- fcst.interval * c(periods, periods.trig)
  
  if (!(all(periods.g %in% periods.st))) {
    for (period.ns in periods.g[!(periods.g %in% periods.st)]) {
      prophet.m <-
        prophet::add_seasonality(prophet.m,
                                 paste0('season_', period.ns),
                                 period.ns,
                                 fourier.order=2)
    }
  }
  
  if (!is.null(x.reg) && !is.null(x.future)) {
    for (x.reg.name in names(x.reg)){
      history[paste0('x.reg.', x.reg.name)] <- x.reg[x.reg.name]
      prophet.m <-
        prophet::add_regressor(prophet.m, paste0('x.reg.', x.reg.name))
      p.fcst.df[paste0('x.reg.', x.reg.name)] <-
        x.future[x.reg.name]
    }
    
  }
  prophet.f <- prophet::fit.prophet(prophet.m, history)
  p.fcst.df <- predict(prophet.f, p.fcst.df)
  
  fcst.df <- data.frame(dt=fcst.dts.seq, forecast=p.fcst.df$yhat,
                        forecast.lower=p.fcst.df$yhat_lower,
                        forecast.upper=p.fcst.df$yhat_upper,
                        se=(p.fcst.df$yhat_upper-p.fcst.df$yhat_lower)/
                          (qnorm(0.5 + pred.level / 2)))
  
  if (!is.null(x.reg) && !is.null(x.future)) {
    model.summary <- regressor_coefficients(prophet.f)
  } else {
    model.summary <- NULL
  }
  model.summary <- rbind(model.summary,
                         data.frame(regressor=c('trend'),
                                    regressor_mode=c('additive'),
                                    center=c(0),
                                    coef_lower=p.fcst.df$trend_lower[2]-
                                      p.fcst.df$trend_lower[1],
                                    coef=p.fcst.df$trend[2]-
                                      p.fcst.df$trend[1],
                                    coef_upper=p.fcst.df$trend_upper[2]-
                                      p.fcst.df$trend_upper[1]))
  model.summary$se <- ((model.summary$coef_upper - model.summary$coef_lower) /
                         (qnorm(0.5 + pred.level / 2)))
  
  prophet.fitted <- predict(prophet.f, history)
  
  fitted.df <- data.frame(dt=history$ds,
                          predicted=prophet.fitted$yhat,
                          residuals=history$y - prophet.fitted$yhat)
  
  fcst.bcktrans <- .BackTransformForecast(fcst.df, transform=transform,
                                          box.cox.lambda=box.cox.lambda)
  
  fcst.bcktrans$forecast$se <- NULL
  
  return(list(forecast=fcst.bcktrans$forecast, trans.forecast=fcst.df,
              model.summary=model.summary, fitted=fitted.df))
  
}



#' Function to get a consensus forecast from an ensemble of forecasts
#'
#' @param forecast.list: list of forecast dataframes (need to contain fields
#'  for dt, forecast, forecast.lower, foreast.upper)
#' @param consensus.method: function to apply to ensemble of forecasts
#'  (pointwise)
#' @return: a dataframe containing the consensus forecast
#' @export
GetConsensusForecast <- function(forecast.list, consensus.method='median') {
  
  if (length(forecast.list) == 1) return(forecast.list[[1]])
  fcst.length <- length(forecast.list[[1]]$dt)
  num.fcsts <- length(forecast.list)
  fcst.mat <- fcst.upper.mat <- fcst.lower.mat <-
    matrix(rep(NA, fcst.length * num.fcsts), ncol=num.fcsts)
  
  for (i in 1:num.fcsts) {
    fcst.mat[, i] <- forecast.list[[i]]$forecast
    fcst.upper.mat[, i] <- forecast.list[[i]]$forecast.upper
    fcst.lower.mat[, i] <- forecast.list[[i]]$forecast.lower
  }
  consensus_fun <- get(consensus.method)
  fcst.consensus <- apply(fcst.mat, 1, consensus_fun, na.rm=T)
  fcst.upper.consensus <- apply(fcst.upper.mat, 1, consensus_fun, na.rm=T)
  fcst.lower.consensus <- apply(fcst.lower.mat, 1, consensus_fun, na.rm=T)
  
  consensus.forecast <- data.frame(dt=forecast.list[[1]]$dt,
                                   forecast=fcst.consensus,
                                   forecast.lower=fcst.lower.consensus,
                                   forecast.upper=fcst.upper.consensus)
  return(consensus.forecast)
}

#' Helper function to match the standard error to the forecast it corresponds
#' to
#'
#' @param se.mat: matrix of standard errors for forecasts of each model in the
#'  ensemble
#' @param fcst: forecast for which to match the corresponding standard error
#' @param fcst.mat matrix of forecasts from each model in the ensemble
#' @noRd
.MatchSe <- function(fcst, se.mat, fcst.mat) {
  se.idx <- which.min(abs(fcst.mat - fcst))
  se.fcst <- se.mat[se.idx]
  if (is.null(se.fcst)) {
    se.fcst <- NA
  }
  return(se.fcst)
}

#' Function to get a forecast sensitivty range from an ensemble of forecasts
#'
#' @param forecast.list: list of forecast dataframes (need to contain fields
#'   for dt, forecast, forecast.lower, foreast.upper)
#' @param lower.method: function to apply to ensemble of forecasts
#'  (pointwise) to get the lower end of the forecast sensitivity range
#' @param upper.method: function to apply to ensemble of forecasts
#'  (pointwise) to get the upper end of the forecast sensitivity range
#' @return: a dataframe containing the range forecast and uncertainty
#' @export
GetForecastRange <- function(forecast.list,
                             lower.method='min', upper.method='max',
                             pred.level=0.8, transform='none',
                             box.cox.lambda=1) {
  
  cols.transf <- c("range.forecast.lower", "range.forecast.upper",
                   "range.uncertainty.lower", "range.uncertainty.upper")
  
  if (length(forecast.list) == 1) {
    forecast.range <- data.frame(dt=forecast.list[[1]]$dt)
    forecast.range$range.forecast.lower <- forecast.list[[1]]$forecast
    forecast.range$range.forecast.upper <- forecast.list[[1]]$forecast
    forecast.range$range.uncertainty.lower <- forecast.list[[1]]$forecast.lower
    forecast.range$range.uncertainty.upper <- forecast.list[[1]]$forecast.upper
    
    forecast.range <-
      .BackTransformForecast(forecast.range, cols.transform=cols.transf,
                             transform=transform,
                             box.cox.lambda=box.cox.lambda)
    return(forecast.range$forecast)
  }
  fcst.length <- length(forecast.list[[1]]$dt)
  num.fcsts <- length(forecast.list)
  fcst.mat <- fcst.se.mat <-
    matrix(rep(NA, fcst.length * num.fcsts), ncol=num.fcsts)
  
  for (i in 1:num.fcsts) {
    fcst.mat[, i] <- forecast.list[[i]]$forecast
    fcst.se.mat[, i] <- forecast.list[[i]]$se
  }
  
  lower_fun <- get(lower.method)
  upper_fun <- get(upper.method)
  range.fcst.lower <- apply(fcst.mat, 1, lower_fun, na.rm=T)
  range.fcst.upper <- apply(fcst.mat, 1, upper_fun, na.rm=T)
  
  fcst.se.lower <- apply(as.matrix(range.fcst.lower, ncol=1), 1,
                         .MatchSe, fcst.se.mat, fcst.mat)
  
  fcst.se.upper <- apply(as.matrix(range.fcst.upper, ncol=1), 1,
                         .MatchSe, fcst.se.mat, fcst.mat)
  forecast.range <- data.frame(dt=forecast.list[[1]]$dt,
                               range.forecast.lower=range.fcst.lower,
                               range.forecast.upper=range.fcst.upper,
                               range.se.lower=fcst.se.lower,
                               range.se.upper=fcst.se.upper)
  
  c.alphas <- apply(forecast.range, 1, GetAdjustedAlpha, alpha=(1-pred.level))
  c.alphas <- as.vector(unlist(c.alphas))
  
  forecast.range$range.uncertainty.lower <-
    (forecast.range$range.forecast.lower -
       c.alphas * forecast.range$range.se.lower)
  forecast.range$range.uncertainty.upper <-
    (forecast.range$range.forecast.upper +
       c.alphas * forecast.range$range.se.upper)
  
  forecast.range <-
    .BackTransformForecast(forecast.range, cols.transform=cols.transf,
                           transform=transform, box.cox.lambda=box.cox.lambda)
  forecast.range <- forecast.range$forecast
  forecast.range$c.alpha <- c.alphas
  
  return(forecast.range)
}
