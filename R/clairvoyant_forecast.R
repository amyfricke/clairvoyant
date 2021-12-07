#' Main forecasting function for clairvoyant package
#'
#' @param history: dataframe containing fields for dt_inde and actual to forecast
#' @param training.end.dt: last dt to be used in training
#' @param forecast.end.dt: dt to forecast until (may be moved back)
#' @param dt.units: units for the datetimes (dts) in the history, string
#' @param dt.format: string giving a formatting function for the dts
#' @param imputation.parameters: list of parameters for imputation obtained
#'  using ImputationParameters()
#' @param disaggregation.parameters: list of parameters for daily disaggregation
#'  obtained using DisaggregationParameters()
#' @param ensemble.parameters: list of parameters for the ensemble of models
#'  (such as the list of models, external regressors, etc) obtained using
#'   EnsembleParameters()
#' @param verbose: boolean for whether to output error or warning messages
#' @return: list of objects for the training data, steps in
#'  processing the training data, steps in obtaining the aggregated/unaggreated
#'  forecast, the final forecast and any available forecast error evaluation
#' @export
Clairvoyant <- function(history,
                        training.end.dt=NULL,
                        forecast.end.dt=NULL,
                        dt.units='days',
                        dt.format='.AsPOSIXlt',
                        imputation.parameters=NULL,
                        aggregation.parameters=NULL,
                        ensemble.parameters=NULL,
                        verbose=TRUE) {
  # Set parameters to their defaults if unset:
  if (is.null(imputation.parameters)) {
    imputation.parameters <- ImputationParameters(events.impute=NULL)
  }
  if (is.null(aggregation.parameters)) {
    aggregation.parameters <- AggregationParameters()
  }
  if (is.null(ensemble.parameters)) {
    ensemble.parameters <- EnsembleParameters()
  }
  # Make sure datetime, dt, is in specified format and ordered according to dt
  require(lubridate)
  dt_format <- get(dt.format)
  dt_units <- get(dt.units)
  history$dt <- dt_format(history$dt)
  history <- history[order(history$dt), ]
  first.dt <- min(history$dt)
  if (is.null(training.end.dt)) {
    training.end.dt <- max(history$dt)
  }

  training.end.dt <- dt_format(training.end.dt)


  # Subset the history to dates that yield an integer number of the largest
  # aggregated buckets for forecasting if desired
  periods.agg <- aggregation.parameters$periods.agg
  periods <- ensemble.parameters$periods
  if (is.null(forecast.end.dt)) {
    forecast.delta <- max(c(periods.agg, periods,
                            ensemble.parameters$periods.trig))
    forecast.end.dt <- training.end.dt + dt_units(forecast.delta + 1)
  }
  forecast.end.dt <- dt_format(forecast.end.dt)
  do.aggregation <- (length(periods.agg) > 0 &&
                     max(periods.agg) > 1)
  aggregate.to.longest <- aggregation.parameters$aggregate.to.longest
  train.dts <- .GetBeginEndDts(first.dt, training.end.dt,
                               dt.units=dt.units,
                               dt.format=dt.format,
                               fix.dt.option="begin.later",
                               need.integer.period=do.aggregation,
                               period=max(periods.agg))
  fcst.dts <- .GetBeginEndDts((train.dts$end.dt + dt_units(1)),
                               forecast.end.dt,
                               dt.units=dt.units,
                               dt.format=dt.format,
                               fix.dt.option="end.later",
                               need.integer.period=do.aggregation,
                               period=max(periods.agg))

  training <- list(unaggregated=list())
  forecast <- list()

  if (!is.null(ensemble.parameters$x.features) &&
      (min(ensemble.parameters$x.features$dt) <= train.dts$begin.dt) &&
      (max(ensemble.parameters$x.features$dt) >= fcst.dts$end.dt)) {
    training$x.features <- list()
    training$x.features$unaggregated <-
      subset(ensemble.parameters$x.features, dt >= train.dts$begin.dt
             & dt <= train.dts$end.dt)
    forecast$x.features <- list()
    forecast$x.features$unaggregated <-
      subset(ensemble.parameters$x.features, dt >= fcst.dts$begin.dt &
             dt <= fcst.dts$end.dt)

  } else if (is.null(ensemble.parameters$x.features)) {
    training$x.features$unaggregated <-
      forecast$x.features$unaggregated <- NULL
  } else {
    stop(paste('x.features must cover both the training and forecast period',
               'and must contain a column for dt'))
  }
  # Fill missing datetimes+ analyze missingness and impute if directed to
  filled.imputed.ob <-
    FillAnalyzeImpute(subset(history, dt >= train.dts$begin.dt
                             & dt <= train.dts$end.dt),
                      imputation.parameters,
                      dt.units,
                      dt.format,
                      training$x.features$unaggregated)

  training$unaggregated$raw <- filled.imputed.ob$filled
  if (imputation.parameters$impute.missing) {
    training$unaggregated$imputed <- filled.imputed.ob$imputed
    transform.tag <- 'imputed'
  } else transform.tag <- 'raw'
  if (imputation.parameters$analyze.missing) {
    training$missingness.analysis <- filled.imputed.ob$missingness.analysis
  }

  # Aggregate the data to weekly to eliminate one level of seasonality
  # then transform the data using the specified transformation (Box cox finds
  # the transformation that stabilizes the variance the best

  agg.transform.ob <-
    .AggregateAndTransform(training$unaggregated[[transform.tag]],
                           periods.agg=periods.agg,
                           agg.fun=aggregation.parameters$agg.fun,
                           periods=periods,
                           dt.format=dt.format,
                           transform=ensemble.parameters$transform)

  training$transformed <- list(df=agg.transform.ob$transformed,
                               transform=ensemble.parameters$transform,
                               box.cox.lambda=agg.transform.ob$box.cox.lambda)

  training$aggregated <- agg.transform.ob$aggregated

  if (!is.null(ensemble.parameters$x.features)) {
    x.cols <- colnames(ensemble.parameters$x.features)
    x.cols <- x.cols[x.cols != 'dt']
    if (aggregation.parameters$aggregate.to.longest &&
        max(periods.agg) > 1) {
      agg.x.ob <-
        AggregateToLongest(training$x.features$unaggregated,
                           cols.agg=x.cols,
                           periods.agg=periods.agg,
                           dt.format=dt.format,
                           agg.fun=aggregation.parameters$agg.fun)
      training$x.features$aggregated <- agg.x.ob
      x.reg <- subset(training$x.features$aggregated[[paste(max(periods.agg))]],
                      select=-c(dt))

      agg.x.ob <-
        AggregateToLongest(forecast$x.features$unaggregated,
                           cols.agg=x.cols,
                           periods.agg=periods.agg,
                           dt.format=dt.format,
                           agg.fun=aggregation.parameters$agg.fun)
      forecast$x.features$aggregated <- agg.x.ob

      x.future <-
        subset(forecast$x.features$aggregated[[paste(max(periods.agg))]],
               select=-c(dt))

    } else {
      x.reg <- subset(training$x.features$unaggregated,
                      select=-c(dt))
      x.future <- subset(forecast$x.features$unaggregated,
                         select=-c(dt))
    }
  } else {
    x.reg <- x.future <- NULL
  }

  # check that there is enough history to fit each seasonal period and remove
  # if not enough
  periods <- agg.transform.ob$periods
  if (length(periods) == 0) periods <- c(1)
  periods[periods < 1] <- 1
  for (p in periods[order(periods, decreasing=T)]) {
    if (nrow(training$transformed$df) <= 2 * p) {
      periods <- c(1, periods[periods != p])
      cat(paste("\n Seasonal period", p,
                "removed due to lack of sufficient history.\n"))
    }
  }

  # Fit an ensemble of models to the (imputed, aggregated) transformed data
  forecast$periods <- periods
  forecast$ensemble <- list()
  forecast$ensemble <- list()
  forecast$ensemble.summary <- list()
  agg.bucket <- max(c(1, periods.agg))
  models <- ensemble.parameters$models
  for (i in 1:length(models)) {
    ForecastMethod <- get(models[i])
    forecast.ob <- try(
          ForecastMethod(history=training$transformed$df,
                         fcst.dts=fcst.dts,
                         dt.units=dt.units, dt.format=dt.format,
                         periods=periods, periods.agg=periods.agg,
                         periods.trig=
                           ensemble.parameters$periods.trig / agg.bucket,
                         pred.level=ensemble.parameters$pred.level,
                         transform=ensemble.parameters$transform,
                         box.cox.lambda=training$transformed$box.cox.lambda,
                         x.reg=x.reg,
                         x.future=x.future,
                         holidays.df=ensemble.parameters$holidays.df),
        silent=!verbose)
    if (all(class(forecast.ob) != 'try-error')) {
      forecast$ensemble[[models[i]]] <- forecast.ob$forecast
      forecast$ensemble.trans[[models[i]]] <- forecast.ob$trans.forecast
      forecast$ensemble.summary[[models[i]]] <- forecast.ob$model.summary

    }
  }

  # Get the consensus forecast using the specified method
  forecast$consensus <-
    GetConsensusForecast(forecast$ensemble,
                         ensemble.parameters$consensus.method)

  range.methods <- ensemble.parameters$range.methods
  range <-
    GetForecastRange(forecast$ensemble.trans, lower.method=range.methods[1],
                     upper.method=range.methods[2],
                     ensemble.parameters$pred.level,
                     ensemble.parameters$transform,
                     box.cox.lambda=training$transformed$box.cox.lambda)
  range$range.se.lower <- range$range.se.upper <- range$c.alpha <- NULL
  forecast$consensus <- merge(forecast$consensus, range, by='dt')

  # Disaggregate the forecast back into unaggregated numbers using multinomial
  #  regression if the training data fed to the ensemble was aggregated

  if (aggregation.parameters$aggregate.to.longest &&
      (length(periods.agg) > 0)) {

    periods.disagg <- periods.agg[periods.agg != max(periods.agg)]
    periods.disagg <- periods.disagg[order(periods.disagg, decreasing=TRUE)]
    periods.disaggr <- c(periods.disagg, 1)
    periods.disagg <- periods.agg[order(periods.agg, decreasing=T)] /
      periods.disaggr
    periods.agg.r <- periods.agg[order(periods.agg, decreasing=TRUE)]
    periods.aggr <- c(max(periods), periods.disagg)

    disaggregate.format <- aggregation.parameters$disaggregate.format
    disaggregate.format <- disaggregate.format[order(periods.aggr,
                                                     decreasing=TRUE)]
    forecast$disaggregated <- list()
    forecast.aggregated <- forecast$consensus
    for (i in 1:length(periods.disaggr)) {

      if (periods.disaggr[i] == 1) {
        unaggregated.history <- training$unaggregated[[transform.tag]]
      } else {
        unaggregated.history <-
          training$aggregated[[paste0(periods.disaggr[i])]]
      }
      aggregated.history <-
        training$aggregated[[paste0(periods.agg.r[i])]]

      forecast$disaggregated[[paste0(periods.disaggr[i])]] <-
        DisaggregateForecast(unaggregated.history,
                             aggregated.history,
                             forecast.aggregated,
                             aggregation.parameters,
                             event.parameters,
                             periods.disagg[i],
                             periods.aggr[i],
                             dt.units=dt.units,
                             dt.format=dt.format,
                             disaggregate.format[i],
                             pred.level=ensemble.parameters$pred.level)
      forecast.aggregated <-
        forecast$disaggregated[[paste0(periods.disaggr[i])]]
    }

    forecast$final <- forecast$disaggregated[['1']]

    forecast$aggregated <- list()
    forecast$aggregated$final <-
      AggregateToLongest(forecast$final,
                         cols.agg=c('forecast',
                                    'forecast.lower',
                                    'forecast.upper',
                                    'range.forecast.lower',
                                    'range.forecast.upper',
                                    'range.uncertainty.lower',
                                    'range.uncertainty.upper'),
                         periods.agg=periods.agg,
                         dt.format=dt.format,
                         agg.fun=aggregation.parameters$agg.fun)


  } else {
    forecast$aggregated <- list(final=NULL)
    forecast$final <-forecast$consensus

  }

  # If there are data not used in training (after the train.end.dt) use this
  # to do out of sample evaluation of the forecasted values
  evaluation <- EvaluateClairvoyant(history, forecast$final,
                                    aggregation.parameters=
                                      aggregation.parameters,
                                    aggregated.histories=
                                      training$aggregated,
                                    aggregated.forecasts=
                                      forecast$aggregated$final,
                                    dt.format=dt.format,
                                    longest.period=max(periods.agg) * period,
                                    get.unagg.pred.int=
                                      aggregation.parameters$get.pred.int)

  # Do separate evaluation for each model in the ensemble
  evaluation$ensemble <- list()
  models <- names(forecast$ensemble)
  for (i in 1:length(models)) {
    if (length(periods.agg) > 0 && max(periods.agg) > 1) {
      aggregated.forecast <- list(
        forecast$ensemble[[models[i]]])
      final.forecast <- forecast$final
      aggregation.parameters$periods.agg <- c(max(periods.agg))
    } else {
      aggregated.forecast <- NULL
      final.forecast <- forecast$ensemble[[models[i]]]
    }
    evaluation.model.i <-
      EvaluateClairvoyant(history, final.forecast,
                         aggregation.parameters=aggregation.parameters,
                         aggregated.histories=training$aggregated,
                         aggregated.forecasts=aggregated.forecast,
                         dt.format=dt.format,
                         longest.period=max(periods.agg) * period,
                         get.unagg.pred.int=
                           aggregation.parameters$get.pred.int)
    if (length(periods.agg) > 0 && max(periods.agg) > 1) {
      evaluation$ensemble[[models[i]]] <-
        evaluation.model.i$aggregated[[paste(max(periods.agg))]]
    } else {
      evaluation$ensemble[[models[i]]] <- evaluation.model.i$unaggregated
    }
  }

  return(list(training=training,
              forecast=forecast,
              evaluation=evaluation))

}


