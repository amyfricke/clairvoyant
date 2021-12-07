
#' Function to create needed hyper parameters for any imputation of
#' missing or other datetimes with unknown values
#'
#' @param missing.fill: dataframe containing columns for dt and event
#' @param analyze.missing: number of days before the event to smooth/estimate effect
#' @param impute.missing: number of days after the event to smooth/estimate effect
#' @param events.impute: number of weeks before the event to use in smoothing/
#'    effect estimation
#' @param method: which forecast function to use to impute if
#'  missing datetimes or other events where the value of the ts is unknown
#' @param transform: which
#' @param periods: which
#' @param pred.level: confidence level for prediction/confidence interval
#' @return: A list containing the hyperparameters needed to smooth events and
#'   estimate their effects
#' @export
ImputationParameters <- function(missing.fill=NA,
                                 analyze.missing=TRUE,
                                 impute.missing=FALSE,
                                 events.impute=NULL,
                                 method='AutoArima',
                                 transform='none',
                                 periods=c(7),
                                 pred.level=0.8,
                                 length.history=10 * max(periods),
                                 use.later.history=TRUE) {
  imputation.parameters <- list()

  stopifnot(is.numeric(missing.fill) || is.na(missing.fill))
  imputation.parameters$missing.fill <- missing.fill

  stopifnot(analyze.missing %in% c(TRUE, FALSE))
  imputation.parameters$analyze.missing <- analyze.missing

  stopifnot(impute.missing %in% c(TRUE, FALSE))
  imputation.parameters$impute.missing <- impute.missing

  stopifnot(c('dt', 'event') %in% colnames(events.impute) ||
            is.null(events.impute))
  imputation.parameters$events.impute <- events.impute

  stopifnot(method %in% c('AutoArima', 'Nnetar',
                          'Mlp', 'Elm', 'Arima011', 'Arima012',
                          'Arima021', 'Bsts'))
  imputation.parameters$method <- method

  stopifnot(tolower(transform) %in% c('none'))
  imputation.parameters$transform <- tolower(transform)

  stopifnot(all(ceiling(periods) %% periods == 0))
  imputation.parameters$periods <- periods

  stopifnot(pred.level >= 0 && pred.level <= 1)
  imputation.parameters$pred.level <- pred.level

  stopifnot(ceiling(length.history) %% length.history == 0)
  imputation.parameters$length.history <- length.history

  stopifnot(use.later.history %in% c(TRUE, FALSE))
  imputation.parameters$use.later.history <- use.later.history

  return(imputation.parameters)
}

#' Function to take data frame with potentially missing values for dts and
#' replace them with fill.default
#'
#' @param history A data frame with two mandatory fields of dt and actual
#' @param fill.default To replace the missing values if fill.default is not
#'   provided NA's will be used.
#' @param dt.units datetime units in the history
#' @param dt.format format function for the datetimes
#' @return The timeseries with missing dates filled by fill.default.
#' @export
FillMissingDts <- function(history, fill.default = NA,
                           dt.units='days',
                           dt.format='.AsPOSIXlt') {
  stopifnot(is.na(fill.default) | is.numeric(fill.default))
  colnames(history) <- c("dt", "actual")
  dt_format <- get(dt.format)
  dt_units <- get(dt.units)
  dt.unit <- substr(dt.units, 1, (nchar(dt.units) - 1))
  history$dt <- dt_format(history$dt)
  history$actual <- as.numeric(history$actual)
  # The starting dt of new.history would be the first dt in history
  start.dt <- min(history$dt, na.rm = TRUE)
  # The end dt of new.history would be the maximum dt in history
  end.dt <- max(history$dt, na.rm = TRUE)
  history <- history[history$dt >= start.dt & history$dt <= end.dt, ]
  all.dts <- seq(from = start.dt, to = end.dt, by = dt.unit)
  # Find missing dates
  missing.dts <- setdiff(as.character(all.dts), as.character(history$dt))
  missing.data <- data.frame(
    dt=dt_format(missing.dts),
    actual=rep(fill.default, length(missing.dts))
  )
  new.history <- rbind(history, missing.data)
  new.history <- new.history[order(new.history$dt), ]
  row.names(new.history) <- 1:length(new.history$actual)
  return(new.history)
}


#' Function to analyze missingness mechanism for patterns in the external
#' regressors, seasonal periods or holidays
#'
#' @param history: timeseries history to analyze missingness pattern on
#' @param x.reg: external regressors to evaluate missingness pattern on
#' @param alpha.sig: the significance level for inferring whether a feature
#'  predicts missingness
#' @return: A model summary indicating whether the missingness depends
#'   on any of the external regressors
#' @export
AnalyzeMissingness <- function(history, x.reg, alpha.sig) {

  if (mean(is.na(history$actual)) < 0.10) {
    return('Not enough missing data to analyze pattern.')
  }
  history$is.missing <- is.na(history$actual)
  comb.history <- merge(history, x.reg, by=dt)
  model.features <- colnames(comb.history)
  not.features <- c('dt', 'actual', 'is.missing', 'actual.upper',
                    'actual.lower')
  model.features <- model.features[!(model.features %in% not.features)]
  model.statement <- paste(model.features, collapse=' + ')
  model.statement <- paste('is.missing ~', model.statement)
  lr.model <- glm(formula(model.statement), family='binomial')
  lr.model.summary <- summary(lr.model)
  features.summary <- lr.model.summary$coeficients
  features.significant <-
    colnames(features.summary)[features.summary[,4] < alpha.sig]
  features.significant <-
    features.significant[features.significant != "(Intercept)"]
  return(list(missingness.summary=lr.model.summary,
              features.significant=features.significant))

}

#' Function to put neighboring missing dts into clumps of training
#' end indices for imputation (reduces the number of imputations to run)
#'
#' @param history: timeseries history to find imputation training indices
#' @param imp.params: an ImputationParameters object giving the parameters
#'   used in imputation
#' @return: A data frame containing the indices that imputation will end on as
#'  well as the missing or unknown dts each end index corresponds to
#' @noRd
GetImputationTrainingIndices <- function(history, imp.params,
                                         dt.format='.AsPOSIXlt') {
  if (is.na(imp.params$missing.fill)) {
    missing.index <- is.na(history$actual)
  } else {
    missing.index <- (history$actual == imp.params$missing.fill)
  }
  if (sum(missing.index) > 0) {
    missing.unknown <-
      data.frame(dt=history$dt[missing.index],
                 event='missing')
  } else missing.unknown <- NULL

  if (!is.null(imp.params$events.impute) &&
      sum(missing.index) > 0) {
    missing.unknown <-
      rbind(subset(missing.unknown,
                   !(dt %in% imp.params$events.impute$dt)),
            subset(imp.params$events.impute,
                   dt <= max(history$dt) &
                   dt >= min(history$dt)))
  } else if (!is.null(imp.params$events.impute)) {
    missing.unknown <- subset(imp.params$events.impute,
                              dt <= max(history$dt) &
                              dt >= min(history$dt))
  }
  missing.unknown <- missing.unknown[order(missing.unknown$dt), ]
  dt_format <- get(dt.format)
  prev.index <-
    dt_format(c(NA, missing.unknown$dt[-c(nrow(missing.unknown))]))
  next.index <-
    dt_format(c(missing.unknown$dt[-c(1)], NA))

  max.train.end <- dt_format(missing.unknown$dt)
  min.train.end <-dt_format(max.train.end)
  train.end <- min(missing.unknown$dt)

  for (i in 1:nrow(missing.unknown)) {

    if (!(is.na(next.index[i]))) {

      if (abs(next.index[i] - missing.unknown$dt[i]) <=
          min(imp.params$periods)) {
        max.train.end[i] <- NA
      } else {
        max.train.end[is.na(max.train.end)] <- missing.unknown$dt[i]
      }

    } else {
      max.train.end[is.na(max.train.end)] <- missing.unknown$dt[i]
    }

    if (!(is.na(prev.index[i]))) {

      if (abs(prev.index[i] - missing.unknown$dt[i]) >=
          min(imp.params$periods)) {
        train.end <- missing.unknown$dt[i]
      } else {
        min.train.end[i] <- train.end
      }

    }
  }

  missing.unknown$impute.from.future <-
    ((min.train.end - min(history$dt)) < imp.params$length.history)
  missing.unknown$train.end <- ifelse(missing.unknown$impute.from.future,
                                      max.train.end,
                                      min.train.end)
  missing.unknown$train.end <- dt_format(missing.unknown$train.end)
  return(missing.unknown)

}

#' Function to impute missing or user specified events to impute
#'
#' @param history: timeseries history to impute missing values for
#' @param imputation.parameters: an ImputationParameters object giving the
#'  parameters used in imputation
#' @param x.reg: the matrix of external regressors to be used in imputing
#'   missing values
#' @return: A list containing the history with imputed missing
#'   values, and statistics summarizing any user specified events to impute
#' @export
ImputeMissingOrUnknown <- function(history, imputation.parameters,
                                   dt.format=".AsPOSIXlt",
                                   x.reg=NULL) {
  impute.train.indices <- GetImputationTrainingIndices(history,
                                                       imputation.parameters,
                                                       dt.format=dt.format)

  # For missing dts with not enough history to impute using the past,
  #  impute using the future and do these last after the other missing
  #  values have been imputed
  impute.train.indices.last <- subset(impute.train.indices,
                                      impute.from.future)
  impute.train.indices.last <-
    impute.train.indices.last[order(impute.train.indices.last$dt,
                                    decreasing=TRUE),]

  impute.train.indices.first <- subset(impute.train.indices,
                                       !(impute.from.future))
  impute.train.indices.first <-
    impute.train.indices.first[order(impute.train.indices.first$dt),]

  impute.train.indices <- rbind(impute.train.indices.first,
                                impute.train.indices.last)
  impute.train.indices.first <- impute.train.indices.last <- NULL

  if (is.null(history$actual.upper)) {
    history$actual.upper <- history$actual
  }
  if (is.null(history$actual.lower)) {
    history$actual.lower <- history$actual
  }
  history.imputed <- history

  ForecastMethod <- get(imp.params$method)

  for (trainend in unique(impute.train.indices$train.end)) {

    impute.train.df <- subset(impute.train.indices, train.end == trainend)

    if (all(!impute.train.df$impute.from.future)) {
      history.impute <- subset(history, dt < trainend &
                               dt >= (trainend -
                                      imputation.parameters$length.history))

      impute.dts <- list(begin.dt=min(impute.train.df$dt),
                           end.dt=max(impute.train.df$dt))
      if (!is.null(x.reg)) {
        x.reg.impute <-
          subset(x.reg, dt < trainend &
                 dt >= (trainend - imputation.parameters$length.history))
        x.reg.impute <- subset(x.reg.impute, select=-c(dt))
        x.future.impute <-
          subset(x.reg, dt >= min(impute.train.df$dt) &
                 dt <= max(impute.train.df$dt))
        x.future.impute <- subset(x.future.impute, select=-c(dt))

      }
      # Forecast each cluster of missing data using recent history
      imputed.history <- ForecastMethod(history.impute, impute.dts,
                                        period=imputation.parameters$periods,
                                        periods.agg=c(1),
                                        pred.level=
                                          imputation.parameters$pred.level,
                                        transform=
                                          imputation.parameters$transform,
                                        x.reg=x.reg.impute,
                                        x.future=x.future.impute)$forecast

      history.imputed$actual[history.imputed$dt %in% impute.train.df$dt] <-
        imputed.history$forecast[imputed.history$dt %in%
                                 impute.train.df$dt]
      history.imputed$actual.lower[history.imputed$dt %in%
                                   impute.train.df$dt] <-
        imputed.history$forecast.lower[imputed.history$dt %in%
                                       impute.train.df$dt]
      history.imputed$actual.upper[history.imputed$dt %in%
                                     impute.train.df$dt] <-
        imputed.history$forecast.upper[imputed.history$dt %in%
                                         impute.train.df$dt]


    } else if (all(impute.train.df$impute.from.future)) {
      # For missing dts without enough prior history to impute use future
      history.impute <- subset(history, dt > trainend &
                               dt <= (trainend +
                                      imputation.parameters$length.history))
      history.impute <- history.impute[order(history.impute$dt,
                                             decreasing=TRUE),]

      impute.dts <- list(begin.dt=min(impute.train.df$dt),
                         end.dt=max(impute.train.df$dt))
      if (!is.null(x.reg)) {
        x.reg.impute <-
          subset(x.reg, dt > trainend &
                   dt <= (trainend + imputation.parameters$length.history))
        x.reg.impute <- x.reg.impute[order(x.reg.impute$dt, decreasing=T),]
        x.reg.impute <- subset(x.reg.impute, select=-c(dt))
        x.future.impute <-
          subset(x.reg, dt >= min(impute.train.df$dt) &
                 dt <= max(impute.train.df$dt))
        x.future.impute <- subset(x.future.impute, select=-c(dt))
      }
      # Forecast early missing or unknown values using the future values
      imputed.history <- ForecastMethod(history.impute, impute.dts,
                                        period=imputation.parameters$periods,
                                        periods.agg=c(1),
                                        pred.level=
                                          imputation.parameters$pred.level,
                                        transform=
                                          imputation.parameters$transform,
                                        x.reg=x.reg.impute,
                                        x.future=x.future.impute)$forecast

      history.imputed$actual[history.imputed$dt %in% impute.train.df$dt] <-
        imputed.history$forecast[imputed.history$dt %in%
                                   impute.train.df$dt]
      history.imputed$actual.lower[history.imputed$dt %in%
                                     impute.train.df$dt] <-
        imputed.history$forecast.lower[imputed.history$dt %in%
                                         impute.train.df$dt]
      history.imputed$actual.upper[history.imputed$dt %in%
                                     impute.train.df$dt] <-
        imputed.history$forecast.upper[imputed.history$dt %in%
                                         impute.train.df$dt]


    }
  }
  events.impute <- imputation.parameters$events.impute
  if (!is.null(events.impute)) {
    events.impute <- subset(events.impute, dt >= min(history$dt) &
                              dt <= max(history$dt))
    events.impute$actual <- history$actual[history$dt %in% events.impute$dt]
    events.impute$imputed <-
      history.imputed$actual[history.imputed$dt %in% events.impute$dt]

    events.impute$imputed.lower <-
      history.imputed$actual.lower[history.imputed$dt %in% events.impute$dt]

    events.impute$imputed.upper <-
      history.imputed$actual.upper[history.imputed$dt %in% events.impute$dt]
    events.impute$impact <- events.impute$actual / events.impute$imputed
    events.impute$impact.lower <-
      events.impute$actual / events.impute$imputed.upper
    events.impute$impact.upper <-
      events.impute$actual / events.impute$imputed.lower
  }

  return(list(imputed=history.imputed, events.impact=events.impute))
}

#' Function to impute missing or user specified events to impute
#'
#' @param history: timeseries history to impute missing values for
#' @param imputation.parameters: an ImputationParameters object giving the
#'  parameters used in imputation
#' @param x.reg: the matrix of external regressors to be used in imputing
#'   missing values
#' @return: A list containing the filled history, missingness analysis results,
#' and the history with imputed missing values (if specified to impute)
#' @export
FillAnalyzeImpute <- function(history, imputation.parameters,
                              dt.units='days', dt.format=".AsPOSIXlt",
                              x.reg=NULL) {

  fill.default <- imputation.parameters$missing.fill
  history.filled <- FillMissingDts(history,
                                   fill.default=fill.default,
                                   dt.units=dt.units,
                                   dt.format=dt.format)
  history.filled$actual.lower <- history.filled$actual
  history.filled$actual.upper <- history.filled$actual

  if (imputation.parameters$analyze.missing) {
    missingness.analysis <- AnalyzeMissingness(history.filled, x.reg)
  } else {
    missingness.analysis <- NULL
  }

  if (imputation.parameters$impute.missing) {
    history.imputed <- ImputeMissingOrUnknown(history.filled,
                                              imputation.parameters,
                                              dt.format=dt.format,
                                              x.reg)
    events.impact <- history.imputed$events.impact
    history.imputed <- history.imputed$imputed
  } else {
    history.imputed <- events.impact <- NULL
  }

  return(list(filled=history.filled, imputed=history.imputed,
              missingness.analysis=missingness.analysis))
}


