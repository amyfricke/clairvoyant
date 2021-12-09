
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
#' Function to get forecast error for the final daily and weekly forecasts and
#'  booleans of whether the prediction interval covers the actual value
#'
#' @param daily.history: the full daily history (including dates that may be in
#'   the forecast
#' @param daily.forecast: the final daily forecast
#' @param weekly.forecast: the final weekly forecast
#' @param get.daily.pred.int: boolean indicating if there are daily prediction
#'   intervals to be evaluated
#' @return: list of dataframes containing forecast evaluation data
#' @export
EvaluateClairvoyant <- function(unaggregated.history, unaggregated.forecast,
                                aggregation.parameters=NULL,
                                aggregated.histories=NULL,
                                aggregated.forecasts=NULL,
                                dt.format='.AsPOSIXlt',
                                longest.period=364,
                                get.unagg.pred.int=TRUE) {

  unaggregated.eval <- subset(unaggregated.history,
                              dt %in% unaggregated.forecast$dt)

  periods.agg <- aggregation.parameters$periods.agg
  agg.fun <- aggregation.parameters$agg.fun

  if (nrow(unaggregated.eval) < min(periods.agg)) return(NULL)

  unaggregated.eval$test <- 1

  if (!is.null(aggregated.histories) && !is.null(aggregated.forecasts)) {
    n.agg <- floor(nrow(unaggregated.eval) / max(periods.agg)) * max(periods.agg)
    aggregated.evals <-
        AggregateToLongest(unaggregated.eval[1:n.agg,],
                           periods.agg=periods.agg,
                           agg.fun=agg.fun,
                           cols.agg=c("actual", "test"),
                           dt.format=dt.format)
    stopifnot(all(aggregated.evals[paste(max(periods.agg))]$test ==
                    max(periods.agg)))

    for (i in 1:length(aggregated.evals)) {

      aggregated.evals[[i]]$test <- NULL
      aggregated.evals[[i]] <- merge(aggregated.evals[[i]],
                                     aggregated.forecasts[[i]], by=c('dt'))
      aggregated.evals[[i]]$error <-
       ((aggregated.evals[[i]]$forecast - aggregated.evals[[i]]$actual) /
        aggregated.evals[[i]]$actual)

      aggregated.evals[[i]]$abs.error <- abs(aggregated.evals[[i]]$error)
      aggregated.evals[[i]]$prediction.interval.coverage <-
        (aggregated.evals[[i]]$actual >= aggregated.evals[[i]]$forecast.lower &
         aggregated.evals[[i]]$actual <= aggregated.evals[[i]]$forecast.upper)

      if (!is.null(aggregated.evals[[i]]$range.forecast.lower)) {

        aggregated.evals[[i]]$range.forecast.coverage <-
          (aggregated.evals[[i]]$actual >=
           aggregated.evals[[i]]$range.forecast.lower &
           aggregated.evals[[i]]$actual <=
           aggregated.evals[[i]]$range.forecast.upper)

        aggregated.evals[[i]]$range.uncertainty.coverage <-
          (aggregated.evals[[i]]$actual >=
           aggregated.evals[[i]]$range.uncertainty.lower &
            aggregated.evals[[i]]$actual <=
           aggregated.evals[[i]]$range.uncertainty.upper)
      }

    }
  } else {
    aggregated.evals <- NULL
  }

  unaggregated.eval$test <- NULL
  unaggregated.eval <- merge(unaggregated.eval,
                             unaggregated.forecast, by=c("dt"))
  unaggregated.eval$error <-
      ((unaggregated.eval$forecast - unaggregated.eval$actual) /
       unaggregated.eval$actual)
  unaggregated.eval$abs.error <- abs(unaggregated.eval$error)

    unaggregated.eval$prediction.interval.coverage <-
      (unaggregated.eval$actual >= unaggregated.eval$forecast.lower &
       unaggregated.eval$actual <= unaggregated.eval$forecast.upper)
  if (!is.null(unaggregated.eval$range.forecast.lower)) {
    unaggregated.eval$range.forecast.coverage <-
      (unaggregated.eval$actual >=
       unaggregated.eval$range.forecast.lower &
       unaggregated.eval$actual <= unaggregated.eval$range.forecast.upper)

    unaggregated.eval$range.uncertainty.coverage <-
      (unaggregated.eval$actual >=
       unaggregated.eval$range.uncertainty.lower &
       unaggregated.eval$actual <= unaggregated.eval$range.uncertainty.upper)
  }
  return(list(unaggregated=unaggregated.eval, aggregated=aggregated.evals))
}

#' Plot the resulting ensemble of models and the consensus forecast
#'
#' @param object: clairvoyant forecast object to plot ensemble
#' @return: NULL< plot producted
#' @export
PlotForecastEnsemble <- function(object, dt.format='as.POSIXct') {

  if (!is.null(object$training$aggregated)) {
    periods.agg <- as.numeric(names(object$training$aggregated))
    training <- object$training$aggregated[[paste(max(periods.agg))]]

  } else if (!is.null(object$training$unaggregated$imputed)) {
    training <- object$training$unaggregated$imputed
  } else {
    training <- object$training$unaggregated$raw
  }
  training$type <- 'training'
  type.levels <- c('training')
  dt_format <- get(dt.format)
  all <- training
  forecast.list <- object$forecast$ensemble
  forecast.list[['consensus']] <-
    object$forecast$consensus[c('dt', 'forecast', 'forecast.lower',
                                'forecast.upper')]
  for (model in ls(forecast.list)) {
    model.forecast <- forecast.list[[model]]
    colnames(model.forecast) <- c('dt', 'actual', 'actual.lower',
                                  'actual.upper')
    model.forecast$type <- paste(model, 'forecast')
    type.levels <- c(type.levels, paste(model, 'forecast'))
    all <- rbind(all, model.forecast)
  }
  all$type <- factor(all$type, levels=type.levels)

  min.year <- as.numeric(substr(min(all$dt), 1, 4))
  min.year <- min.year + (min(all$dt) >=
                          dt_format(paste(min.year, "-07-01", sep="")))
  min.year <- dt_format(paste(min.year, "-01-01", sep=""))
  max.year <- as.numeric(substr(max(all$dt), 1, 4))
  max.year <- max.year + (max(all$dt) >=
                          dt_format(paste(max.year, "-07-01", sep="")))
  max.year <- dt_format(paste(max.year, "-01-01", sep=""))
  dtbreaks <- seq(min.year, max.year, by="year")
  all$line.width <- 1 + (all$type == "consensus")
  all$dt <- dt_format(all$dt)

  title.name <-
      "Training and forecast model ensemble"
  p <- ggplot2::ggplot(data=all, ggplot2::aes(x=dt, y=actual, color=type)) +
       ggplot2::geom_line(size=all$line.width)
  p <- p + ggplot2::scale_x_datetime(breaks=dtbreaks,
                                     labels=scales::date_format("%Y-%m-%d")) +
      ggplot2::xlab("") + ggplot2::ylab("") +
      ggplot2::ggtitle(title.name) +
      ggplot2::theme(axis.text.x=ggplot2::element_text(size=15, face="bold"),
                     axis.text.y=ggplot2::element_text(size=15, face="bold"),
                     plot.title=ggplot2::element_text(size=20,
                                                      vjust=2, face="bold"),
                     legend.title=ggplot2::element_blank(),
                     legend.text=ggplot2::element_text(size=12),
                     legend.background=
                     ggplot2::element_rect(fill = "transparent"),
                     legend.position="bottom",
                     legend.key=ggplot2::element_rect(fill = "transparent"),
                      legend.key.width=grid::unit(1.5, "cm")) +
      ggplot2::theme(panel.border=ggplot2::element_blank(),
                     panel.background=ggplot2::element_blank(),
                     axis.line=ggplot2::element_line(color="grey"),
                     panel.grid.major=
                     ggplot2::element_line(color="#CCCCCC",
                                           linetype="dashed"),
                     panel.grid.minor=ggplot2::element_line(color="#DDDDDD",
                                                            linetype="dashed"))
  print(p)
}

#' Plot the resulting ensemble of models and the consensus forecast
#'
#' @param object: clairvoyant forecast object to plot ensemble
#' @return: NULL< plot producted
#' @export
PlotEnsembleError <- function(object, dt.format='as.POSIXct') {


  dt_format <- get(dt.format)
  error.list <- object$evaluation$ensemble
  cols.incl <- colnames(error.list[[1]])
  if (length(object$evaluation$aggregated) > 0) {
    period.agg <- max(as.numeric(ls(object$evaluation$aggregated)))
    error.list[['consensus']] <-
      object$evaluation$aggregated[[paste(period.agg)]][cols.incl]
  } else {
    error.list[['consensus']] <- object$evaluation$unaggregated[cols.incl]
  }
  type.levels <- NULL
  all <- NULL
  for (model in ls(error.list)) {
    model.error <- error.list[[model]]
    colnames(model.error) <- cols.incl
    model.error$type <- paste(model, 'error')
    type.levels <- c(type.levels, paste(model, 'error'))
    all <- rbind(all, model.error)
  }
  all$type <- factor(all$type, levels=type.levels)

  min.year <- as.numeric(substr(min(all$dt), 1, 4))
  min.year <- min.year + (min(all$dt) >=
                            dt_format(paste(min.year, "-07-01", sep="")))
  min.year <- dt_format(paste(min.year, "-01-01", sep=""))
  max.year <- as.numeric(substr(max(all$dt), 1, 4))
  max.year <- max.year + (max(all$dt) >=
                            dt_format(paste(max.year, "-07-01", sep="")))
  max.year <- dt_format(paste(max.year, "-01-01", sep=""))
  dtbreaks <- seq(min.year, max.year, by="year")
  all$line.width <- 1 + (all$type == "consensus")
  all$dt <- dt_format(all$dt)

  title.name <-
    "Training and forecast model ensemble"
  p <- ggplot2::ggplot(data=all, ggplot2::aes(x=dt, y=error, color=type)) +
    ggplot2::geom_line(size=all$line.width)
  p <- p + ggplot2::scale_x_datetime(breaks=dtbreaks,
                                     labels=scales::date_format("%Y-%m-%d")) +
    ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::ggtitle(title.name) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(size=15, face="bold"),
                   axis.text.y=ggplot2::element_text(size=15, face="bold"),
                   plot.title=ggplot2::element_text(size=20,
                                                    vjust=2, face="bold"),
                   legend.title=ggplot2::element_blank(),
                   legend.text=ggplot2::element_text(size=12),
                   legend.background=
                     ggplot2::element_rect(fill = "transparent"),
                   legend.position="bottom",
                   legend.key=ggplot2::element_rect(fill = "transparent"),
                   legend.key.width=grid::unit(1.5, "cm")) +
    ggplot2::theme(panel.border=ggplot2::element_blank(),
                   panel.background=ggplot2::element_blank(),
                   axis.line=ggplot2::element_line(color="grey"),
                   panel.grid.major=
                     ggplot2::element_line(color="#CCCCCC",
                                           linetype="dashed"),
                   panel.grid.minor=ggplot2::element_line(color="#DDDDDD",
                                                          linetype="dashed"))
  print(p)
}


#' Function to plot the training, forecast and actual values
#'
#' @param object: clairvoyant forecast object with steps in process
#' @param aggregated: whether to plot aggregated forecast
#' @param freq: frequency forecast should be aggregated to for the plot
#' @param plot.range: whether the range forecast should be plotted
#' @param short.term: boolean indicating if should be truncated to short term
#' @param n.history: number of training dts to include in short.term plot
#' @param n.forecast: number of forecast dts to include in short.term plot
#' @param freq: frequency of error (weekly aggregate or daily)
#' @param outl.dir: output directory
#' @return: none (plot producted)
#' @export
PlotTrainingActualAndForecast <- function(object,
                                          aggregated=FALSE,
                                          freq=7,
                                          plot.range=TRUE,
                                          short.term=FALSE,
                                          dt.format='.AsPOSIXlt',
                                          dt.units='days',
                                          n.history=28, n.forecast=14) {
  if (aggregated) {
    training <- object$training$aggregated[[paste(freq)]]
    forecast <- object$forecast$aggregated$final[[paste(freq)]]
    evaluation <- object$evaluation$aggregated[[paste(freq)]]
  } else {
    training <- object$training$unaggregated$raw
    forecast <- object$forecast$final
    evaluation <- object$evaluation$unaggregated
  }
  dt_format <- get(dt.format)
  colnames(forecast) <- c('dt', 'value', 'value.lower', 'value.upper',
                          'range.value.lower', 'range.value.upper',
                          'range.uncertainty.lower', 'range.uncertainty.upper')
  training$range.value.lower <- NA
  training$range.value.upper <- NA
  training$range.uncertainty.lower <- NA
  training$range.uncertainty.upper <- NA
  if (is.null(training$actual.lower)) {
    training$actual.lower <- NA
  }
  if (is.null(training$actual.upper)) {
    training$actual.upper <- NA
  }
  if (short.term) {
    training <- tail(training, n.history)
    forecast <- head(forecast, n.forecast)
    if (!is.null(evaluation)) {
      evaluation <- head(evaluation, n.forecast)
      evaluation$dt <- dt_format(evaluation$dt)
    }
  }

  colnames(training) <- c('dt', 'value', 'value.lower', 'value.upper',
                          'range.value.lower', 'range.value.upper',
                          'range.uncertainty.lower', 'range.uncertainty.upper')
  training$type <- "training"
  color.idx <- c(3, 1)
  
  forecast$type <- "forecast"
  all <- rbind(training, forecast)
  all$type <- factor(all$type, levels=c("training", "forecast"))
  all$dt <- dt_format(all$dt)
  all$idx <- c(1:nrow(all))
  forecast$dt <- dt_format(forecast$dt)
  forecast$idx  <- all$idx[all$dt %in% forecast$dt]
  
  if (!is.null(evaluation)) {
    evaluation$type <- "actuals"
    color.idx <- c(2, 1, 3)
    evaluation$idx <- all$idx[all$dt %in% evaluation$dt]
  
  }
  
  if (!short.term) {
    # decide x-axis start year, end year and line width
    min.year <- as.numeric(substr(dt_format(min(all$dt)), 1, 4))
    min.year <- min.year + (dt_format(min(all$dt)) >=
                            dt_format(paste0(min.year, "-07-01")))
    min.year <- dt_format(paste0(min.year, "-01-01"))
    max.year <- as.numeric(substr(dt_format(max(all$dt)), 1, 4))
    max.year <- max.year + (dt_format(max(all$dt)) >=
                            dt_format(paste0(max.year, "-07-01")))
    max.year <- dt_format(paste0(max.year, "-01-01"))
    dtbreaks <- seq(min.year, max.year, by="year")
  } else {
    dtbreaks <- seq(min(training$dt), 
                    max(forecast$dt), by=paste(freq, dt.units))
  }
  
  line.width <- 0.5 + (aggregated) * 0.5

  freq.tmp <- ifelse(aggregated, paste("Aggregated to period", freq),
                     "Unaggregated")
  
  dtbreaks.idx <- NULL
  for (brk in dtbreaks) {
    dtbreaks.idx <- c(dtbreaks.idx,
                      which.min(abs(all$dt - dt_format(brk))))
  }
  
  
  substr(freq.tmp, 1, 1) <- toupper(substr(freq.tmp, 1, 1))
  title.name <- paste(freq.tmp, "Training and Forecast")
  
  
  p <- ggplot2::ggplot(data=all, ggplot2::aes(x=idx, y=value, color=type))
  
  p <- p + ggplot2::geom_ribbon(data=forecast,
                                ggplot2::aes(x=idx, ,
                                             ymax=range.uncertainty.upper,
                                             ymin=range.uncertainty.lower,
                                             color=NULL),
                                 alpha=0.2, color='transparent',
                                size=line.width, fill='red')
   
  p <- p + ggplot2::geom_ribbon(data=forecast,
                                ggplot2::aes(x=idx,
                                             ymax=range.value.upper,
                                             ymin=range.value.lower,
                                             color=NULL),
                                alpha=0.4, color='transparent',
                                size=line.width, fill='red')
  p <- p + ggplot2::geom_line(size=line.width)
  if (!is.null(evaluation)) {
    p <- p + ggplot2::geom_line(data=evaluation, ggplot2::aes(
      x=idx, y=actual),
      size=line.width, alpha=0.6)
  }
  p <- p +
    ggplot2::scale_color_manual(values=scales::hue_pal()(3)[color.idx],
                                name="") +
    ggplot2::scale_x_continuous(breaks=dtbreaks.idx,
                                labels=dtbreaks) +
    
    ggplot2::coord_cartesian(ylim=c(0.9 * min(training$value),
                                    2* max(forecast$value))) +
    ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::ggtitle(title.name) +
    ggplot2::theme(axis.text.x=ggplot2::element_text(size=15, face="bold"),
                   axis.text.y=ggplot2::element_text(size=15, face="bold"),
                   plot.title=ggplot2::element_text(size=20,
                                                    vjust=2, face="bold"),
                   legend.position="bottom",
                   legend.text=ggplot2::element_text(size=15),
                   legend.background=
                     ggplot2::element_rect(fill = "transparent"),
                   legend.key=ggplot2::element_rect(fill = "transparent"),
                   legend.key.width=grid::unit(1.5, "cm")) +
    ggplot2::theme(panel.border=ggplot2::element_blank(),
                   panel.background=ggplot2::element_blank(),
                   axis.line=ggplot2::element_line(color="grey"),
                   panel.grid.major=
                     ggplot2::element_line(color="#CCCCCC",
                                           linetype="dashed"),
                   panel.grid.minor=
                     ggplot2::element_line(color="#DDDDDD",
                                           linetype="dashed"))
  print(p)

}

#' Function to plot error for two different forecast specifications
#'
#' @param object1: first clairvoyant forecast object to compare
#' @param object2: second clairvoyant forecast object to compare
#' @param label1: label for the first forecast specification
#' @param label2: label for the 2nd forecast specification
#' @param freq: the frequency of the error to plot (daily or weekly)
#' @param short.term: whether to truncate to the first n.error measurements
#' @param n.error: number of measurements to plot if truncating
#' @param log.scale: boolean of whether to plot the error on the log scale
#' @param dts.annotations: dataframe of events to annotate the plot
#' @return: null, plot is produced
#' @export
PlotErrorComparison <- function(object1, object2=NULL, label1='', label2='',
                                aggregated=FALSE,
                                freq=7, dt.format='as.POSIXct',
                                short.term=T, n.error=14,
                                log.scale=T, dts.annotations=NULL) {
  if (aggregated) {
    error1 <- object1$evaluation$aggregated[[paste(freq)]]
    if (!is.null(object2)) {
      error2 <- object2$evaluation$aggregated[[paste(freq)]]
      error2$type <- label2
    }
  } else {
    error1 <- object1$evaluation$unaggregated
    if (!is.null(object2)) {
      error2 <- object2$evaluation$unaggregated
      error2$type <- label2
    } else error2 <- NULL
  }

  error1$type <- label1


  if (short.term) {
    error1 <- head(error1, n.error)
    error2 <- head(error2, n.error)
  }
  if (!is.null(dts.annotations)) {
    if (class(dts.annotations) == 'character' &&
        dts.annotations == 'DOW') {
      dts.annotations <- data.frame(dt=error1$dt,
                                    event=format(dt_format(error1$dt), '%a'))
    }
    dts.annotations <-
        subset(dts.annotations, dt %in% c(error1$dt, error2$dt))
    dts.annotations <- dts.annotations[order(dts.annotations$dt), ]
    dts.annotations$error <-
        subset(error1, dt %in% dts.annotations$dt)$error

  }

  error.all <- rbind(error1, error2)
  if (!is.null(error2)) {
    error.all$type <- factor(error.all$type, levels=c(label1, label2))
  } else error.all$type <- factor(error.all$type, levels=c(label1))

  freq.tmp <- freq
  substr(freq.tmp, 1, 1) <- toupper(substr(freq.tmp, 1, 1))
  title.name <-
      paste(freq.tmp, "Relative Forecast Error:",
           "(forecast - truth) / truth \n")
   if (label1 != '' || label2 != '') {
     title.name <- paste(title.name, label1, 'vs', label2)
  }
  p <- ggplot2::ggplot(data = error.all,
                       ggplot2::aes(x=dt, y=error, color=type)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey", size = 1) +
    ggplot2::geom_line(size = 1) +
    ggplot2::scale_x_datetime(labels = scales::date_format("%Y-%m-%d")) +
    ggplot2::xlab("") + ggplot2::ylab("") +
    ggplot2::ggtitle(title.name) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 15, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 15, face = "bold"),
      plot.title = ggplot2::element_text(size = 20, vjust = 2, face = "bold"),
      legend.position = 'bottom',
      legend.text = ggplot2::element_text(size = 15),
      legend.title = ggplot2::element_blank(),
      legend.background = ggplot2::element_rect(fill = "transparent"),
      legend.key = ggplot2::element_rect(fill = "transparent"),
      legend.key.width = grid::unit(1.5, "cm")
    ) +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "grey"),
      panel.grid.major = ggplot2::element_line(color = "#CCCCCC",
                                               linetype = "dashed"),
      panel.grid.minor = ggplot2::element_line(color = "#DDDDDD",
                                               linetype = "dashed")
    )
  if (label1 == '' && label2 == '') {
    p <- p + ggplot2::theme(legend.position = 'none')
  }
  if (!is.null(dts.annotations)) {
    p <- p +
        ggplot2::geom_text(data=dts.annotations, size=6,
                           vjust=-1 ^ c(1:nrow(dts.annotations)) * 0.5,
                           ggplot2::aes(x=dt, y=error,
                                        label=event, color=label1))
  }
  if (log.scale) {
    custom_log_y_trans <- function() {
    scales::trans_new(
        "custom_log_y",
         transform = function (x) (sign(x) * log(abs(x) + 1, 10)),
         inverse = function (y) (sign(y) * (10 ^ (abs(y)) - 1)),
         domain = c(-Inf, Inf))
    }
    errs.na <- error.all$error[!is.na(error.all$error)]
    if (min(errs.na) < 0) {
      y.min.less0 <- min(c(-1e-1, min(errs.na[errs.na < 0])))
      y.max.less0 <- min(c(-1e-1, max(errs.na[errs.na < 0])))
      y.breaks.less0 <-
          c(round(log(-y.min.less0, 10)):round(log(-y.max.less0, 10)))
      y.breaks.less0 <- - 10 ^ (y.breaks.less0)
    } else y.breaks.less0 <- NULL
    if (max(errs.na) > 0) {
      y.min.greater0 <- max(c(1e-1, min(error.all$error[error.all$error > 0])))
      y.max.greater0 <- max(c(1e-1, max(error.all$error[error.all$error > 0])))
      y.breaks.great0 <-
          c(round(log(y.min.greater0, 10)):round(log(y.max.greater0, 10)))
      y.breaks.great0 <- 10 ^ (y.breaks.great0)
    } else y.breaks.great0 <- NULL
    y.breaks <- c(y.breaks.less0, y.breaks.great0)

    p <- p + ggplot2::scale_y_continuous(
        trans=custom_log_y_trans(),
         limits=c(min(error.all$error, 0),
         max(error.all$error, 0)), breaks=y.breaks)
  } else {
    p <- p + ylim(min(error.all$error, 0), max(error.all$error, 0))
  }
  print(p)

}
