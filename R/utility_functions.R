#' Internal method for formatting datetimes with the origin
#' always supplied
#' @return datetime object with origin so its compatible with
#'  numeric datetime objects
.AsPOSIXlt <- function(dt) {
  return(as.POSIXct(dt, origin='1970-01-01 00:00:00', tz='GMT'))
}

#' Internal method for getting begin.date and end.date. Dates
#' are adjusted so that they span an integer number of weeks.
#'
#' @param begin.dt Candidate start datetime
#' @param end.dt Candidate end datetime (inclusive)
#' @param dt.units String giving the units of the datetimes
#' @param fix.date.option whether to cut data to the right or
#'   left ('end.later' vs 'begin.later')
#' @param need.integer.period whether an integer number of
#'   periods are needed (if aggregating)
#' @param period length of the period for which
#' @return A list containing a begin.date and an end.date
#' @noRd
.GetBeginEndDts <- function(begin.dt, end.dt,
                            dt.units='days',
                            dt.format='.AsPOSIXlt',
                            fix.dt.option,
                            need.integer.period=FALSE, period=7) {
  stopifnot((!is.null(begin.dt) && !is.null(end.dt)))

  if (!need.integer.period) {
    return(list(begin.dt=begin.dt, end.dt=end.dt))
  }
  dt_units <- get(dt.units)
  dt.unit <- substr(dt.units, 1, (nchar(dt.units) - 1))
  dt_format <- get(dt.format)

  # Make sure training data covers an integer # of periods
  # if we are aggregating over seasonal periods
  l.d <- length(seq(begin.dt, end.dt, by=dt.unit))
  seq.d <- rep(c(1:period), ceiling(l.d / period))
  begin.d <- seq.d[1]
  end.d <- seq.d[l.d] + 1
  if (begin.d != end.d) {
    delta <- end.d - begin.d
    if (delta < 0) {
      delta <- delta + period
    }

    if (fix.dt.option == "begin.later") {
      original.begin.dt <- begin.dt
      begin.dt <- begin.dt + dt_units(delta)
    } else {
      original.end.dt <- end.dt
      end.dt <- end.dt + dt_units(period - delta)
    }
  }

  stopifnot(length(seq(begin.dt, end.dt, by=dt.unit)) %% period == 0)

  return(list(begin.dt=dt_format(begin.dt), end.dt=dt_format(end.dt)))
}


#' Internal function to generate a transformed time series
#'
#' @param history vector of values to transform
#' @param transform transformation to apply, options are log, none or box cox
#' @param box.cox.lambda lambda value for the box cox transformation if applied
#' @return A vector of the transformed series
#' @noRd
.TransformTimeSeries <- function(history, transform='box_cox',
                                 box.cox.lambda=1) {

  stopifnot(transform %in% c('log', 'box_cox', 'none'))
  if (transform == 'log' || box.cox.lambda == 0) {
    t.history <- log(history)
  } else if (transform == 'box_cox' && box.cox.lambda != 1) {
    t.history <- (history ^ box.cox.lambda - 1) / box.cox.lambda
  } else {
    t.history <- history
  }
  return(t.history)
}

#' Internal function to back transform a time series
#'
#' @param history vector of values to back transform
#' @param transform transformation applied, options are log, none or Box Cox
#' @param box.cox.lambda lambda value for the Box Cox transformation if applied
#' @return A vector of the back transformed series
#' @noRd
.BackTransformTimeSeries <- function(history, transform='box_cox',
                                     box.cox.lambda=1) {
  # Internal function to back transform a time series
  stopifnot(transform %in% c('log', 'box_cox', 'none'))
  if (transform == 'log' || box.cox.lambda == 0) {
    b.history <- exp(history)
  } else if (transform == 'box_cox' && box.cox.lambda != 1) {
    b.history <- (box.cox.lambda * history + 1) ^ (1 / box.cox.lambda)
  } else {
    b.history <- history
  }
  return(b.history)
}

#' Helper (internal) function to divide vectors that might have 0s at the same
#'  places
#'
#' @param x: vector in the numerator
#' @param y: vector in the denominator
#' @return: x / y
#' @noRd
.Divide0 <- function(x, y) {
  # Internal function to yield x / y for x and y vectors, yields 1 when xi=yi=0
  z <- x / y
  z[which(x == 0 & y == 0)] <- 1
  return(z)
}

#' Function to get lower quartile
#'
#' @param x: vector to get lower quartile of
#' @return: lower quartile of x
#' @export
LowerQuartile <- function(x, na.rm=TRUE) {
  return(quantile(x, probs=c(0.25), na.rm=na.rm))
}

#' Function to get upper quartile
#'
#' @param x: vector to get lower quartile of
#' @return: lower quartile of x
#' @export
UpperQuartile <- function(x, na.rm=TRUE) {
  return(quantile(x, probs=c(0.75), na.rm=na.rm))
}


#' Fix for prophet model regression coefficients extraction
#'
#' @param m: fitted prophet model
#' @return: dataframe containing the regression coefficients
#' @noRd
regressor_coefficients <- function(m){
  if (length(m$extra_regressors) == 0) {
    stop("No extra regressors found.")
  }
  regr_names <- names(m$extra_regressors)
  regr_modes <- unlist(lapply(m$extra_regressors, function(x) x$mode))
  regr_mus <- unlist(lapply(m$extra_regressors, function (x) x$mu))
  regr_stds <- unlist(lapply(m$extra_regressors, function(x) x$std))

  beta_indices <- which(m$train.component.cols[, regr_names, drop = FALSE] == 1,
                        arr.ind = TRUE)[, "row"]
  betas <- m$params$beta[, beta_indices, drop = FALSE]
  # If regressor is additive, multiply by the scale factor to put coefficients 
  # on the original training data scale.
  y_scale_indicator <- matrix(
    data = ifelse(regr_modes == "additive", m$y.scale, 1),
    nrow = nrow(betas),
    ncol = ncol(betas),
    byrow = TRUE
  )
  coefs <- betas * y_scale_indicator  / regr_stds

  percentiles = c((1 - m$interval.width) / 2, 1 - (1 - m$interval.width) / 2)
  bounds <- apply(coefs, 2, stats::quantile, probs = percentiles)

  df <- data.frame(
    regressor = regr_names,
    regressor_mode = regr_modes,
    center = regr_mus,
    coef_lower = bounds[1, ],
    coef = apply(coefs, 2, mean),
    coef_upper = bounds[2, ],
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  return(df)
}
