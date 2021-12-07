#' Helper to find the optimal z alpha score for computing uncertainty regions
#'  for range forecasts - see "Vansteelandt, Stijn, et al. "Ignorance and
#'  uncertainty regions as inferential tools in a sensitivity analysis."
#'  Statistica Sinica (2006): 953-979.
#'
#' @param c.alpha.star: candidate z score for uncertainty region
#' @param fcst.lower: the lower end of the range forecast
#' @param fcst.upper: the upper end of the range forecast
#' @param se.lower: the standard error of the lower end range forecast
#' @param se.lower: the standard error of the upper end range forecast
#' @return The absolute difference between the desired confidence level and
#'  that for the candidate range z score
AdjustedAlphaHelper <- function(c.alpha.star, fcst.lower, fcst.upper,
                                se.lower, se.upper, alpha) {
  a <- -1 * (c.alpha.star)
  b <-  a - ((fcst.upper - fcst.lower) / se.upper)
  value1 <- (pnorm(-1 * a) - pnorm(b))
  c <- -1 * a + ((fcst.upper - fcst.lower) / se.lower)
  value2 <- (pnorm(c) - pnorm(a))
  minimum.value <- min(value1, value2)

  return(abs(minimum.value - (1 - alpha)))

}

#' Function to find the optimal z alpha score for computing uncertainty regions
#'  for range forecasts - see "Vansteelandt, Stijn, et al. "Ignorance and
#'  uncertainty regions as inferential tools in a sensitivity analysis."
#'  Statistica Sinica (2006): 953-979.
#'
#' @param range.fcst.df: dataframe or named vector containing the forecast
#'  range limits and the corresponding standard errors
#' @param alpha: the nominal type I error for the uncertainty regions
#' @return The adjusted Z score given the length of the range forecast
#'  and the standard errors of the limits of the range forecast
#'  @export
GetAdjustedAlpha <- function(range.fcst.df, alpha) {
  names.df <- names(range.fcst.df)
  range.fcst.df <- as.numeric(range.fcst.df)
  names(range.fcst.df) <- names.df
  c.alpha <- -1 * qnorm(seq(alpha / 2, alpha, by=0.0005))

  delta.c.alphas <- lapply(c.alpha, AdjustedAlphaHelper,
                           fcst.lower=range.fcst.df['range.forecast.lower'],
                           fcst.upper=range.fcst.df['range.forecast.upper'],
                           se.lower=range.fcst.df['range.se.lower'],
                           se.upper=range.fcst.df['range.se.lower'],
                           alpha=alpha)
  delta.c.alphas <- as.vector(unlist(delta.c.alphas))

  c.alpha.adjusted <- c.alpha[which.min(delta.c.alphas)]
  return(c.alpha.adjusted)
}
