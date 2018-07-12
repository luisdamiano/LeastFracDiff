#' Finds the least fractional differentiation order subject to fun returning less than threshold.
#'
#' @param x A time series. Either a numberic vector or an `xts` object.
#' @param fun A function taking as only argument the time series and returning a p-value.
#' @param threshold A numerical threshold for the p-value.
#' @param dFrom Lower bound for `d`.
#' @param dTo Upper bound for `d`.
#' @param dBy Increment for each in the grid.
#'
#' @return The least value of `d`.
#' @export
#'
#' @examples \dontrun{min_d(x, adf, 0.05)}
min_d  <- function(x, fun, threshold, dFrom = 0, dTo = 1, dBy = 0.1) {
  dSeq  <- seq(from = dFrom, to = dTo, by = dBy)
  for (i in 2:length(dSeq)) {
    if (fun(diff_series(x, dSeq[i])) < threshold)
      return(dSeq[i])
  }
}

#' Returns a grid for different values of the fractional differentiation order parameter `d`.
#'
#' @param x A time series. Either a numberic vector or an `xts` object.
#' @param fun A function taking as only argument the time series and returning a p-value.
#' @param dFrom Lower bound for `d`.
#' @param dTo Upper bound for `d`.
#' @param dBy Increment for each in the grid.
#'
#' @return A grid with four columns: the value of the fractional differentiation order `d`, the p-value returned by `fun`, the Pearson correlation coefficient between the original and the transform series, and the Spearman counterpart.
#' @export
#' @importFrom stats cor
#'
#' @examples \dontrun{grid_d(x, adf)}
grid_d <- function(x, fun, dFrom = 0, dTo = 1, dBy = 0.1) {
  dSeq    <- seq(from = dFrom, to = dTo, by = dBy)
  out     <- sapply(
    dSeq,
    function(d) {
      dx <- diff_series(x, d)
      c(
        pValue              = fun(dx),
        pearsonCor          = cor(x, dx),
        spearmanCorrelation = cor(x, dx, method = "spearman")
      )
    }
  )

  cbind(d = dSeq, t(out))
}

#' Finds the least fractional differentiation order subject to fun returning less than threshold for a rolling window.
#'
#' @param x A time series. Either a numberic vector or an `xts` object.
#' @param width The number of observations included in each window.
#' @param by The number of observations between each run.
#' @param ... Other parameters to be passed to zoo::rollapply and min_d.
#'
#' @return An `xts` object with the time series of `d`.
#' @importFrom zoo rollapply
#' @export
#'
#' @examples \dontrun{roll_min_d(x, width = 252, by = 1)}
roll_min_d  <- function(x, width, by, ...) {
  zoo::rollapply(
    data = x, FUN = min_d ,
    width = width, by = by,
    ...
  )
}

#' Applies a fractional differentiation to a time series.
#'
#' @param x A time series. Either a numberic vector or an `xts` object.
#' @param d Order of the fractional differentiation.
#'
#' @return A time series object.
#' @export
#'
#' @examples \dontrun{diff_series(x, 0.5)}
diff_series <- function(x, d) {
  fracdiff::diffseries(x, d)
}

#' Return the p-value of the Augmented Dickey-Fuller test.
#'
#' @param x A time series. Either a numberic vector or an `xts` object.
#'
#' @return A p-value.
#' @export
#'
#' @examples \dontrun{adf(x)}
adf <- function(x) {
  suppressWarnings(tseries::adf.test(x)$p.value)
}

#' Return the p-value of the Kwiatkowski–Phillips–Schmidt–Shin test.
#'
#' @param x A time series. Either a numberic vector or an `xts` object.
#'
#' @return A p-value.
#' @export
#'
#' @examples \dontrun{kpss(x)}
kpss <- function(x) {
  suppressWarnings(tseries::kpss.test(x)$p.value)
}
