# ---------------------------------------------------------------------------
# Validation against measured runoff, soil loss and nutrient loss
# ---------------------------------------------------------------------------

#' Validate predictions against measurements
#'
#' Goodness-of-fit statistics for comparing runoff, soil-loss or
#' nutrient-loss predictions with runoff-plot or gauging measurements:
#' bias, percent bias (PBIAS), root mean square error (RMSE), mean absolute
#' error, coefficient of determination, Nash-Sutcliffe efficiency (NSE) and
#' Kling-Gupta efficiency (KGE). For predictive distributions it also reports
#' the continuous ranked probability score (CRPS) and the coverage of central
#' prediction intervals, which show whether the stated uncertainty is honest.
#'
#' @param obs Measured values.
#' @param pred Predictions: a numeric vector or a `dm_predictive` object.
#' @param level Nominal coverage levels of the central prediction intervals.
#' @return An object of class `dm_validation`: a list with `stats` (a named
#'   vector) and `coverage` (observed against nominal coverage, for
#'   predictive distributions).
#' @references Moriasi, D. N. et al. (2007). Model evaluation guidelines for
#'   systematic quantification of accuracy in watershed simulations.
#'   *Transactions of the ASABE*, 50, 885-900.
#'
#'   Gupta, H. V., Kling, H., Yilmaz, K. K. and Martinez, G. F. (2009).
#'   Decomposition of the mean squared error and NSE performance criteria.
#'   *Journal of Hydrology*, 377, 80-91. \doi{10.1016/j.jhydrol.2009.08.003}
#'
#'   Gneiting, T. and Raftery, A. E. (2007). Strictly proper scoring rules,
#'   prediction, and estimation. *Journal of the American Statistical
#'   Association*, 102, 359-378. \doi{10.1198/016214506000001437}
#' @examples
#' ws <- dm_simulate_watershed(n_side = 20, seed = 1)
#' g <- ws$grid
#' soil <- dm_rusle(R = dm_uncertain(g$R, cv = 0.15, lower = 0), K = g$K, LS = g$ls,
#'                  C = g$C, n_draws = 300, seed = 1)
#' dm_validate(g$soil_loss, soil)
#' @export
dm_validate <- function(obs, pred, level = c(0.5, 0.9)) {
  obs <- as.numeric(obs)
  mu <- if (inherits(pred, "dm_predictive")) pred$mean else as.numeric(pred)
  if (length(mu) != length(obs)) stop("'obs' and 'pred' must have the same length.", call. = FALSE)
  ok <- is.finite(obs) & is.finite(mu)
  o <- obs[ok]; p <- mu[ok]
  if (length(o) < 3) stop("At least 3 complete pairs are needed.", call. = FALSE)
  r <- if (stats::sd(p) > 0 && stats::sd(o) > 0) stats::cor(o, p) else NA_real_
  st <- c(n = length(o), bias = mean(p - o), pbias = 100 * sum(p - o) / sum(o),
          rmse = sqrt(mean((p - o)^2)), mae = mean(abs(p - o)), r2 = r^2,
          nse = 1 - sum((p - o)^2) / sum((o - mean(o))^2),
          kge = 1 - sqrt((r - 1)^2 + (stats::sd(p) / stats::sd(o) - 1)^2 +
                           (mean(p) / mean(o) - 1)^2))
  cov_tab <- NULL
  if (inherits(pred, "dm_predictive")) {
    if (!is.null(pred$draws)) {
      D <- pred$draws[ok, , drop = FALSE]
      crps <- vapply(seq_along(o), function(i) {
        x <- sort(D[i, ]); M <- length(x)
        mean(abs(x - o[i])) - sum((2 * seq_len(M) - M - 1) * x) / M^2
      }, numeric(1))
      cover <- vapply(level, function(a) {
        lo <- apply(D, 1, stats::quantile, (1 - a) / 2, names = FALSE)
        hi <- apply(D, 1, stats::quantile, (1 + a) / 2, names = FALSE)
        mean(o >= lo & o <= hi)
      }, numeric(1))
    } else {
      s <- pmax(pred$sd[ok], 1e-12); z <- (o - p) / s
      crps <- s * (z * (2 * stats::pnorm(z) - 1) + 2 * stats::dnorm(z) - 1 / sqrt(pi))
      cover <- vapply(level, function(a) mean(abs(z) <= stats::qnorm((1 + a) / 2)), numeric(1))
    }
    st <- c(st, crps = mean(crps))
    cov_tab <- data.frame(nominal = level, observed = cover)
  }
  structure(list(stats = st, coverage = cov_tab), class = "dm_validation")
}

#' @export
print.dm_validation <- function(x, ...) {
  s <- x$stats
  cat("Validation against", s[["n"]], "measurements\n")
  print(round(s[-1], 3))
  nse <- s[["nse"]]; pb <- abs(s[["pbias"]])
  cat("Rating (Moriasi et al. 2007 sediment guide): ",
      if (nse > 0.75 && pb < 15) "very good" else if (nse > 0.65 && pb < 30) "good" else
        if (nse > 0.5 && pb < 55) "satisfactory" else "unsatisfactory", "\n", sep = "")
  if (!is.null(x$coverage)) {
    cat("Prediction-interval coverage (nominal -> observed):\n")
    for (i in seq_len(nrow(x$coverage))) {
      cat(sprintf("  %3.0f%% -> %3.0f%%\n", 100 * x$coverage$nominal[i], 100 * x$coverage$observed[i]))
    }
  }
  invisible(x)
}
