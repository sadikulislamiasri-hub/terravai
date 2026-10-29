# ---------------------------------------------------------------------------
# AI and machine-learning predictive distributions
# ---------------------------------------------------------------------------

# Correlation matrix for a spatial Gaussian copula
.copula_cor <- function(copula, coords, n) {
  if (is.null(copula)) return(NULL)
  if (inherits(copula, "dm_covmodel")) {
    if (is.null(coords)) stop("A covariance-model copula needs 'coords'.", call. = FALSE)
    R <- dm_cov(coords, model = copula) / (copula$psill + copula$nugget)
    diag(R) <- 1
    return(R)
  }
  R <- as.matrix(copula)
  if (!all(dim(R) == n) || any(abs(diag(R) - 1) > 1e-8)) {
    stop("'copula' must be a dm_covmodel or an n x n correlation matrix.", call. = FALSE)
  }
  R
}

# Uniform draws (n x M), independent across draws, spatially dependent
# through an optional Gaussian copula
.copula_uniforms <- function(n, M, R = NULL) {
  Z <- matrix(stats::rnorm(n * M), n, M)
  if (!is.null(R)) {
    L <- chol(R + diag(1e-8, n))
    Z <- crossprod(L, Z)
  }
  stats::pnorm(Z)
}

#' Predictive distribution from predictive quantiles (with a spatial copula)
#'
#' Converts predictive quantiles, for example from a quantile regression
#' forest or a quantile neural network, into Monte Carlo draws. Each
#' location's quantile function is interpolated linearly between the supplied
#' quantiles (and held constant beyond them, so supply quantiles that cover the
#' tails you care about).
#'
#' Machine-learning models give a predictive distribution for each location
#' separately, but say nothing about how errors at nearby locations move
#' together. That joint behaviour matters for value of information and for
#' aggregates. `copula` supplies it through a spatial Gaussian copula: draws
#' keep each location's machine-learning distribution exactly, while their
#' ranks are correlated in space. Without a copula, locations are treated as
#' independent.
#'
#' @param q Numeric matrix of quantiles, one row per location and one column
#'   per probability in `probs`.
#' @param probs Strictly increasing probabilities in (0, 1) matching the
#'   columns of `q`.
#' @param coords Optional coordinates (required for a covariance-model copula).
#' @param n_draws Number of draws per location.
#' @param copula `NULL` (independent locations), a `dm_covmodel` whose
#'   correlation function links locations, or an `n x n` correlation matrix.
#' @param seed Optional random seed; the global random state is restored.
#' @return An object of class `dm_predictive` (type `"draws"`).
#' @examples
#' q <- rbind(c(4.2, 4.9, 5.6), c(5.5, 6.0, 6.6), c(5.0, 5.5, 6.1))
#' xy <- cbind(c(0, 0.1, 0.9), 0)
#' p <- dm_quantiles(q, probs = c(0.05, 0.5, 0.95), coords = xy,
#'                   copula = dm_covmodel("exponential", 1, 0.3), seed = 1)
#' round(cor(t(p$draws)), 2)
#' @export
dm_quantiles <- function(q, probs, coords = NULL, n_draws = 200, copula = NULL, seed = NULL) {
  q <- as.matrix(q)
  if (ncol(q) != length(probs)) stop("ncol(q) must equal length(probs).", call. = FALSE)
  if (is.unsorted(probs, strictly = TRUE) || any(probs <= 0 | probs >= 1)) {
    stop("'probs' must be strictly increasing within (0, 1).", call. = FALSE)
  }
  if (n_draws < 2) stop("'n_draws' must be at least 2.", call. = FALSE)
  n <- nrow(q)
  R <- .copula_cor(copula, coords, n)
  U <- .with_seed(seed, .copula_uniforms(n, n_draws, R))
  dr <- matrix(0, n, n_draws)
  for (i in seq_len(n)) {
    dr[i, ] <- stats::approx(probs, cummax(q[i, ]), xout = U[i, ], rule = 2,
                             ties = "ordered")$y
  }
  dm_draws(dr, coords = coords)
}

#' Conformal predictive distribution for any point predictor
#'
#' Wraps any AI or machine-learning point predictor (gradient boosting, a
#' neural network, a random forest) into a predictive distribution using
#' split-conformal calibration: the prediction errors on a held-out
#' calibration set define the spread of the distribution at new locations.
#' Errors can be scaled by a local difficulty measure (for example the spread
#' of an ensemble), which gives wider distributions where the model is less
#' sure. A spatial copula links errors at nearby locations.
#'
#' For spatial data, choose calibration locations that are spatially
#' separated from the training data (for example spatial blocks), otherwise
#' the distribution will be too narrow.
#'
#' @param pred Point predictions at the new locations.
#' @param calib_obs,calib_pred Observed values and the model's predictions at
#'   calibration locations not used for training.
#' @param scale,calib_scale Optional positive local scale (difficulty) at the
#'   new and calibration locations. Errors are standardised by it.
#' @inheritParams dm_quantiles
#' @return An object of class `dm_predictive` (type `"draws"`).
#' @references Vovk, V., Gammerman, A. and Shafer, G. (2005). *Algorithmic
#'   Learning in a Random World*. Springer, New York.
#' @examples
#' set.seed(1)
#' calib_obs <- rnorm(100); calib_pred <- calib_obs + rnorm(100, 0, 0.3)
#' p <- dm_conformal(pred = c(0.2, 1.1), calib_obs, calib_pred, seed = 2)
#' p
#' @export
dm_conformal <- function(pred, calib_obs, calib_pred, scale = NULL, calib_scale = NULL,
                         coords = NULL, n_draws = 200, copula = NULL, seed = NULL) {
  pred <- as.numeric(pred)
  if (length(calib_obs) != length(calib_pred) || length(calib_obs) < 10) {
    stop("Need at least 10 calibration pairs of equal length.", call. = FALSE)
  }
  if (xor(is.null(scale), is.null(calib_scale))) {
    stop("Supply both 'scale' and 'calib_scale', or neither.", call. = FALSE)
  }
  r <- as.numeric(calib_obs) - as.numeric(calib_pred)
  if (!is.null(scale)) {
    if (any(!(scale > 0)) || any(!(calib_scale > 0))) stop("Scales must be positive.", call. = FALSE)
    r <- r / calib_scale
  } else scale <- rep(1, length(pred))
  if (length(scale) == 1L) scale <- rep(scale, length(pred))
  n <- length(pred)
  R <- .copula_cor(copula, coords, n)
  U <- .with_seed(seed, .copula_uniforms(n, n_draws, R))
  E <- matrix(stats::quantile(r, U, type = 8, names = FALSE), n, n_draws)
  dm_draws(pred + scale * E, coords = coords)
}

#' Quantile regression forest with a spatial copula
#'
#' Fits a quantile regression forest (Meinshausen, 2006) with the 'ranger'
#' package and returns a spatially coherent predictive distribution at new
#' locations. Spatial dependence of the forest's errors is estimated by
#' fitting a covariance model to its out-of-bag residuals, and is imposed on
#' the forest's predictive quantiles through a spatial Gaussian copula
#' ([dm_quantiles()]). The result can be passed directly to [dm_decide()],
#' [dm_voi()] (use `method = "ensemble"`) and [dm_evaluate()].
#'
#' @param formula Model formula, e.g. `soc ~ ndvi + elevation`.
#' @param data Training data frame, including coordinate columns.
#' @param newdata Data frame of prediction locations with the same predictors
#'   and coordinate columns.
#' @param coords Names of the two coordinate columns.
#' @param probs Quantile levels predicted by the forest.
#' @param copula `TRUE` to estimate the spatial copula from out-of-bag
#'   residuals, `FALSE` for independent locations, or a `dm_covmodel`.
#' @param n_draws Number of draws per location.
#' @param num.trees Number of trees.
#' @param seed Random seed for the forest and the draws.
#' @param ... Further arguments for `ranger::ranger()`.
#' @return A `dm_predictive` object with attributes `forest` (the fitted
#'   'ranger' model) and `copula` (the covariance model used, if any).
#' @references Meinshausen, N. (2006). Quantile regression forests. *Journal
#'   of Machine Learning Research*, 7, 983-999.
#' @examples
#' if (requireNamespace("ranger", quietly = TRUE)) {
#'   sim <- dm_simulate(n_side = 12, n_obs = 80, seed = 3)
#'   sim$grid$cov1 <- sim$grid$z + rnorm(nrow(sim$grid), 0, 0.5)
#'   train <- cbind(sim$obs, cov1 = sim$grid$cov1[sim$obs$cell])
#'   p <- dm_fit_qrf(z ~ cov1, train, sim$grid, num.trees = 100, seed = 1)
#'   p
#' }
#' @export
dm_fit_qrf <- function(formula, data, newdata, coords = c("x", "y"),
                       probs = c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99),
                       copula = TRUE, n_draws = 200, num.trees = 500, seed = NULL, ...) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    stop("Package 'ranger' is needed for dm_fit_qrf().", call. = FALSE)
  }
  if (!all(coords %in% names(data)) || !all(coords %in% names(newdata))) {
    stop("Coordinate columns not found in 'data' and 'newdata'.", call. = FALSE)
  }
  rf <- ranger::ranger(formula, data = data, quantreg = TRUE, num.trees = num.trees,
                       keep.inbag = TRUE, seed = seed, ...)
  q <- stats::predict(rf, data = newdata, type = "quantiles", quantiles = probs)$predictions
  yname <- all.vars(formula)[1]
  cm <- NULL
  if (isTRUE(copula)) {
    res <- data[[yname]] - rf$predictions
    cm <- tryCatch(dm_fit_cov(data[, coords], res), error = function(e) NULL)
    if (is.null(cm)) warning("Could not fit a copula to out-of-bag residuals; using independence.",
                             call. = FALSE)
  } else if (inherits(copula, "dm_covmodel")) cm <- copula
  nc <- nrow(newdata)
  if (!is.null(cm) && nc > 5000) {
    warning("More than 5000 prediction locations: copula skipped (independent locations).",
            call. = FALSE)
    cm <- NULL
  }
  out <- dm_quantiles(q, probs, coords = newdata[, coords], n_draws = n_draws, copula = cm,
                      seed = seed)
  attr(out, "forest") <- rf
  attr(out, "copula") <- cm
  out
}

# ---------------------------------------------------------------------------
# Plain-language summaries
# ---------------------------------------------------------------------------

.where <- function(x, y, xr, yr) {
  fx <- (mean(x) - xr[1]) / max(diff(xr), 1e-12)
  fy <- (mean(y) - yr[1]) / max(diff(yr), 1e-12)
  ns <- if (fy > 2 / 3) "north" else if (fy < 1 / 3) "south" else ""
  ew <- if (fx > 2 / 3) "east" else if (fx < 1 / 3) "west" else ""
  w <- paste0(ns, ew)
  if (w == "") "central part" else paste(w, "part")
}

#' Plain-language summary of a decision map
#'
#' Writes a short, plain-language summary of a decision map for
#' non-specialists such as farmers, extension staff or district officers:
#' where action is advised, how confident the advice is, where it is
#' uncertain, how robust it is to disagreement about costs and, optionally,
#' where new samples would help most. The text is generated deterministically
#' from the numbers, so it can be checked and reproduced.
#'
#' @param x A `dm_decision` object.
#' @param voi Optional `dm_voi` or `dm_design` object.
#' @param area Name of the mapped area, used in the text.
#' @param uncertain Probability of a wrong decision above which a location is
#'   called uncertain.
#' @return A character string of class `dm_narrative`.
#' @examples
#' sim <- dm_simulate(n_side = 15, n_obs = 30, seed = 2)
#' fit <- dm_krige(sim$obs[, c("x", "y")], sim$obs$z, sim$grid[, c("x", "y")],
#'                 model = sim$model)
#' dec <- dm_decide(fit, dm_loss(0.5, "above", ratio = 3, action = "apply lime"))
#' dm_narrate(dec, area = "the block")
#' @export
dm_narrate <- function(x, voi = NULL, area = "the area", uncertain = 0.2) {
  if (!inherits(x, "dm_decision")) stop("'x' must be a dm_decision object.", call. = FALSE)
  loss <- attr(x, "loss")
  pct <- function(v) paste0(round(100 * v), "%")
  act <- x$act
  unc <- x$p_wrong > uncertain
  robust <- x$breakeven < loss$ratio / 4 | x$breakeven > loss$ratio * 4
  s <- character(0)
  s <- c(s, sprintf("Across %s, the recommended action \"%s\" applies to %s of the %d mapped locations.",
                    area, loss$action, pct(mean(act)), nrow(x)))
  s <- c(s, sprintf("Advice is %s: the chance that the recommendation is wrong is below %s at %s of locations.",
                    if (mean(!unc) > 0.8) "mostly firm" else if (mean(!unc) > 0.5) "partly firm" else "often uncertain",
                    pct(uncertain), pct(mean(!unc))))
  has_xy <- all(c("x", "y") %in% names(x))
  if (any(unc)) {
    where <- if (has_xy) paste0(", concentrated in the ",
                                .where(x$x[unc], x$y[unc], range(x$x), range(x$y)), " of ", area)
    else ""
    s <- c(s, sprintf("Uncertain advice covers %s of locations%s, where values lie close to the threshold of %s.",
                      pct(mean(unc)), where, format(loss$threshold)))
  }
  s <- c(s, sprintf("The advice at %s of locations would not change even if a missed action were judged four times more or less costly than the current ratio of %s.",
                    pct(mean(robust)), format(loss$ratio)))
  if (!is.null(voi)) {
    if (inherits(voi, "dm_design") && nrow(voi)) {
      rl <- attr(voi, "remaining_loss")
      s <- c(s, sprintf("The proposed %d new samples are expected to reduce decision losses by %s.",
                        nrow(voi), pct(1 - rl[length(rl)] / rl[1])))
    } else if (inherits(voi, "dm_voi") && all(c("x", "y") %in% names(voi))) {
      top <- voi$voi >= stats::quantile(voi$voi, 0.9)
      s <- c(s, sprintf("New samples would be most valuable in the %s of %s.",
                        .where(voi$x[top], voi$y[top], range(x$x), range(x$y)), area))
    }
  }
  structure(paste(s, collapse = " "), class = "dm_narrative")
}

#' @export
print.dm_narrative <- function(x, ...) {
  cat(strwrap(unclass(x), width = min(getOption("width"), 80)), sep = "\n")
  invisible(x)
}
