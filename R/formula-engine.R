# ---------------------------------------------------------------------------
# Shared engine: apply a formula to numbers, rasters or uncertain inputs
# ---------------------------------------------------------------------------

#' Uncertain input for formula-based models
#'
#' Describes an input that is not known exactly, such as a curve number,
#' an erodibility value or a rainfall total, so that its uncertainty can be
#' carried through the runoff, soil-loss and nutrient-loss formulas by
#' Monte Carlo simulation. The result of any formula with an uncertain input
#' is a `dm_predictive` object that the decision, value-of-information and
#' design functions accept directly.
#'
#' @param mean Mean value(s): a number, a vector (one per location) or a
#'   single-layer `SpatRaster`.
#' @param sd Standard deviation(s), or use `cv`.
#' @param cv Coefficient of variation (sd / mean), used when `sd` is `NULL`.
#' @param dist `"normal"` (truncated at `lower` and `upper`) or
#'   `"lognormal"` (always positive; `mean` and `sd` on the original scale).
#' @param lower,upper Physical bounds; normal draws are truncated to them.
#' @param coords Optional coordinates (taken from the raster when `mean` is a
#'   `SpatRaster`).
#' @return An object of class `dm_uncertain`.
#' @examples
#' cn <- dm_uncertain(78, sd = 4, lower = 30, upper = 100)
#' dm_runoff_cn(P = 80, CN = cn, n_draws = 500, seed = 1)
#' @export
dm_uncertain <- function(mean, sd = NULL, cv = NULL, dist = c("normal", "lognormal"),
                         lower = -Inf, upper = Inf, coords = NULL) {
  dist <- match.arg(dist)
  template <- cells <- NULL
  if (inherits(mean, "SpatRaster")) {
    r <- .from_raster(mean)
    mean <- r$values; coords <- r$coords; template <- r$template; cells <- r$cells
    if (inherits(sd, "SpatRaster")) sd <- .from_raster(sd)$values
  }
  if (is.null(sd)) {
    if (is.null(cv)) stop("Supply 'sd' or 'cv'.", call. = FALSE)
    sd <- abs(mean) * cv
  }
  n <- max(length(mean), length(sd))
  mean <- rep_len(mean, n); sd <- rep_len(sd, n)
  if (any(!is.finite(mean)) || any(!is.finite(sd)) || any(sd < 0)) {
    stop("'mean' and 'sd' must be finite with sd >= 0.", call. = FALSE)
  }
  if (dist == "lognormal" && any(mean <= 0)) stop("Lognormal inputs need mean > 0.", call. = FALSE)
  structure(list(mean = mean, sd = sd, dist = dist, lower = lower, upper = upper,
                 coords = if (!is.null(coords)) .as_coords(coords, n) else NULL,
                 template = template, cells = cells), class = "dm_uncertain")
}

#' @export
print.dm_uncertain <- function(x, ...) {
  cat("<dm_uncertain> ", length(x$mean), " value(s), ", x$dist, ", mean ",
      format(signif(min(x$mean), 4)), " to ", format(signif(max(x$mean), 4)),
      ", sd ", format(signif(min(x$sd), 3)), " to ", format(signif(max(x$sd), 3)), "\n", sep = "")
  invisible(x)
}

.draw_uncertain <- function(u, M) {
  n <- length(u$mean)
  if (u$dist == "lognormal") {
    s2 <- log(1 + (u$sd / u$mean)^2)
    mu <- log(u$mean) - s2 / 2
    return(matrix(exp(stats::rnorm(n * M, mu, sqrt(s2))), n, M))
  }
  # truncated normal by inversion
  a <- stats::pnorm((u$lower - u$mean) / pmax(u$sd, 1e-300))
  b <- stats::pnorm((u$upper - u$mean) / pmax(u$sd, 1e-300))
  U <- matrix(stats::runif(n * M), n, M)
  P <- a + U * (b - a)
  P <- pmin(pmax(P, 1e-12), 1 - 1e-12)
  out <- u$mean + u$sd * stats::qnorm(P)
  out[u$sd == 0, ] <- u$mean[u$sd == 0]
  pmin(pmax(out, u$lower), u$upper)
}

.is_uncertain <- function(a) inherits(a, c("dm_uncertain", "dm_predictive"))

# Evaluate `fun` (vectorised in its arguments) on any mix of numbers,
# SpatRasters, dm_uncertain and dm_predictive inputs.
.apply_formula <- function(fun, args, n_draws = 500, seed = NULL) {
  unc <- vapply(args, .is_uncertain, logical(1))
  rst <- vapply(args, function(a) inherits(a, "SpatRaster"), logical(1))
  if (any(unc)) {
    lens <- vapply(args[unc], function(a) length(a$mean), integer(1))
    n <- max(lens)
    first <- args[unc][[which.max(lens)]]
    M <- n_draws
    res <- .with_seed(seed, {
      mats <- lapply(names(args), function(nm) {
        a <- args[[nm]]
        if (inherits(a, "dm_uncertain")) {
          D <- .draw_uncertain(a, M)
        } else if (inherits(a, "dm_predictive")) {
          D <- if (!is.null(a$draws)) {
            a$draws[, rep_len(seq_len(ncol(a$draws)), M), drop = FALSE]
          } else {
            matrix(stats::rnorm(length(a$mean) * M, a$mean, a$sd), length(a$mean), M)
          }
        } else if (inherits(a, "SpatRaster")) {
          stop("Do not mix SpatRaster and uncertain inputs; wrap the raster with ",
               "dm_uncertain() or extract its values.", call. = FALSE)
        } else {
          D <- matrix(rep_len(a, n), n, M)
        }
        if (nrow(D) == 1L && n > 1L) D <- D[rep(1L, n), , drop = FALSE]
        if (nrow(D) != n) stop("Uncertain inputs have different lengths.", call. = FALSE)
        D
      })
      names(mats) <- names(args)
      do.call(fun, mats)
    })
    out <- dm_draws(matrix(res, n, M), coords = first$coords)
    out$template <- first$template; out$cells <- first$cells
    return(out)
  }
  if (any(rst)) {
    if (!requireNamespace("terra", quietly = TRUE)) stop("Package 'terra' is needed.", call. = FALSE)
    tmpl <- args[rst][[1]][[1]]
    vals <- lapply(args, function(a) if (inherits(a, "SpatRaster")) terra::values(a[[1]], mat = FALSE) else a)
    v <- do.call(fun, vals)
    out <- terra::rast(tmpl)
    terra::values(out) <- as.numeric(v)
    return(out)
  }
  do.call(fun, args)
}

.check_range <- function(x, nm, lower = -Inf, upper = Inf) {
  v <- if (.is_uncertain(x)) x$mean else if (inherits(x, "SpatRaster")) NULL else x
  if (!is.null(v)) {
    v <- v[is.finite(v)]
    if (any(v < lower) || any(v > upper)) {
      stop("'", nm, "' must lie between ", lower, " and ", upper, ".", call. = FALSE)
    }
  }
  invisible(TRUE)
}
