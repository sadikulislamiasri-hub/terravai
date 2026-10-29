# ---------------------------------------------------------------------------
# Predictive distributions
# ---------------------------------------------------------------------------

.as_coords <- function(coords, n) {
  if (is.null(coords)) return(NULL)
  coords <- as.data.frame(coords)
  if (ncol(coords) < 2) stop("'coords' must have two columns (x, y).", call. = FALSE)
  coords <- coords[, 1:2]
  names(coords) <- c("x", "y")
  if (nrow(coords) != n) {
    stop("'coords' has ", nrow(coords), " rows but there are ", n,
         " predictions.", call. = FALSE)
  }
  rownames(coords) <- NULL
  coords
}

.new_predictive <- function(mean, sd, cov = NULL, draws = NULL, coords = NULL,
                            type, template = NULL, cells = NULL) {
  structure(list(mean = as.numeric(mean), sd = as.numeric(sd), cov = cov,
                 draws = draws, coords = coords, type = type,
                 template = template, cells = cells),
            class = "dm_predictive")
}

.from_raster <- function(r) {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is needed for SpatRaster input.", call. = FALSE)
  }
  tab <- terra::as.data.frame(r, xy = TRUE, cells = TRUE, na.rm = TRUE)
  list(values = tab[[4]], coords = tab[, c("x", "y")], cells = tab$cell,
       template = terra::rast(r[[1]]))
}

#' Gaussian predictive distribution
#'
#' Describes predictions by their means and standard deviations, optionally
#' with the full predictive covariance matrix. This is the natural output of
#' kriging, Gaussian-process and many Bayesian models. The covariance matrix is
#' needed for value-of-information calculations in which one sample informs
#' several locations ([dm_voi()]).
#'
#' @param mean Numeric vector of predictive means, or a single-layer
#'   `SpatRaster` (package `terra`).
#' @param sd Numeric vector of predictive standard deviations, or a
#'   `SpatRaster`. Ignored when `cov` is supplied.
#' @param cov Optional predictive covariance matrix.
#' @param coords Optional two-column matrix or data frame of coordinates.
#'   Taken from the raster when `mean` is a `SpatRaster`.
#' @return An object of class `dm_predictive`.
#' @seealso [dm_draws()], [dm_quantiles()], [dm_krige()]
#' @examples
#' p <- dm_gaussian(mean = c(4.8, 5.4, 6.1), sd = c(0.3, 0.4, 0.2))
#' p
#' @export
dm_gaussian <- function(mean, sd = NULL, cov = NULL, coords = NULL) {
  template <- cells <- NULL
  if (inherits(mean, "SpatRaster")) {
    rm <- .from_raster(mean)
    if (!inherits(sd, "SpatRaster")) stop("'sd' must also be a SpatRaster.", call. = FALSE)
    rs <- .from_raster(sd)
    if (!identical(rm$cells, rs$cells)) {
      stop("'mean' and 'sd' rasters must have identical non-missing cells.", call. = FALSE)
    }
    mean <- rm$values; sd <- rs$values; coords <- rm$coords
    template <- rm$template; cells <- rm$cells
  }
  n <- length(mean)
  if (!is.null(cov)) {
    cov <- as.matrix(cov)
    if (!all(dim(cov) == n)) stop("'cov' must be an n x n matrix.", call. = FALSE)
    sd <- sqrt(pmax(diag(cov), 0))
  }
  if (is.null(sd)) stop("Supply 'sd' or 'cov'.", call. = FALSE)
  if (length(sd) == 1L) sd <- rep(sd, n)
  if (length(sd) != n) stop("'sd' must have the same length as 'mean'.", call. = FALSE)
  if (any(!is.finite(mean)) || any(!is.finite(sd)) || any(sd < 0)) {
    stop("'mean' and 'sd' must be finite, with 'sd' >= 0.", call. = FALSE)
  }
  .new_predictive(mean, sd, cov = cov, coords = .as_coords(coords, n),
                  type = "gaussian", template = template, cells = cells)
}

#' Predictive distribution from simulation draws
#'
#' Describes predictions by Monte Carlo draws: conditional simulations from
#' geostatistics, posterior samples from Bayesian models, or ensemble members.
#' Non-Gaussian shapes (for example skewed contaminant concentrations) are
#' preserved exactly in decision maps; value-of-information calculations use a
#' Gaussian approximation based on the draws' means and covariance.
#'
#' @param draws Numeric matrix with one row per location and one column per
#'   draw.
#' @param coords Optional two-column matrix or data frame of coordinates.
#' @return An object of class `dm_predictive`.
#' @examples
#' set.seed(1)
#' d <- matrix(rlnorm(3 * 200, log(c(8, 15, 40)), 0.5), nrow = 3)
#' dm_draws(d)
#' @export
dm_draws <- function(draws, coords = NULL) {
  draws <- as.matrix(draws)
  if (!is.numeric(draws) || any(!is.finite(draws))) {
    stop("'draws' must be a finite numeric matrix.", call. = FALSE)
  }
  if (ncol(draws) < 2) stop("'draws' needs at least two columns (draws).", call. = FALSE)
  n <- nrow(draws)
  .new_predictive(rowMeans(draws), apply(draws, 1, sd), draws = draws,
                  coords = .as_coords(coords, n), type = "draws")
}

#' @export
print.dm_predictive <- function(x, ...) {
  cat("<dm_predictive> ", length(x$mean), " location(s), type: ", x$type, "\n", sep = "")
  cat("  mean range : ", format(signif(min(x$mean), 4)), " to ",
      format(signif(max(x$mean), 4)), "\n", sep = "")
  cat("  sd range   : ", format(signif(min(x$sd), 4)), " to ",
      format(signif(max(x$sd), 4)), "\n", sep = "")
  cat("  covariance : ", if (!is.null(x$cov)) "full matrix" else if (!is.null(x$draws))
    paste(ncol(x$draws), "draws") else "diagonal only", "\n", sep = "")
  cat("  coordinates: ", if (is.null(x$coords)) "none" else "yes", "\n", sep = "")
  invisible(x)
}

#' @export
length.dm_predictive <- function(x) length(x$mean)

# Covariance used for value of information -------------------------------------
.Sigma <- function(p, warn = TRUE) {
  if (!is.null(p$cov)) return(p$cov)
  if (!is.null(p$draws)) {
    if (nrow(p$draws) > 5000) {
      stop("Too many locations to estimate a covariance from draws; ",
           "use 'candidates' or aggregate first.", call. = FALSE)
    }
    return(stats::cov(t(p$draws)))
  }
  if (warn) {
    warning("No covariance available: assuming independent locations, so a sample ",
            "informs only its own location.", call. = FALSE)
  }
  diag(p$sd^2, nrow = length(p$sd))
}
