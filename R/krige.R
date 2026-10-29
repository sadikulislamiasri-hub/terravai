# ---------------------------------------------------------------------------
# A small, dependency-free kriging engine with full predictive covariance
# ---------------------------------------------------------------------------

#' Covariance model
#'
#' Defines a stationary, isotropic covariance model for the latent spatial
#' process. The nugget is treated as measurement error (and micro-scale
#' variation): it is added to observations but not to the latent field that
#' decisions concern.
#'
#' @param model One of `"exponential"`, `"spherical"`, `"gaussian"` or
#'   `"matern"`.
#' @param psill Partial sill (variance of the latent process).
#' @param range Range parameter, in coordinate units. For the exponential
#'   model the practical range is about `3 * range`.
#' @param nugget Nugget variance.
#' @param kappa Smoothness for the Matern model.
#' @return An object of class `dm_covmodel`.
#' @examples
#' dm_covmodel("exponential", psill = 1, range = 0.2, nugget = 0.1)
#' @export
dm_covmodel <- function(model = c("exponential", "spherical", "gaussian", "matern"),
                        psill = 1, range = 1, nugget = 0, kappa = 0.5) {
  model <- match.arg(model)
  if (!(psill > 0) || !(range > 0) || nugget < 0 || !(kappa > 0)) {
    stop("Require psill > 0, range > 0, nugget >= 0 and kappa > 0.", call. = FALSE)
  }
  structure(list(model = model, psill = psill, range = range, nugget = nugget,
                 kappa = kappa), class = "dm_covmodel")
}

#' @export
print.dm_covmodel <- function(x, ...) {
  cat("<dm_covmodel> ", x$model, ": psill ", format(signif(x$psill, 4)), ", range ",
      format(signif(x$range, 4)), ", nugget ", format(signif(x$nugget, 4)),
      if (x$model == "matern") paste0(", kappa ", format(x$kappa)) else "", "\n", sep = "")
  if (!is.null(x$loglik)) {
    cat("  fitted by ", x$method, ": log-likelihood ", format(round(x$loglik, 3)),
        "; trend coefficients ", paste(format(signif(x$beta, 4)), collapse = ", "), "\n", sep = "")
  }
  invisible(x)
}

.cov_d <- function(d, cm) {
  h <- d / cm$range
  switch(cm$model,
    exponential = cm$psill * exp(-h),
    gaussian = cm$psill * exp(-h^2),
    spherical = cm$psill * ifelse(h < 1, 1 - 1.5 * h + 0.5 * h^3, 0),
    matern = {
      k <- cm$kappa
      out <- cm$psill * (2^(1 - k) / gamma(k)) * h^k * besselK(pmax(h, 1e-10), k)
      out[h < 1e-10] <- cm$psill
      out
    })
}

.crossdist <- function(a, b) {
  a <- as.matrix(a); b <- as.matrix(b)
  sqrt(pmax(outer(a[, 1], b[, 1], "-")^2 + outer(a[, 2], b[, 2], "-")^2, 0))
}

#' Covariance matrix between two sets of locations
#'
#' @param coords1,coords2 Two-column matrices or data frames of coordinates.
#'   If `coords2` is `NULL`, the covariance of `coords1` with itself.
#' @param model A `dm_covmodel`.
#' @param nugget Add the nugget on the diagonal (for observations).
#' @return A covariance matrix of the latent process.
#' @examples
#' xy <- cbind(c(0, 0.1, 0.5), c(0, 0, 0))
#' dm_cov(xy, model = dm_covmodel("exponential", 1, 0.2))
#' @export
dm_cov <- function(coords1, coords2 = NULL, model, nugget = FALSE) {
  if (!inherits(model, "dm_covmodel")) stop("'model' must be a dm_covmodel.", call. = FALSE)
  same <- is.null(coords2)
  C <- .cov_d(.crossdist(coords1, if (same) coords1 else coords2), model)
  if (same && nugget) diag(C) <- diag(C) + model$nugget
  C
}

.safe_chol <- function(K) {
  L <- tryCatch(chol(K), error = function(e) NULL)
  if (!is.null(L)) return(L)
  jit <- 1e-10 * mean(diag(K))
  for (i in 1:8) {
    L <- tryCatch(chol(K + diag(jit, nrow(K))), error = function(e) NULL)
    if (!is.null(L)) {
      warning("Covariance matrix was not positive definite (duplicate locations or a zero ",
              "nugget?); added ", format(signif(jit, 2)), " to the diagonal.", call. = FALSE)
      return(L)
    }
    jit <- jit * 100
  }
  stop("Covariance matrix is not positive definite.", call. = FALSE)
}

.design_X <- function(X, n) {
  if (is.null(X)) return(matrix(1, n, 1, dimnames = list(NULL, "(Intercept)")))
  X <- as.matrix(X)
  if (nrow(X) != n) stop("Covariate matrix has the wrong number of rows.", call. = FALSE)
  cbind("(Intercept)" = 1, X)
}

.negll <- function(par, D, y, X, model, kappa, reml) {
  cm <- dm_covmodel(model, psill = exp(par[1]), range = exp(par[2]),
                    nugget = exp(par[3]), kappa = kappa)
  K <- .cov_d(D, cm); diag(K) <- diag(K) + cm$nugget
  L <- tryCatch(chol(K), error = function(e) NULL)
  if (is.null(L)) return(1e10)
  Ki_X <- backsolve(L, forwardsolve(t(L), X))
  Ki_y <- backsolve(L, forwardsolve(t(L), y))
  XtKX <- crossprod(X, Ki_X)
  beta <- tryCatch(solve(XtKX, crossprod(X, Ki_y)), error = function(e) NULL)
  if (is.null(beta)) return(1e10)
  r <- y - X %*% beta
  Ki_r <- backsolve(L, forwardsolve(t(L), r))
  val <- sum(log(diag(L))) + 0.5 * sum(r * Ki_r)
  if (reml) val <- val + 0.5 * as.numeric(determinant(XtKX, logarithm = TRUE)$modulus)
  val
}

#' Fit a covariance model by (restricted) maximum likelihood
#'
#' Estimates the partial sill, range and nugget of a covariance model, with a
#' linear trend on optional covariates, by restricted maximum likelihood
#' (default) or maximum likelihood. Suitable for up to a few thousand
#' observations.
#'
#' @param coords Two-column matrix or data frame of observation coordinates.
#' @param y Observed values.
#' @param X Optional covariate matrix (an intercept is always added).
#' @param model Covariance model type (see [dm_covmodel()]).
#' @param kappa Matern smoothness (fixed).
#' @param reml Use restricted maximum likelihood.
#' @param start Optional starting values `c(psill, range, nugget)`.
#' @param n_starts Number of starting values tried (1 to 4); the fit with the
#'   highest likelihood is kept.
#' @return A fitted `dm_covmodel` with extra elements `beta` (trend
#'   coefficients), `loglik`, `method`, `convergence` (0 = converged) and
#'   `at_bound` (parameters at the search limits). A warning is given when
#'   the optimiser does not converge or a parameter ends at its limit.
#' @examples
#' sim <- dm_simulate(n_side = 10, n_obs = 60, seed = 5)
#' dm_fit_cov(sim$obs[, c("x", "y")], sim$obs$z)
#' @export
dm_fit_cov <- function(coords, y, X = NULL, model = "exponential", kappa = 0.5,
                       reml = TRUE, start = NULL, n_starts = 4) {
  y <- as.numeric(y)
  n <- length(y)
  if (n < 10) stop("At least 10 observations are needed.", call. = FALSE)
  coords <- as.matrix(coords)[, 1:2, drop = FALSE]
  if (nrow(coords) != n) stop("'coords' and 'y' differ in length.", call. = FALSE)
  Xd <- .design_X(X, n)
  D <- .crossdist(coords, coords)
  v <- stats::var(as.numeric(stats::lm.fit(Xd, y)$residuals))
  if (!(v > 0)) stop("The data have no variation.", call. = FALSE)
  start <- start %||% c(0.7 * v, max(D) / 6, 0.3 * v)
  model <- match.arg(model, c("exponential", "spherical", "gaussian", "matern"))
  lo <- log(c(v * 1e-4, max(D) * 1e-3, v * 1e-6))
  hi <- log(c(v * 20, max(D) * 5, v * 20))
  # multi-start optimisation: the likelihood surface of covariance models is
  # often flat or multimodal in the range, so several starts are compared
  starts <- list(start, c(0.5 * v, max(D) / 20, 0.5 * v), c(0.9 * v, max(D) / 3, 0.1 * v),
                 c(0.3 * v, max(D) / 10, 0.7 * v))[seq_len(max(1, min(n_starts, 4)))]
  fits <- lapply(starts, function(st) {
    tryCatch(stats::optim(pmin(pmax(log(st), lo), hi), .negll, D = D, y = y, X = Xd,
                          model = model, kappa = kappa, reml = reml, method = "L-BFGS-B",
                          lower = lo, upper = hi), error = function(e) NULL)
  })
  fits <- Filter(Negate(is.null), fits)
  if (!length(fits)) stop("Covariance fitting failed from every starting value.", call. = FALSE)
  opt <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  at_bound <- c(psill = FALSE, range = FALSE, nugget = FALSE)
  at_bound[] <- (opt$par - lo) < 1e-3 | (hi - opt$par) < 1e-3
  at_bound[["nugget"]] <- (hi[3] - opt$par[3]) < 1e-3     # a zero nugget is legitimate
  if (opt$convergence != 0) {
    warning("Covariance fitting did not converge (optim code ", opt$convergence, ").",
            call. = FALSE)
  }
  if (any(at_bound)) {
    warning("Fitted ", paste(names(at_bound)[at_bound], collapse = ", "),
            " at the search limit: spatial dependence is weakly identified ",
            "(e.g. range larger than the study area). Interpret with care or fix it via 'start'.",
            call. = FALSE)
  }
  cm <- dm_covmodel(model, psill = exp(opt$par[1]), range = exp(opt$par[2]),
                    nugget = exp(opt$par[3]), kappa = kappa)
  K <- dm_cov(coords, model = cm, nugget = TRUE)
  L <- chol(K)
  Ki_X <- backsolve(L, forwardsolve(t(L), Xd))
  Ki_y <- backsolve(L, forwardsolve(t(L), y))
  beta <- drop(solve(crossprod(Xd, Ki_X), crossprod(Xd, Ki_y)))
  names(beta) <- colnames(Xd) %||% paste0("b", seq_along(beta))
  cm$beta <- beta
  cm$loglik <- -opt$value - (n - if (reml) ncol(Xd) else 0) / 2 * log(2 * pi)
  cm$method <- if (reml) "REML" else "ML"
  cm$convergence <- opt$convergence
  cm$at_bound <- at_bound
  cm
}

#' Kriging with full predictive covariance
#'
#' Universal kriging (ordinary kriging when `X` is `NULL`) of the latent
#' process at new locations. Returns the full predictive covariance, including
#' the uncertainty of the estimated trend, which is what value-of-information
#' and aggregation calculations need.
#'
#' @param coords Observation coordinates.
#' @param y Observed values.
#' @param newcoords Prediction coordinates.
#' @param model A `dm_covmodel` (for example from [dm_fit_cov()]).
#' @param X,newX Optional covariates at observation and prediction locations.
#' @param full_cov Return the full predictive covariance matrix (needed for
#'   [dm_voi()]); if `FALSE`, only standard deviations.
#' @return A `dm_predictive` object of type `"gaussian"`.
#' @examples
#' sim <- dm_simulate(n_side = 12, n_obs = 30, seed = 5)
#' fit <- dm_krige(sim$obs[, c("x", "y")], sim$obs$z, sim$grid[, c("x", "y")],
#'                 model = sim$model)
#' fit
#' @export
dm_krige <- function(coords, y, newcoords, model, X = NULL, newX = NULL, full_cov = TRUE) {
  if (!inherits(model, "dm_covmodel")) stop("'model' must be a dm_covmodel.", call. = FALSE)
  y <- as.numeric(y); n <- length(y)
  coords <- as.matrix(coords)[, 1:2, drop = FALSE]; newcoords <- as.matrix(newcoords)[, 1:2, drop = FALSE]
  m <- nrow(newcoords)
  if (xor(is.null(X), is.null(newX))) stop("Supply both 'X' and 'newX', or neither.", call. = FALSE)
  Xd <- .design_X(X, n); X0 <- .design_X(newX, m)
  K <- dm_cov(coords, model = model, nugget = TRUE)
  L <- .safe_chol(K)
  k0 <- dm_cov(coords, newcoords, model)
  Ki_k0 <- backsolve(L, forwardsolve(t(L), k0))
  Ki_X <- backsolve(L, forwardsolve(t(L), Xd))
  Ki_y <- backsolve(L, forwardsolve(t(L), y))
  M <- solve(crossprod(Xd, Ki_X))
  beta <- M %*% crossprod(Xd, Ki_y)
  r <- y - Xd %*% beta
  mu <- as.numeric(X0 %*% beta + crossprod(k0, backsolve(L, forwardsolve(t(L), r))))
  Q <- X0 - crossprod(k0, Ki_X)
  if (full_cov) {
    S <- dm_cov(newcoords, model = model) - crossprod(k0, Ki_k0) + Q %*% M %*% t(Q)
    S <- (S + t(S)) / 2
    diag(S) <- pmax(diag(S), 0)
    out <- dm_gaussian(mu, cov = S, coords = newcoords)
  } else {
    v <- model$psill - colSums(k0 * Ki_k0) + rowSums((Q %*% M) * Q)
    out <- dm_gaussian(mu, sd = sqrt(pmax(v, 0)), coords = newcoords)
  }
  out
}
