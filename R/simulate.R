# ---------------------------------------------------------------------------
# Simulated landscapes with known truth (for examples, teaching, testing)
# ---------------------------------------------------------------------------

.with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = globalenv()), add = TRUE)
  }
  set.seed(seed)
  expr
}

#' Simulate a spatial field with known truth and noisy observations
#'
#' Generates a Gaussian random field on a regular grid over the unit square,
#' plus a random sample of noisy observations. Because the true field is
#' known everywhere, the quality of maps and decisions can be measured
#' exactly. Used in examples, tests and the vignette.
#'
#' @param n_side Number of grid cells along each side (the grid has
#'   `n_side^2` cells).
#' @param n_obs Number of observed locations (sampled from grid cells).
#' @param model A `dm_covmodel`; its nugget is used as observation noise.
#' @param mean Constant mean of the field, or a function of `x` and `y`.
#' @param seed Optional random seed. The global random state is restored
#'   afterwards.
#' @return A list of class `dm_sim` with `grid` (columns `x`, `y`, `z` =
#'   true value), `obs` (columns `x`, `y`, `z` = observed value, `cell`) and
#'   `model`.
#' @examples
#' sim <- dm_simulate(n_side = 20, n_obs = 50, seed = 1)
#' sim
#' @export
dm_simulate <- function(n_side = 30, n_obs = 60,
                        model = dm_covmodel("exponential", psill = 1, range = 0.15,
                                            nugget = 0.05),
                        mean = 0, seed = NULL) {
  if (n_side < 3) stop("'n_side' must be at least 3.", call. = FALSE)
  N <- n_side^2
  if (n_obs < 1 || n_obs > N) stop("'n_obs' must be between 1 and n_side^2.", call. = FALSE)
  .with_seed(seed, {
    g <- (seq_len(n_side) - 0.5) / n_side
    grid <- expand.grid(x = g, y = g)
    C <- dm_cov(grid, model = model)
    L <- chol(C + diag(1e-8 * model$psill, N))
    mu <- if (is.function(mean)) mean(grid$x, grid$y) else rep(mean, N)
    grid$z <- as.numeric(mu + crossprod(L, stats::rnorm(N)))
    cells <- sort(sample.int(N, n_obs))
    obs <- data.frame(x = grid$x[cells], y = grid$y[cells],
                      z = grid$z[cells] + stats::rnorm(n_obs, 0, sqrt(model$nugget)),
                      cell = cells)
    structure(list(grid = grid, obs = obs, model = model), class = "dm_sim")
  })
}

#' @export
print.dm_sim <- function(x, ...) {
  cat("<dm_sim> ", nrow(x$grid), " grid cells, ", nrow(x$obs), " observations\n", sep = "")
  cat("  true field range: ", format(signif(min(x$grid$z), 3)), " to ",
      format(signif(max(x$grid$z), 3)), "\n", sep = "")
  print(x$model)
  invisible(x)
}
