# ---------------------------------------------------------------------------
# Transferring nutrient-loss relations to a new site with local calibration
# ---------------------------------------------------------------------------

#' Nutrient loss at a new site, calibrated with a few local plots
#'
#' Nutrient loss per tonne of eroded soil differs strongly between sites
#' (soil nutrient content, enrichment, fertiliser history), so a relation
#' fitted elsewhere transfers poorly. This function fits the multi-site
#' model
#' \deqn{\log Y_{ij} = \mu_j + b \log A_{ij} + \varepsilon_{ij}, \quad
#' \mu_j \sim N(m, \tau^2), \quad \varepsilon_{ij} \sim N(0, \sigma^2),}
#' with \eqn{Y} nutrient loss, \eqn{A} soil loss and \eqn{j} the site: a slope
#' \eqn{b} shared by all sites (estimated within sites) and a site effect
#' \eqn{\mu_j}. For a new site the site effect is unknown; it is updated from
#' any local measurements by a conjugate normal step. The predictive draws
#' share one site effect, so their errors are correlated across the site,
#' which is what [dm_voi()] and [dm_design()] need to value local plots.
#'
#' @param soil_loss,nutrient_loss,site Training data from several sites:
#'   soil loss (t/ha), nutrient loss (kg/ha) and a site identifier. Records
#'   with soil loss at or below `min_soil` or non-positive nutrient loss are
#'   not used.
#' @param new_soil_loss Soil loss at the new site's prediction locations (a
#'   number, vector or `dm_predictive`; its mean is used).
#' @param local_soil_loss,local_nutrient Optional local calibration plots at
#'   the new site.
#' @param min_soil Smallest soil loss used in fitting (t/ha).
#' @param n_draws,seed Monte Carlo settings.
#' @return A `dm_predictive` of nutrient loss (kg/ha) with draws, and the
#'   attributes `params` (`b`, `m`, `tau`, `sigma`, between-site share of
#'   variance) and `site_effect` (posterior mean and variance of the new site
#'   effect). Use `log` of the draws, or [dm_gaussian()] on the log scale with
#'   covariance \eqn{\tau_{post}^2 J + \sigma^2 I}, for closed-form value of
#'   information.
#' @examples
#' set.seed(1)
#' site <- rep(c("A", "B", "C", "D"), each = 30)
#' A <- rlnorm(120, 1.5, 0.8)
#' mu <- c(A = 0, B = 0.8, C = -0.6, D = 0.3)[site]
#' N <- exp(mu + 0.5 * log(A) + rnorm(120, 0, 0.5))
#' p0 <- dm_nutrient_transfer(A, N, site, new_soil_loss = c(2, 5, 10))
#' p3 <- dm_nutrient_transfer(A, N, site, new_soil_loss = c(2, 5, 10),
#'                            local_soil_loss = c(3, 6, 9),
#'                            local_nutrient = c(4, 6, 8))
#' attr(p0, "params")
#' rbind(no_local = p0$sd, three_local = p3$sd)
#' @export
dm_nutrient_transfer <- function(soil_loss, nutrient_loss, site, new_soil_loss,
                                 local_soil_loss = NULL, local_nutrient = NULL,
                                 min_soil = 0.5, n_draws = 1000, seed = NULL) {
  if (!(length(soil_loss) == length(nutrient_loss) && length(site) == length(soil_loss))) {
    stop("'soil_loss', 'nutrient_loss' and 'site' must have the same length.", call. = FALSE)
  }
  ok <- is.finite(soil_loss) & is.finite(nutrient_loss) & soil_loss > min_soil &
    nutrient_loss > 0 & !is.na(site)
  site <- factor(site[ok])
  if (nlevels(site) < 3) stop("At least 3 sites are needed to estimate between-site variation.",
                              call. = FALSE)
  y <- log(nutrient_loss[ok]); la <- log(soil_loss[ok])
  fe <- stats::lm(y ~ la + site)
  b <- unname(stats::coef(fe)["la"])
  s2 <- sum(stats::residuals(fe)^2) / fe$df.residual
  mc <- tapply(y - b * la, site, mean); nc <- tapply(y, site, length)
  m <- mean(mc); tau2 <- max(stats::var(mc) - mean(s2 / nc), 1e-4)
  post <- c(m, tau2)
  if (!is.null(local_nutrient)) {
    if (length(local_soil_loss) != length(local_nutrient)) {
      stop("'local_soil_loss' and 'local_nutrient' must have the same length.", call. = FALSE)
    }
    keep <- local_soil_loss > 0 & local_nutrient > 0
    rl <- log(local_nutrient[keep]) - b * log(local_soil_loss[keep])
    if (length(rl)) {
      v <- 1 / (1 / tau2 + length(rl) / s2)
      post <- c(v * (m / tau2 + sum(rl) / s2), v)
    }
  }
  A0 <- if (inherits(new_soil_loss, "dm_predictive")) new_soil_loss$mean else new_soil_loss
  A0 <- pmax(A0, 1e-3)
  n <- length(A0)
  D <- .with_seed(seed, {
    mu <- stats::rnorm(n_draws, post[1], sqrt(post[2]))
    eps <- matrix(stats::rnorm(n * n_draws, 0, sqrt(s2)), n, n_draws)
    exp(b * log(A0) + eps + matrix(mu, n, n_draws, byrow = TRUE))
  })
  out <- dm_draws(D)
  attr(out, "params") <- c(b = b, m = m, tau = sqrt(tau2), sigma = sqrt(s2),
                           between_share = tau2 / (tau2 + s2))
  attr(out, "site_effect") <- c(mean = post[1], var = post[2])
  attr(out, "log_mean") <- b * log(A0) + post[1]
  out
}
