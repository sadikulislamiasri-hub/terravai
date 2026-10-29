# ---------------------------------------------------------------------------
# Value of information and sampling design
# ---------------------------------------------------------------------------

# Gauss-Hermite nodes and weights for a standard normal (Golub-Welsch)
.gh_normal <- function(k) {
  if (k < 2) stop("'n_quad' must be at least 2.", call. = FALSE)
  J <- matrix(0, k, k)
  off <- sqrt(seq_len(k - 1))
  J[cbind(1:(k - 1), 2:k)] <- off
  J[cbind(2:k, 1:(k - 1))] <- off
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = e$vectors[1, ]^2)
}

.min_loss <- function(m, s, loss) {
  el <- .exp_loss_gauss(m, s, loss)
  pmin(el$act, el$no)
}

.resolve_candidates <- function(candidates, n) {
  if (is.null(candidates)) return(seq_len(n))
  if (is.logical(candidates)) {
    if (length(candidates) != n) stop("Logical 'candidates' must have length n.", call. = FALSE)
    return(which(candidates))
  }
  candidates <- as.integer(candidates)
  if (any(candidates < 1 | candidates > n)) stop("'candidates' out of range.", call. = FALSE)
  unique(candidates)
}

# Decision quantities: Q = H Z. Returns mean, variance and cross-covariance
# with candidate observations.
.decision_moments <- function(mu, Sigma, target, cand) {
  if (is.null(target)) {
    list(m = mu, v = diag(Sigma), C = Sigma[, cand, drop = FALSE])
  } else {
    H <- as.matrix(target)
    if (ncol(H) != length(mu)) stop("'target' must have one column per location.", call. = FALSE)
    HS <- H %*% Sigma
    list(m = as.numeric(H %*% mu), v = rowSums(HS * H), C = HS[, cand, drop = FALSE])
  }
}

# Ensemble Kalman preposterior value of information ---------------------------
.min_loss_draws <- function(D, loss) {
  el <- .exp_loss_draws(D, loss)
  pmin(el$act, el$no)
}

.voi_ensemble <- function(D, loss, cand, noise, target, n_outer, localize = 2) {
  Q <- if (is.null(target)) D else as.matrix(target) %*% D
  M <- ncol(D)
  L0 <- .min_loss_draws(Q, loss)
  live <- L0 > 1e-9 * max(L0, 1e-300)
  Qc <- Q - rowMeans(Q)
  sdq <- sqrt(rowSums(Qc^2) / (M - 1))
  u <- stats::ppoints(n_outer)
  out <- numeric(length(cand))
  for (i in seq_along(cand)) {
    s <- cand[i]
    y <- D[s, ] + stats::rnorm(M, 0, sqrt(noise))   # perturbed observations
    yc <- y - mean(y)
    vy <- sum(yc^2)
    if (vy <= 0) next
    cq <- as.numeric(Qc %*% yc)
    K <- cq / vy                                        # Kalman gains
    # Localisation: ignore gains whose ensemble correlation is not
    # distinguishable from zero (|r| < localize / sqrt(M)).
    r <- cq / ((M - 1) * pmax(sdq, 1e-300) * sqrt(vy / (M - 1)))
    rows <- which(live & abs(r) > localize / sqrt(M) & abs(K) > 0)
    if (!length(rows)) next
    Qr <- Q[rows, , drop = FALSE]; Kr <- K[rows]
    yk <- stats::quantile(y, u, type = 8, names = FALSE)
    after <- numeric(length(rows))
    for (k in seq_along(yk)) {
      upd <- Qr + outer(Kr, yk[k] - y)
      after <- after + .min_loss_draws(upd, loss) / length(yk)
    }
    out[i] <- sum(L0[rows] - after)
  }
  list(voi = pmax(out, 0), evpi = sum(L0))
}

# Core: expected loss reduction for each candidate (columns of C)
.voi_core <- function(m, v, C, sv_cand, noise, loss, gh, chunk = 256) {
  L0 <- .min_loss(m, sqrt(pmax(v, 0)), loss)
  # A decision's contribution to VOI is bounded by its current expected loss,
  # so decisions that are already (almost) certain can be skipped.
  live <- L0 > 1e-9 * max(L0, 1e-300)
  m <- m[live]; v <- v[live]; L0 <- L0[live]; C <- C[live, , drop = FALSE]
  nc <- ncol(C)
  out <- numeric(nc)
  if (!length(L0)) return(out)
  for (s0 in seq(1, nc, by = chunk)) {
    jj <- s0:min(nc, s0 + chunk - 1)
    Cj <- C[, jj, drop = FALSE]
    VR <- sweep(Cj^2, 2, sv_cand[jj] + noise, "/")
    VR <- pmin(VR, v * (1 - 1e-12))
    keep <- rowSums(VR > 1e-8 * pmax(v, 1e-300)) > 0
    if (!any(keep)) next
    VRk <- VR[keep, , drop = FALSE]
    mk <- m[keep]; vk <- v[keep]
    sd_post <- sqrt(pmax(vk - VRk, 0))
    sd_pre <- sqrt(VRk)
    after <- matrix(0, nrow(VRk), ncol(VRk))
    for (q in seq_along(gh$x)) {
      after <- after + gh$w[q] * .min_loss(mk + sd_pre * gh$x[q], sd_post, loss)
    }
    out[jj] <- colSums(L0[keep] - after)
  }
  pmax(out, 0)
}

#' Value-of-information map
#'
#' Computes, for every candidate location, the expected reduction in decision
#' loss from measuring it before deciding (the expected value of sample
#' information). A measurement at one location sharpens predictions at every
#' correlated location, so its value is summed over all decisions it can
#' change. Value is high only where a new measurement could flip a decision,
#' and near zero where decisions are already clear, however uncertain the
#' prediction is.
#'
#' Two methods are available.
#'
#' * `"gaussian"` assumes Gaussian updating: the predictive mean after
#'   sampling is Gaussian around the current mean and the posterior variance is
#'   known in advance. Expectations over the future mean use Gauss-Hermite
#'   quadrature, so fixed and linear losses are handled exactly under this
#'   assumption. For draws, their means and covariance are used.
#' * `"ensemble"` works directly on Monte Carlo draws, such as AI and
#'   machine-learning ensembles, quantile-forest draws or posterior samples. For
#'   each candidate it simulates possible measurements and updates the whole
#'   ensemble with an ensemble Kalman (perturbed-observation) update (Evensen,
#'   1994). Losses are then evaluated on the updated ensemble, so skewed and
#'   multimodal predictive distributions keep their shape. Results carry Monte
#'   Carlo error, which falls as the number of draws grows.
#'
#' `"auto"` uses `"ensemble"` for draws-based predictions with at most 1,000
#' candidates and `"gaussian"` otherwise.
#'
#' @param pred A `dm_predictive` object, ideally with a covariance matrix.
#' @param loss A `dm_loss` object.
#' @param candidates Indices (or a logical vector) of locations that could be
#'   sampled. Default: all locations.
#' @param noise Measurement-error variance of a new observation (for example
#'   the nugget of a fitted covariance model).
#' @param target Optional matrix `H` (one row per decision unit, one column
#'   per location) for decisions about aggregates `Q = H Z`, such as district
#'   means ([dm_block_matrix()]) or watershed loads ([dm_flow_matrix()]). If
#'   `NULL`, one decision is made per location.
#' @param n_quad Number of quadrature nodes (Gaussian method).
#' @param method `"auto"`, `"gaussian"` or `"ensemble"` (see Details).
#' @param n_outer Number of simulated measurement values per candidate
#'   (ensemble method).
#' @param localize Ensemble method: gains are used only where the ensemble
#'   correlation between a location and the candidate measurement exceeds
#'   `localize / sqrt(number of draws)` in absolute value. This removes the
#'   spurious correlations of finite ensembles (covariance localisation). Use
#'   `0` to switch it off.
#' @param seed Optional random seed for the ensemble method; the global random
#'   state is restored.
#' @references Evensen, G. (1994). Sequential data assimilation with a
#'   nonlinear quasi-geostrophic model using Monte Carlo methods to forecast
#'   error statistics. *Journal of Geophysical Research*, 99(C5), 10143-10162.
#'   \doi{10.1029/94JC00572}
#' @return A data frame of class `dm_voi` with columns `site`, `x`, `y` (if
#'   known) and `voi` (expected loss reduction, in loss units), and the
#'   attribute `evpi` (current total expected loss, the upper bound on any
#'   sampling value).
#' @examples
#' sim <- dm_simulate(n_side = 15, n_obs = 30, seed = 2)
#' fit <- dm_krige(sim$obs[, c("x", "y")], sim$obs$z, sim$grid[, c("x", "y")],
#'                 model = sim$model)
#' v <- dm_voi(fit, dm_loss(0.5, "above", ratio = 3), noise = sim$model$nugget)
#' head(v[order(-v$voi), ])
#' plot(v)
#' @export
dm_voi <- function(pred, loss, candidates = NULL, noise = 0, target = NULL, n_quad = 20,
                   method = c("auto", "gaussian", "ensemble"), n_outer = 15, localize = 2,
                   seed = NULL) {
  .check_pred(pred); .check_loss(loss)
  method <- match.arg(method)
  n <- length(pred$mean)
  cand <- .resolve_candidates(candidates, n)
  if (method == "auto") {
    method <- if (!is.null(pred$draws) && length(cand) <= 1000) "ensemble" else "gaussian"
  }
  if (method == "ensemble") {
    if (is.null(pred$draws)) stop("method = 'ensemble' needs draws (see dm_draws()).", call. = FALSE)
    res <- .with_seed(seed, .voi_ensemble(pred$draws, loss, cand, noise, target, n_outer,
                                          localize))
    voi <- res$voi; evpi <- res$evpi
  } else {
    Sigma <- .Sigma(pred)
    mom <- .decision_moments(pred$mean, Sigma, target, cand)
    gh <- .gh_normal(n_quad)
    voi <- .voi_core(mom$m, mom$v, mom$C, diag(Sigma)[cand], noise, loss, gh)
    evpi <- sum(.min_loss(mom$m, sqrt(pmax(mom$v, 0)), loss))
  }
  out <- data.frame(site = cand, voi = voi)
  if (!is.null(pred$coords)) out <- cbind(out[1], pred$coords[cand, , drop = FALSE], out[2])
  rownames(out) <- NULL
  attr(out, "method") <- method
  attr(out, "evpi") <- evpi
  attr(out, "loss") <- loss
  attr(out, "all_coords") <- pred$coords
  class(out) <- c("dm_voi", "data.frame")
  out
}

#' @export
print.dm_voi <- function(x, n = 6, ...) {
  cat("<dm_voi> ", nrow(x), " candidate location(s), ", attr(x, "method") %||% "gaussian",
      " method\n", sep = "")
  cat(sprintf("  current expected loss (EVPI) : %.4g\n", attr(x, "evpi")))
  cat(sprintf("  maximum single-sample value  : %.4g (site %d)\n", max(x$voi),
              x$site[which.max(x$voi)]))
  cat(sprintf("  candidates with value > 1%% of max: %.1f%%\n",
              100 * mean(x$voi > 0.01 * max(x$voi))))
  top <- utils::head(as.data.frame(unclass(x))[order(-x$voi), , drop = FALSE], n)
  top$voi <- signif(top$voi, 4)
  print(top, row.names = FALSE)
  invisible(x)
}

#' @rdname dm_voi
#' @param x A `dm_voi` object.
#' @param main Plot title.
#' @param ... Passed to the internal map plotting function.
#' @export
plot.dm_voi <- function(x, main = "Value of information of one new sample", ...) {
  if (!all(c("x", "y") %in% names(x))) stop("No coordinates to map.", call. = FALSE)
  .map_plot(x[, c("x", "y")], x$voi, main = main, palette = "Inferno", rev = TRUE, ...)
  invisible(x)
}

#' Budget-constrained sampling design by value of information
#'
#' Selects new sampling locations one at a time, each time taking the
#' location with the largest value of information (or value per unit cost),
#' then updating the predictive covariance as if it had been measured. The
#' design stops when the requested number of samples is reached, the budget
#' is spent, or the value of the best remaining sample falls below
#' `stop_below` (for example the cost of one sample expressed in loss units).
#'
#' The design is non-adaptive: it uses the current predictive means
#' throughout, as is standard for planning a survey before any new data are
#' collected. Re-run it after new data arrive.
#'
#' @inheritParams dm_voi
#' @param n Maximum number of samples to select.
#' @param cost Optional cost of sampling each candidate (same length as the
#'   candidates). When given, locations are ranked by value per unit cost.
#' @param budget Optional total budget in the units of `cost`.
#' @param stop_below Stop when the best remaining value (per unit cost, if
#'   `cost` is given) is below this number.
#' @param lazy Use lazy greedy evaluation: values computed earlier are reused
#'   as upper bounds and only the leading candidate is re-evaluated. This is
#'   much faster and exact whenever the value of a candidate does not increase
#'   as other samples are added, which is the usual case; set `FALSE` to
#'   re-evaluate every candidate at every step.
#' @return A data frame of class `dm_design` with the selected sites in order,
#'   their coordinates, marginal value at selection (`voi`), cumulative value
#'   and cost, and the attribute `remaining_loss` (expected loss after each
#'   step).
#' @examples
#' sim <- dm_simulate(n_side = 15, n_obs = 30, seed = 2)
#' fit <- dm_krige(sim$obs[, c("x", "y")], sim$obs$z, sim$grid[, c("x", "y")],
#'                 model = sim$model)
#' des <- dm_design(fit, dm_loss(0.5, "above", ratio = 3), n = 8,
#'                  noise = sim$model$nugget)
#' des
#' plot(des)
#' @export
dm_design <- function(pred, loss, n = 10, candidates = NULL, noise = 0, cost = NULL,
                      budget = Inf, stop_below = 0, target = NULL, n_quad = 20,
                      lazy = TRUE) {
  .check_pred(pred); .check_loss(loss)
  N <- length(pred$mean)
  cand <- .resolve_candidates(candidates, N)
  if (!is.null(cost)) {
    if (length(cost) == 1L) cost <- rep(cost, length(cand))
    if (length(cost) != length(cand) || any(!(cost > 0))) {
      stop("'cost' must be positive with one value per candidate.", call. = FALSE)
    }
  } else cost <- rep(1, length(cand))
  names(cost) <- cand
  Sigma <- .Sigma(pred)
  gh <- .gh_normal(n_quad)
  sel <- integer(0); vals <- numeric(0); spent <- 0
  mom0 <- .decision_moments(pred$mean, Sigma, target, 1L)
  rem <- sum(.min_loss(mom0$m, sqrt(pmax(mom0$v, 0)), loss))
  avail <- cand
  voi_of <- function(ids) {
    mom <- .decision_moments(pred$mean, Sigma, target, ids)
    .voi_core(mom$m, mom$v, mom$C, diag(Sigma)[ids], noise, loss, gh)
  }
  bound <- NULL
  if (lazy) {
    bound <- stats::setNames(voi_of(avail), avail)
  }
  for (k in seq_len(n)) {
    afford <- avail[cost[as.character(avail)] <= budget - spent]
    if (!length(afford)) break
    if (lazy) {
      # Lazy greedy: re-evaluate the candidate with the largest stale value
      # until it stays on top (assumes values only shrink as data are added).
      repeat {
        sc <- bound[as.character(afford)] / cost[as.character(afford)]
        o <- order(sc, decreasing = TRUE)
        s <- afford[o[1]]
        vs <- voi_of(s)
        bound[as.character(s)] <- vs
        nxt <- if (length(o) > 1) sc[o[2]] else -Inf
        if (vs / cost[[as.character(s)]] >= nxt) break
      }
      vb <- vs; sb <- vs / cost[[as.character(s)]]
    } else {
      v <- voi_of(afford)
      score <- v / cost[as.character(afford)]
      b <- which.max(score)
      s <- afford[b]; vb <- v[b]; sb <- score[[b]]
    }
    if (sb <= stop_below) break
    sel <- c(sel, s); vals <- c(vals, vb); spent <- spent + cost[[as.character(s)]]
    cs <- Sigma[, s]
    Sigma <- Sigma - tcrossprod(cs) / (Sigma[s, s] + noise)
    avail <- setdiff(avail, s)
    mom2 <- .decision_moments(pred$mean, Sigma, target, s)
    rem <- c(rem, sum(.min_loss(mom2$m, sqrt(pmax(mom2$v, 0)), loss)))
  }
  out <- data.frame(order = seq_along(sel), site = sel, voi = vals,
                    cum_voi = cumsum(vals), cum_cost = cumsum(cost[as.character(sel)]))
  if (!is.null(pred$coords) && length(sel)) {
    out <- cbind(out[1:2], pred$coords[sel, , drop = FALSE], out[3:5])
  }
  rownames(out) <- NULL
  attr(out, "remaining_loss") <- unname(rem)
  attr(out, "loss") <- loss
  attr(out, "pred") <- pred
  attr(out, "target") <- target
  class(out) <- c("dm_design", "data.frame")
  out
}

#' @export
print.dm_design <- function(x, ...) {
  rl <- attr(x, "remaining_loss")
  cat("<dm_design> ", nrow(x), " sample(s) selected\n", sep = "")
  cat(sprintf("  expected loss before sampling : %.4g\n", rl[1]))
  cat(sprintf("  expected loss after design    : %.4g (%.1f%% reduction)\n",
              rl[length(rl)], 100 * (1 - rl[length(rl)] / rl[1])))
  df <- as.data.frame(unclass(x))
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], signif, 4)
  print(df, row.names = FALSE)
  invisible(x)
}

#' @rdname dm_design
#' @param x A `dm_design` object.
#' @param main Plot title.
#' @param ... Passed to the internal map plotting function.
#' @export
plot.dm_design <- function(x, main = "Value-of-information sampling plan", ...) {
  p <- attr(x, "pred")
  if (is.null(p$coords) || !nrow(x)) stop("Nothing to map.", call. = FALSE)
  if (is.null(attr(x, "target"))) {
    dec <- dm_decide(p, attr(x, "loss"))
    z <- dec$exp_loss; lab <- "expected loss before sampling"
  } else {
    z <- p$sd; lab <- "predictive sd"
  }
  .map_plot(p$coords, z, main = main, palette = "Rocket", rev = TRUE,
            pts = x[, c("x", "y")], pts_labels = x$order,
            legend_note = paste("background:", lab, "| numbers: sampling order"), ...)
  invisible(x)
}
