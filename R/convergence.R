# ---------------------------------------------------------------------------
# Monte Carlo convergence diagnostics
# ---------------------------------------------------------------------------

#' Check whether enough Monte Carlo draws were used
#'
#' Results from uncertain inputs ([dm_uncertain()]), simulations or AI
#' ensembles are Monte Carlo estimates. This function reports their Monte
#' Carlo standard error (MCSE) for the mean and, if a decision is given, for
#' the probability that action is needed, shows running estimates as draws
#' accumulate, identifies locations whose decision could still change
#' because of Monte Carlo noise, and estimates how many draws would reach a
#' target precision. The MCSE of a probability \eqn{p} from \eqn{M}
#' independent draws is \eqn{\sqrt{p(1-p)/M}}; of a mean, \eqn{s/\sqrt{M}}.
#'
#' @param x A `dm_predictive` object with draws.
#' @param loss Optional `dm_loss`; adds the exceedance probability and
#'   decision stability.
#' @param target_mcse_p Target MCSE for probabilities (default 0.01).
#' @param target_rel_mean Target MCSE of the mean, relative to the
#'   predictive standard deviation (default 0.05).
#' @return An object of class `dm_convergence` (print and plot it) with
#'   per-location MCSEs, running estimates and the recommended number of
#'   draws.
#' @examples
#' A <- dm_rusle(R = dm_uncertain(800, cv = 0.2, lower = 0), K = 0.03, LS = 2,
#'               C = dm_uncertain(0.25, sd = 0.07, lower = 0, upper = 1),
#'               n_draws = 400, seed = 1)
#' cv <- dm_convergence(A, dm_loss_tolerance())
#' cv
#' plot(cv)
#' @export
dm_convergence <- function(x, loss = NULL, target_mcse_p = 0.01, target_rel_mean = 0.05) {
  .check_pred(x)
  if (is.null(x$draws)) stop("'x' has no draws; convergence applies to Monte Carlo results.",
                             call. = FALSE)
  D <- x$draws
  M <- ncol(D)
  s <- apply(D, 1, stats::sd)
  mcse_mean <- s / sqrt(M)
  steps <- unique(round(exp(seq(log(min(20, M)), log(M), length.out = 30))))
  cm <- t(apply(D, 1, cumsum))
  if (nrow(D) == 1L) cm <- matrix(cm, nrow = 1)
  run_mean <- sweep(cm[, steps, drop = FALSE], 2, steps, "/")
  out <- list(M = M, mcse_mean = mcse_mean, rel_mcse_mean = mcse_mean / pmax(s, 1e-300),
              steps = steps, run_mean = run_mean, loss = loss)
  need_mean <- ceiling(M * max(out$rel_mcse_mean / target_rel_mean)^2)
  need_p <- NA_integer_
  if (!is.null(loss)) {
    .check_loss(loss)
    el <- .exp_loss_draws(D, loss)
    p <- el$p_need
    mcse_p <- sqrt(p * (1 - p) / M)
    ind <- if (loss$direction == "above") D > loss$threshold else D < loss$threshold
    cp <- t(apply(ind, 1, cumsum))
    if (nrow(D) == 1L) cp <- matrix(cp, nrow = 1)
    out$run_p <- sweep(cp[, steps, drop = FALSE], 2, steps, "/")
    out$p <- p; out$mcse_p <- mcse_p
    if (loss$type == "fixed") {
      cut <- dm_cutoff(loss$ratio)
      out$unstable <- abs(p - cut) < 2 * mcse_p & mcse_p > 0
    } else {
      diff_l <- el$act - el$no
      se_d <- sqrt(apply(.realised_diff(D, loss), 1, stats::var) / M)
      out$unstable <- abs(diff_l) < 2 * se_d & se_d > 0
    }
    need_p <- ceiling(M * max(mcse_p / target_mcse_p)^2)
  }
  out$recommended <- max(need_mean, need_p, na.rm = TRUE)
  out$targets <- c(p = target_mcse_p, rel_mean = target_rel_mean)
  class(out) <- "dm_convergence"
  out
}

# Per-draw difference between the losses of acting and not acting
.realised_diff <- function(D, loss) {
  t <- loss$threshold
  need <- if (loss$direction == "above") D > t else D < t
  w <- abs(D - t)
  loss$cost_fp * w * !need - loss$cost_fn * w * need
}

#' @export
print.dm_convergence <- function(x, ...) {
  cat("<dm_convergence> ", x$M, " draws, ", length(x$mcse_mean), " location(s)\n", sep = "")
  cat(sprintf("  MCSE of mean / predictive sd : max %.3f (target %.3f)\n",
              max(x$rel_mcse_mean), x$targets[["rel_mean"]]))
  if (!is.null(x$loss)) {
    cat(sprintf("  MCSE of P(action needed)     : max %.4f (target %.4f)\n",
                max(x$mcse_p), x$targets[["p"]]))
    cat(sprintf("  decisions within Monte Carlo noise of the switch point: %d (%.1f%%)\n",
                sum(x$unstable), 100 * mean(x$unstable)))
  }
  ok <- x$recommended <= x$M
  cat("  draws recommended for the targets: ", x$recommended,
      if (ok) " -> current draws are sufficient" else " -> increase n_draws", "\n", sep = "")
  invisible(x)
}

#' @rdname dm_convergence
#' @param ... Unused.
#' @export
plot.dm_convergence <- function(x, ...) {
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  has_p <- !is.null(x$run_p)
  graphics::par(mfrow = c(1, if (has_p) 2 else 1), mar = c(4, 4, 3, 1), mgp = c(2, 0.6, 0))
  pick <- if (has_p) order(-x$mcse_p)[seq_len(min(5, length(x$mcse_p)))] else
    order(-x$rel_mcse_mean)[seq_len(min(5, length(x$mcse_mean)))]
  cols <- grDevices::hcl.colors(length(pick), "Dark 3")
  rm <- x$run_mean[pick, , drop = FALSE]
  graphics::matplot(x$steps, t(rm), type = "l", lty = 1, lwd = 2, col = cols, log = "x",
                    xlab = "Number of draws", ylab = "Running mean",
                    main = "Running mean (least stable locations)")
  if (has_p) {
    rp <- x$run_p[pick, , drop = FALSE]
    graphics::matplot(x$steps, t(rp), type = "l", lty = 1, lwd = 2, col = cols, log = "x",
                      ylim = c(0, 1), xlab = "Number of draws", ylab = "Running P(action needed)",
                      main = "Running exceedance probability")
    if (x$loss$type == "fixed") graphics::abline(h = dm_cutoff(x$loss$ratio), lty = 2,
                                                 col = "grey40")
  }
  invisible(x)
}
