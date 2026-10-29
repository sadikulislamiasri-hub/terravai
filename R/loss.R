# ---------------------------------------------------------------------------
# Losses for threshold decisions
# ---------------------------------------------------------------------------

#' Define a threshold decision and its losses
#'
#' A threshold decision chooses between acting (lime the field, target the
#' village, declare the well unsafe) and not acting. Action is *needed* when
#' the true value lies above (or below) a threshold. Two errors are possible:
#' acting when it was not needed (a false positive, cost `cost_fp`) and not
#' acting when it was needed (a false negative, cost `ratio * cost_fp`). The
#' only judgement required from a stakeholder is the cost ratio: how many
#' times worse a missed action is than an unneeded one.
#'
#' With `type = "fixed"` each wrong decision costs a fixed amount. With
#' `type = "linear"` the cost grows with the distance between the true value
#' and the threshold (for example yield foregone or contaminant dose).
#' Losses are computed on the scale of the predictions: if a model works on
#' log values, give the threshold on the log scale.
#'
#' @param threshold Action threshold, on the scale of the predictions.
#' @param direction `"above"` if action is needed when the value exceeds the
#'   threshold (e.g. arsenic), `"below"` if it is needed when the value is under
#'   the threshold (e.g. soil pH for liming).
#' @param ratio Cost of a missed action divided by the cost of an unneeded
#'   action (lambda). Values above 1 mean missing a needed action is worse.
#' @param type `"fixed"` or `"linear"` loss.
#' @param cost_fp Cost of an unneeded action, in any unit (default 1, so
#'   results are in units of this cost).
#' @param action Short label for the action, used in plots and printing.
#' @return An object of class `dm_loss`.
#' @examples
#' # Declare a well unsafe if arsenic > 10 ug/L; missing an unsafe well is
#' # judged 5 times worse than replacing a safe one.
#' dm_loss(10, "above", ratio = 5, action = "declare unsafe")
#' @export
dm_loss <- function(threshold, direction = c("above", "below"), ratio = 1,
                    type = c("fixed", "linear"), cost_fp = 1, action = "act") {
  direction <- match.arg(direction); type <- match.arg(type)
  if (!is.numeric(threshold) || length(threshold) != 1 || !is.finite(threshold)) {
    stop("'threshold' must be a single finite number.", call. = FALSE)
  }
  if (!is.numeric(ratio) || length(ratio) != 1 || !(ratio > 0) || !is.finite(ratio)) {
    stop("'ratio' must be a single positive number.", call. = FALSE)
  }
  if (!(cost_fp > 0)) stop("'cost_fp' must be positive.", call. = FALSE)
  structure(list(threshold = threshold, direction = direction, ratio = ratio,
                 type = type, cost_fp = cost_fp, cost_fn = ratio * cost_fp,
                 action = action), class = "dm_loss")
}

#' @export
print.dm_loss <- function(x, ...) {
  cat("<dm_loss> ", x$action, " when value is ", x$direction, " ", format(x$threshold),
      "\n", sep = "")
  cat("  loss type     : ", x$type, "\n", sep = "")
  cat("  cost ratio    : ", format(x$ratio), " (missed action vs unneeded action)\n", sep = "")
  cat("  act when P(action needed) > ", format(signif(dm_cutoff(x$ratio), 3)),
      if (x$type == "linear") " (fixed-loss equivalent)" else "", "\n", sep = "")
  invisible(x)
}

#' Probability cut-off implied by a cost ratio
#'
#' Under a fixed loss, the expected-loss-minimising rule acts when the
#' probability that action is needed exceeds `1 / (1 + ratio)`. This helper
#' converts between the two, which is a simple way to discuss losses with
#' stakeholders.
#'
#' @param ratio Cost ratio (missed action / unneeded action).
#' @param cutoff Probability cut-off. If supplied, the matching ratio is
#'   returned instead.
#' @return A numeric vector.
#' @examples
#' dm_cutoff(c(1, 4, 9))          # 0.5, 0.2, 0.1
#' dm_cutoff(cutoff = 0.2)        # ratio 4
#' @export
dm_cutoff <- function(ratio = NULL, cutoff = NULL) {
  if (!is.null(cutoff)) return((1 - cutoff) / cutoff)
  1 / (1 + ratio)
}

# Expected losses ----------------------------------------------------------------

# Partial moments for a Gaussian: E[(Z - t)+], E[(t - Z)+]
.pm_gauss <- function(m, s, t) {
  s <- pmax(s, 1e-12)
  z <- (m - t) / s
  list(p_above = stats::pnorm(z),
       up = s * stats::dnorm(z) + (m - t) * stats::pnorm(z),
       dn = s * stats::dnorm(z) + (t - m) * stats::pnorm(-z))
}

# Expected losses of acting / not acting. Returns list(act, no, p_need).
.exp_loss_gauss <- function(m, s, loss) {
  pm <- .pm_gauss(m, s, loss$threshold)
  above <- loss$direction == "above"
  p_need <- if (above) pm$p_above else 1 - pm$p_above
  if (loss$type == "fixed") {
    act <- loss$cost_fp * (1 - p_need)
    no <- loss$cost_fn * p_need
  } else {
    e_need <- if (above) pm$up else pm$dn
    e_not <- if (above) pm$dn else pm$up
    act <- loss$cost_fp * e_not
    no <- loss$cost_fn * e_need
  }
  list(act = act, no = no, p_need = p_need)
}

.exp_loss_draws <- function(D, loss) {
  t <- loss$threshold
  need <- if (loss$direction == "above") D > t else D < t
  p_need <- rowMeans(need)
  if (loss$type == "fixed") {
    act <- loss$cost_fp * (1 - p_need)
    no <- loss$cost_fn * p_need
  } else {
    dist <- abs(D - t)
    act <- loss$cost_fp * rowMeans(dist * !need)
    no <- loss$cost_fn * rowMeans(dist * need)
  }
  list(act = act, no = no, p_need = p_need)
}

.exp_loss <- function(p, loss) {
  if (!is.null(p$draws)) .exp_loss_draws(p$draws, loss) else
    .exp_loss_gauss(p$mean, p$sd, loss)
}

# Realised loss of an action given the true value (perfect information loss = 0)
.realised_loss <- function(y, act, loss) {
  t <- loss$threshold
  need <- if (loss$direction == "above") y > t else y < t
  w <- if (loss$type == "fixed") 1 else abs(y - t)
  loss$cost_fp * w * (act & !need) + loss$cost_fn * w * (!act & need)
}
