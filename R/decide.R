# ---------------------------------------------------------------------------
# Decision maps
# ---------------------------------------------------------------------------

#' Map the optimal decision, its expected loss and its robustness
#'
#' For every location, chooses the action that minimises expected loss under
#' the predictive distribution and reports how confident and how costly that
#' choice is.
#'
#' The expected loss of the optimal action equals the expected value of
#' perfect information at that location: what a perfect measurement would
#' save. The break-even ratio is the cost ratio at which the decision would
#' flip; a decision-maker whose own ratio is above it should act. Mapping it
#' shows how robust each decision is to disagreement about costs.
#'
#' @param pred A `dm_predictive` object.
#' @param loss A `dm_loss` object.
#' @return A data frame of class `dm_decision` with one row per location and
#'   columns `x`, `y` (if coordinates are known), `mean`, `sd`, `p_need`
#'   (probability that action is needed), `act` (optimal action),
#'   `loss_act`, `loss_no` (expected losses of each action), `exp_loss`
#'   (expected loss of the optimal action), `p_wrong` (probability that the
#'   optimal action is wrong) and `breakeven` (cost ratio at which the
#'   decision flips).
#' @examples
#' p <- dm_gaussian(mean = c(4.8, 5.4, 5.6, 6.3), sd = 0.3)
#' dm_decide(p, dm_loss(5.5, "below", ratio = 3, action = "lime"))
#' @export
dm_decide <- function(pred, loss) {
  .check_pred(pred); .check_loss(loss)
  el <- .exp_loss(pred, loss)
  act <- el$act < el$no
  a_unit <- el$act / loss$cost_fp
  n_unit <- el$no / loss$cost_fn
  be <- ifelse(n_unit > 0, a_unit / n_unit, Inf)
  out <- data.frame(mean = pred$mean, sd = pred$sd, p_need = el$p_need, act = act,
                    loss_act = el$act, loss_no = el$no,
                    exp_loss = pmin(el$act, el$no),
                    p_wrong = ifelse(act, 1 - el$p_need, el$p_need),
                    breakeven = be)
  if (!is.null(pred$coords)) out <- cbind(pred$coords, out)
  attr(out, "loss") <- loss
  attr(out, "pred") <- pred[c("coords", "template", "cells")]
  class(out) <- c("dm_decision", "data.frame")
  out
}

#' @export
summary.dm_decision <- function(object, uncertain = 0.2, ...) {
  loss <- attr(object, "loss")
  res <- c(locations = nrow(object),
           share_act = mean(object$act),
           mean_exp_loss = mean(object$exp_loss),
           total_exp_loss = sum(object$exp_loss),
           share_uncertain = mean(object$p_wrong > uncertain),
           share_robust = mean(object$breakeven < loss$ratio / 4 |
                                 object$breakeven > loss$ratio * 4))
  structure(res, uncertain = uncertain, class = c("summary.dm_decision", "numeric"))
}

#' @export
print.summary.dm_decision <- function(x, ...) {
  cat("Decision summary\n")
  cat(sprintf("  locations                       : %d\n", as.integer(x[["locations"]])))
  cat(sprintf("  share where action is optimal   : %.1f%%\n", 100 * x[["share_act"]]))
  cat(sprintf("  mean expected loss (EVPI)       : %.4g\n", x[["mean_exp_loss"]]))
  cat(sprintf("  share with P(wrong) > %.2f      : %.1f%%\n", attr(x, "uncertain"),
              100 * x[["share_uncertain"]]))
  cat(sprintf("  share robust to 4x cost changes : %.1f%%\n", 100 * x[["share_robust"]]))
  invisible(x)
}

#' @export
print.dm_decision <- function(x, n = 6, ...) {
  loss <- attr(x, "loss")
  cat("<dm_decision> ", loss$action, " when value is ", loss$direction, " ",
      format(loss$threshold), " (ratio ", format(loss$ratio), ", ", loss$type,
      " loss)\n", sep = "")
  print(summary(x))
  cat("\n")
  df <- as.data.frame(unclass(x))[seq_len(min(n, nrow(x))), , drop = FALSE]
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], signif, 4)
  print(df, row.names = FALSE)
  if (nrow(x) > n) cat("... ", nrow(x) - n, " more rows\n", sep = "")
  invisible(x)
}

#' Plot a decision map
#'
#' @param x A `dm_decision` object.
#' @param what Which layer to map: `"act"` (optimal action), `"p_need"`,
#'   `"p_wrong"`, `"exp_loss"`, `"breakeven"` (on a log scale) or `"mean"`.
#' @param main Plot title.
#' @param ... Passed to the internal map plotting function.
#' @return The input, invisibly.
#' @examples
#' sim <- dm_simulate(n_side = 20, n_obs = 40, seed = 1)
#' fit <- dm_krige(sim$obs[, c("x", "y")], sim$obs$z, sim$grid[, c("x", "y")],
#'                 model = sim$model)
#' dec <- dm_decide(fit, dm_loss(0.5, "above", ratio = 3))
#' plot(dec, "p_wrong")
#' @export
plot.dm_decision <- function(x, what = c("act", "p_need", "p_wrong", "exp_loss",
                                         "breakeven", "mean"), main = NULL, ...) {
  what <- match.arg(what)
  if (!all(c("x", "y") %in% names(x))) stop("No coordinates to map.", call. = FALSE)
  loss <- attr(x, "loss")
  coords <- x[, c("x", "y")]
  if (what == "act") {
    .map_plot(coords, as.integer(x$act), discrete = TRUE,
              labels = c(paste("no", loss$action), loss$action),
              cols = c("#D9D9D9", "#B23A48"),
              main = main %||% paste("Optimal decision:", loss$action), ...)
  } else if (what == "breakeven") {
    v <- log10(pmin(pmax(x$breakeven, 1e-3), 1e3))
    centre <- log10(loss$ratio)
    .map_plot(coords, v, main = main %||% "Break-even cost ratio (log10)",
              palette = "Blue-Red 3", rev = TRUE, zlim = centre + c(-2, 2),
              legend_note = paste0("act where your log10(ratio) exceeds the value; current ",
                                   format(round(centre, 2)), " (white)"), ...)
  } else {
    titles <- c(p_need = "P(action needed)", p_wrong = "P(decision is wrong)",
                exp_loss = "Expected loss of optimal decision", mean = "Predicted mean")
    pal <- c(p_need = "YlOrRd", p_wrong = "Reds", exp_loss = "Rocket", mean = "Viridis")
    rev <- what %in% c("exp_loss", "p_need", "p_wrong")
    zl <- if (what %in% c("p_need", "p_wrong")) c(0, if (what == "p_wrong") 0.5 else 1) else NULL
    .map_plot(coords, x[[what]], main = main %||% titles[[what]], palette = pal[[what]],
              rev = rev, zlim = zl, ...)
  }
  invisible(x)
}

.check_pred <- function(p) {
  if (!inherits(p, "dm_predictive")) {
    stop("'pred' must be a dm_predictive object (see dm_gaussian(), dm_draws()).",
         call. = FALSE)
  }
}
.check_loss <- function(l) {
  if (!inherits(l, "dm_loss")) stop("'loss' must be a dm_loss object (see dm_loss()).",
                                    call. = FALSE)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
