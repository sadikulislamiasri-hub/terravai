# ---------------------------------------------------------------------------
# Evaluating maps by accuracy and by decision quality
# ---------------------------------------------------------------------------

.actions_for <- function(p, loss, rule) {
  if (is.numeric(p)) {
    t <- loss$threshold
    return(if (loss$direction == "above") p > t else p < t)
  }
  if (rule == "plugin") {
    t <- loss$threshold
    return(if (loss$direction == "above") p$mean > t else p$mean < t)
  }
  el <- .exp_loss(p, loss)
  el$act < el$no
}

.point_of <- function(p) if (is.numeric(p)) p else p$mean

#' Compare maps by accuracy and by decision regret
#'
#' Evaluates competing predictions against validation data with standard
#' accuracy metrics and with realised decision regret: the loss incurred by
#' acting on each map, relative to acting with perfect knowledge of the true
#' values. Repeating this over a range of cost ratios shows whether the most
#' accurate map is also the best map for the decision.
#'
#' Point predictions are turned into actions by the plug-in rule (act if the
#' prediction is beyond the threshold), which is how accuracy-optimised maps
#' are normally used. Probabilistic predictions (`dm_predictive`) use the
#' expected-loss-minimising rule unless `rule = "plugin"`.
#'
#' @param obs True (validation) values.
#' @param predictions A named list. Each element is either a numeric vector of
#'   point predictions or a `dm_predictive` object, aligned with `obs`.
#' @param loss A `dm_loss` object.
#' @param ratios Cost ratios to evaluate. Defaults to the ratio in `loss`.
#' @param rule `"bayes"` (expected-loss-minimising action) or `"plugin"`, for
#'   probabilistic predictions.
#' @return A data frame of class `dm_evaluation` with one row per method and
#'   ratio: `method`, `ratio`, `rmse`, `mae`, `bias`, `regret` (mean realised
#'   loss per location), `wrong` (share of wrong decisions), `missed`
#'   (share of missed actions), `unneeded` (share of unneeded actions),
#'   `rank_rmse` and `rank_regret`.
#' @examples
#' sim <- dm_simulate(n_side = 20, n_obs = 40, seed = 6)
#' fit <- dm_krige(sim$obs[, c("x", "y")], sim$obs$z, sim$grid[, c("x", "y")],
#'                 model = sim$model)
#' smooth <- rep(mean(sim$obs$z), nrow(sim$grid))   # a naive map
#' ev <- dm_evaluate(sim$grid$z, list(kriging = fit, mean_only = smooth),
#'                   dm_loss(0.5, "above"), ratios = c(1, 3, 10))
#' summary(ev)
#' @export
dm_evaluate <- function(obs, predictions, loss, ratios = NULL, rule = c("bayes", "plugin")) {
  .check_loss(loss); rule <- match.arg(rule)
  if (!is.list(predictions) || is.null(names(predictions)) || any(names(predictions) == "")) {
    stop("'predictions' must be a named list.", call. = FALSE)
  }
  obs <- as.numeric(obs)
  for (nm in names(predictions)) {
    if (length(.point_of(predictions[[nm]])) != length(obs)) {
      stop("Prediction '", nm, "' does not match the length of 'obs'.", call. = FALSE)
    }
  }
  ratios <- ratios %||% loss$ratio
  rows <- list()
  for (r in ratios) {
    lr <- dm_loss(loss$threshold, loss$direction, ratio = r, type = loss$type,
                  cost_fp = loss$cost_fp, action = loss$action)
    need <- if (lr$direction == "above") obs > lr$threshold else obs < lr$threshold
    for (nm in names(predictions)) {
      p <- predictions[[nm]]
      pt <- .point_of(p)
      act <- .actions_for(p, lr, rule)
      rows[[length(rows) + 1]] <- data.frame(
        method = nm, ratio = r,
        rmse = sqrt(mean((pt - obs)^2)), mae = mean(abs(pt - obs)), bias = mean(pt - obs),
        regret = mean(.realised_loss(obs, act, lr)),
        wrong = mean(act != need), missed = mean(!act & need), unneeded = mean(act & !need))
    }
  }
  out <- do.call(rbind, rows)
  out$rank_rmse <- stats::ave(out$rmse, out$ratio, FUN = function(v) rank(v, ties.method = "min"))
  out$rank_regret <- stats::ave(out$regret, out$ratio, FUN = function(v) rank(v, ties.method = "min"))
  attr(out, "loss") <- loss
  class(out) <- c("dm_evaluation", "data.frame")
  out
}

#' @export
summary.dm_evaluation <- function(object, ...) {
  s <- lapply(split(object, object$ratio), function(d) {
    best_acc <- d$method[which.min(d$rmse)]
    best_dec <- d$method[which.min(d$regret)]
    tau <- if (nrow(d) > 2) suppressWarnings(stats::cor(d$rmse, d$regret, method = "kendall")) else NA_real_
    data.frame(ratio = d$ratio[1], best_by_rmse = best_acc, best_by_regret = best_dec,
               reversal = best_acc != best_dec,
               extra_regret_pct = 100 * (d$regret[d$method == best_acc] / min(d$regret) - 1),
               kendall_tau = tau)
  })
  out <- do.call(rbind, s)
  rownames(out) <- NULL
  class(out) <- c("summary.dm_evaluation", "data.frame")
  out
}

#' @export
print.summary.dm_evaluation <- function(x, ...) {
  cat("Accuracy versus decision quality\n")
  df <- as.data.frame(unclass(x))
  df$extra_regret_pct <- round(df$extra_regret_pct, 1)
  df$kendall_tau <- round(df$kendall_tau, 2)
  print(df, row.names = FALSE)
  cat(sprintf("\nThe most accurate map was not the best decision map at %d of %d cost ratio(s).\n",
              sum(df$reversal), nrow(df)))
  invisible(x)
}

#' @export
print.dm_evaluation <- function(x, ...) {
  cat("<dm_evaluation> ", length(unique(x$method)), " method(s) x ",
      length(unique(x$ratio)), " cost ratio(s)\n", sep = "")
  df <- as.data.frame(unclass(x))
  num <- vapply(df, is.numeric, logical(1))
  df[num] <- lapply(df[num], signif, 4)
  print(df, row.names = FALSE)
  invisible(x)
}

#' Plot an evaluation: accuracy versus regret across cost ratios
#'
#' @param x A `dm_evaluation` object.
#' @param ... Unused.
#' @return The input, invisibly.
#' @export
plot.dm_evaluation <- function(x, ...) {
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  graphics::par(mfrow = c(1, 2), mar = c(4, 7, 3, 1), mgp = c(2, 0.6, 0), cex.main = 1)
  meth <- unique(x$method)
  cols <- stats::setNames(grDevices::hcl.colors(length(meth), "Dark 3"), meth)
  acc <- x[x$ratio == x$ratio[1], ]
  acc <- acc[order(acc$rmse, decreasing = TRUE), ]
  graphics::barplot(acc$rmse, names.arg = acc$method, horiz = TRUE, las = 1,
                    col = cols[acc$method], border = NA, xlab = "RMSE",
                    main = "Accuracy (lower is better)")
  graphics::par(mar = c(4, 4, 3, 1))
  rs <- sort(unique(x$ratio))
  R <- sapply(meth, function(m) x$regret[x$method == m][match(rs, x$ratio[x$method == m])])
  R <- matrix(R, nrow = length(rs))
  rel <- R / apply(R, 1, min)
  graphics::matplot(rs, rel, type = "b", pch = 16, lty = 1, col = cols[meth], log = if (length(rs) > 1 && min(rs) > 0) "x" else "",
                    xlab = "Cost ratio (missed / unneeded action)", ylab = "Regret relative to best",
                    main = "Decision regret (1 = best)")
  graphics::abline(h = 1, lty = 3, col = "grey50")
  graphics::legend("topleft", legend = meth, col = cols[meth], pch = 16, lty = 1,
                   bty = "n", cex = 0.8)
  invisible(x)
}

#' Decision curves across thresholds
#'
#' Computes the realised regret of each map over a range of action
#' thresholds. For a fixed cost ratio this traces how the ranking of maps
#' depends on the threshold, in the spirit of Murphy diagrams for threshold
#' decisions (Ehm et al., 2016).
#'
#' @inheritParams dm_evaluate
#' @param thresholds Numeric vector of thresholds.
#' @param direction,ratio,type Passed to [dm_loss()].
#' @return A data frame of class `dm_curve` with columns `method`,
#'   `threshold` and `regret`.
#' @references Ehm, W., Gneiting, T., Jordan, A. and Krueger, F. (2016). Of
#'   quantiles and expectiles: consistent scoring functions, Choquet
#'   representations and forecast rankings. *Journal of the Royal Statistical
#'   Society Series B*, 78(3). \doi{10.1111/rssb.12154}
#' @examples
#' sim <- dm_simulate(n_side = 15, n_obs = 30, seed = 7)
#' fit <- dm_krige(sim$obs[, c("x", "y")], sim$obs$z, sim$grid[, c("x", "y")],
#'                 model = sim$model)
#' cur <- dm_threshold_curve(sim$grid$z, list(kriging = fit, plugin = fit$mean),
#'                           thresholds = seq(-1, 1, by = 0.25), ratio = 3)
#' plot(cur)
#' @export
dm_threshold_curve <- function(obs, predictions, thresholds, direction = "above",
                               ratio = 1, type = "fixed", rule = c("bayes", "plugin")) {
  rule <- match.arg(rule)
  rows <- list()
  for (t in thresholds) {
    ev <- dm_evaluate(obs, predictions, dm_loss(t, direction, ratio = ratio, type = type),
                      rule = rule)
    rows[[length(rows) + 1]] <- data.frame(method = ev$method, threshold = t, regret = ev$regret)
  }
  out <- do.call(rbind, rows)
  attr(out, "ratio") <- ratio
  class(out) <- c("dm_curve", "data.frame")
  out
}

#' @rdname dm_threshold_curve
#' @param x A `dm_curve` object.
#' @param ... Unused.
#' @export
plot.dm_curve <- function(x, ...) {
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op))
  graphics::par(mar = c(4, 4, 3, 1), mgp = c(2, 0.6, 0))
  meth <- unique(x$method)
  cols <- grDevices::hcl.colors(length(meth), "Dark 3")
  th <- sort(unique(x$threshold))
  R <- sapply(meth, function(m) x$regret[x$method == m][match(th, x$threshold[x$method == m])])
  R <- matrix(R, nrow = length(th))
  graphics::matplot(th, R, type = "l", lty = 1, lwd = 2, col = cols,
                    xlab = "Action threshold", ylab = "Mean regret per location",
                    main = paste0("Decision curves (cost ratio ", attr(x, "ratio"), ")"))
  graphics::legend("topright", legend = meth, col = cols, lwd = 2, bty = "n", cex = 0.8)
  invisible(x)
}
