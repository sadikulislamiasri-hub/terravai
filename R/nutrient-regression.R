# ---------------------------------------------------------------------------
# Empirical N, P, K loss equations (piecewise linear), ICAR-IISWC data
# Coefficients from the SoilRunNutri Python package (MIT licence;
# see inst/COPYRIGHTS).
# ---------------------------------------------------------------------------

.nutrient_coef_table <- function() {
  data.frame(
    nutrient = rep(c("N", "P", "K"), each = 8),
    basis = rep(rep(c("soil", "runoff"), each = 4), 3),
    lower = rep(c(0, 5, 10, 20, 0, 10, 20, 30), 3),
    upper = rep(c(5, 10, 20, Inf, 10, 20, 30, Inf), 3),
    slope = c(2.3647, 1.2927, 0.7909, 0.8012,  4.285, 0.415, 0.0679, 0.0662,
              0.2096, 0.0853, 0.1048, 0.0622,  0.189, 0.0478, 0.027, 0.0262,
              3.1155, 1.2188, 1.0204, 0.5914,  2.289, 0.4252, 0.4512, 0.4361),
    intercept = c(0, 0.270, -0.2909, 1.1358,  0, 0.246, -0.0175, 0.0144,
                  0, -0.0103, 0.0291, 0.0794,  0, 0.018, -0.0066, 0.0046,
                  0, 0.2178, 2.5149, 0.4662,  0, -0.1895, -0.3018, 0.262),
    stringsAsFactors = FALSE)
}

#' Coefficients of the empirical nutrient-loss equations
#'
#' Returns the piecewise linear coefficients used by
#' [dm_nutrient_empirical()], with the input range each segment applies to.
#' For each segment, nutrient loss (kg ha\eqn{^{-1}}) = slope x input +
#' intercept, where the input is soil loss (t ha\eqn{^{-1}}) or runoff
#' (\% of rainfall).
#'
#' @return A data frame with columns `nutrient`, `basis`, `lower`, `upper`,
#'   `slope`, `intercept`, and the attribute `jumps` giving the size of the
#'   discontinuity at each internal breakpoint.
#' @examples
#' dm_nutrient_coef()
#' @export
dm_nutrient_coef <- function() {
  tab <- .nutrient_coef_table()
  jumps <- do.call(rbind, lapply(split(tab, paste(tab$nutrient, tab$basis)), function(d) {
    d <- d[order(d$lower), ]
    b <- d$upper[-nrow(d)]
    left <- d$slope[-nrow(d)] * b + d$intercept[-nrow(d)]
    right <- d$slope[-1] * b + d$intercept[-1]
    data.frame(nutrient = d$nutrient[1], basis = d$basis[1], breakpoint = b,
               left = left, right = right, ratio = left / right)
  }))
  rownames(jumps) <- NULL
  attr(tab, "jumps") <- jumps
  tab
}

#' Empirical nitrogen, phosphorus and potassium loss
#'
#' Estimates erosion-induced loss of N, P or K (kg ha\eqn{^{-1}}) from soil
#' loss (t ha\eqn{^{-1}}) or from runoff (\% of rainfall), using piecewise
#' linear regressions fitted to long-term runoff-plot experiments of the
#' ICAR-Indian Institute of Soil and Water Conservation, Dehradun, across
#' several agro-ecological regions of India. The equations were first
#' released in the Python package 'SoilRunNutri' (Pant, Islam, Kumar and
#' Madhu); see [dm_nutrient_coef()] for the coefficients.
#'
#' The function is vectorised, works on rasters, and propagates uncertainty:
#' if the input is a [dm_uncertain()] or `dm_predictive` object (for example
#' soil loss from [dm_rusle()] or an AI model), the output is a predictive
#' distribution of nutrient loss that can be passed to [dm_decide()] and
#' [dm_voi()].
#'
#' @section Caution:
#' The published segments are not continuous at their breakpoints: the
#' predicted loss can change abruptly when the input crosses 5, 10 or 20
#' t ha\eqn{^{-1}} of soil loss, or 10, 20 or 30\% runoff (for nitrogen with
#' runoff, by about a factor of ten at 10\%). Inspect
#' `attr(dm_nutrient_coef(), "jumps")` and treat estimates near breakpoints
#' with care; formula-based alternatives are [dm_nutrient_sediment()] and
#' [dm_nutrient_runoff()].
#'
#' @param x Soil loss (t ha\eqn{^{-1}}) or runoff (\%): numeric, `SpatRaster`,
#'   [dm_uncertain()] or `dm_predictive`.
#' @param nutrient `"N"`, `"P"` or `"K"`.
#' @param basis `"soil"` (input is soil loss) or `"runoff"` (input is runoff
#'   percentage).
#' @param n_draws,seed Monte Carlo settings used when the input is uncertain.
#' @return Nutrient loss (kg ha\eqn{^{-1}}), in the form of the input.
#' @examples
#' dm_nutrient_empirical(c(3, 8, 15, 30), "N", "soil")
#' dm_nutrient_empirical(25, "P", "runoff")
#' A <- dm_uncertain(14, sd = 3, lower = 0)
#' dm_nutrient_empirical(A, "K", "soil", seed = 1)
#' @export
dm_nutrient_empirical <- function(x, nutrient = c("N", "P", "K"),
                                  basis = c("soil", "runoff"), n_draws = 500, seed = NULL) {
  nutrient <- match.arg(nutrient); basis <- match.arg(basis)
  .check_range(x, "x", 0, if (basis == "runoff") 100 else Inf)
  tab <- .nutrient_coef_table()
  tab <- tab[tab$nutrient == nutrient & tab$basis == basis, ]
  brk <- c(tab$lower, Inf)
  f <- function(x) {
    if (any(x < 0, na.rm = TRUE)) stop("Inputs cannot be negative.", call. = FALSE)
    i <- findInterval(x, brk, rightmost.closed = FALSE)
    i[i < 1] <- 1
    out <- tab$slope[i] * x + tab$intercept[i]
    dim(out) <- dim(x)
    pmax(out, 0)
  }
  .apply_formula(f, list(x = x), n_draws, seed)
}
