# ---------------------------------------------------------------------------
# Indian standards: particle-size systems, soil-test ratings, erosivity,
# slope classes and erosion classes
# ---------------------------------------------------------------------------

#' Convert particle-size data from any system to USDA sand, silt and clay
#'
#' Soil laboratories report particle-size fractions under different systems:
#' USDA (clay < 0.002 mm, silt 0.002-0.05 mm, sand 0.05-2 mm, as used by
#' USDA Soil Taxonomy), the International (Atterberg/ISSS) system (silt
#' 0.002-0.02 mm, sand 0.02-2 mm) used by many laboratories, or the Indian
#' Standard engineering classification IS 1498 (silt 0.002-0.075 mm, fine
#' sand 0.075-0.425 mm). USDA texture classes ([dm_texture_class()]) and the
#' K factor ([dm_k_epic()]) need USDA fractions. This function converts any
#' cumulative particle-size distribution by log-linear interpolation of the
#' cumulative percentage finer between the measured diameters.
#'
#' @param cum_pct Matrix or data frame of cumulative percentage finer than
#'   each diameter in `diam_mm` (one row per sample, columns in increasing
#'   diameter). The largest diameter should be 2 mm (100 \% for the fine
#'   earth).
#' @param diam_mm Diameters (mm) matching the columns of `cum_pct`; must
#'   include 0.002.
#' @return A data frame with USDA `sand`, `silt` and `clay` (\%).
#' @seealso [dm_psd_isss()] and [dm_psd_is1498()] for the common Indian
#'   laboratory layouts.
#' @examples
#' # cumulative % finer than 0.002, 0.02, 0.2 and 2 mm
#' dm_psd_to_usda(rbind(c(20, 45, 80, 100)), c(0.002, 0.02, 0.2, 2))
#' @export
dm_psd_to_usda <- function(cum_pct, diam_mm) {
  P <- as.matrix(cum_pct)
  if (ncol(P) != length(diam_mm)) stop("ncol(cum_pct) must equal length(diam_mm).", call. = FALSE)
  if (is.unsorted(diam_mm, strictly = TRUE)) stop("'diam_mm' must be increasing.", call. = FALSE)
  if (!any(abs(diam_mm - 0.002) < 1e-9)) stop("'diam_mm' must include 0.002 mm (clay limit).",
                                              call. = FALSE)
  if (any(P < 0 | P > 100.5, na.rm = TRUE)) stop("Percentages must lie between 0 and 100.",
                                                 call. = FALSE)
  if (any(t(apply(P, 1, diff)) < -1e-8, na.rm = TRUE)) {
    stop("Cumulative percentages must not decrease with diameter.", call. = FALSE)
  }
  if (max(diam_mm) < 0.05) stop("Need a diameter of at least 0.05 mm.", call. = FALSE)
  lx <- log10(diam_mm)
  p50 <- apply(P, 1, function(p) stats::approx(lx, p, xout = log10(0.05))$y)
  clay <- P[, which.min(abs(diam_mm - 0.002))]
  top <- P[, ncol(P)]
  out <- data.frame(sand = top - p50, silt = p50 - clay, clay = clay) * (100 / top)
  rownames(out) <- NULL
  out
}

#' Particle-size conversion for International (ISSS) laboratory data
#'
#' Converts International-system fractions (clay < 0.002 mm, silt
#' 0.002-0.02 mm, fine sand 0.02-0.2 mm, coarse sand 0.2-2 mm) to USDA
#' fractions. If only total sand is known, interpolation spans 0.02-2 mm and
#' is less accurate (a warning is given).
#'
#' @param clay,silt Clay and International silt (\%).
#' @param fine_sand,coarse_sand Fine sand (0.02-0.2 mm) and coarse sand
#'   (0.2-2 mm) (\%), or supply total `sand`.
#' @param sand Total sand (0.02-2 mm) when the split is not available.
#' @return A data frame with USDA `sand`, `silt` and `clay` (\%).
#' @examples
#' dm_psd_isss(clay = 25, silt = 15, fine_sand = 35, coarse_sand = 25)
#' @export
dm_psd_isss <- function(clay, silt, fine_sand = NULL, coarse_sand = NULL, sand = NULL) {
  if (!is.null(fine_sand) && !is.null(coarse_sand)) {
    cum <- cbind(clay, clay + silt, clay + silt + fine_sand, clay + silt + fine_sand + coarse_sand)
    return(dm_psd_to_usda(cum, c(0.002, 0.02, 0.2, 2)))
  }
  if (is.null(sand)) stop("Supply 'fine_sand' and 'coarse_sand', or 'sand'.", call. = FALSE)
  warning("Without the fine/coarse sand split, the 0.02-0.05 mm fraction is interpolated ",
          "across 0.02-2 mm and is less accurate.", call. = FALSE)
  dm_psd_to_usda(cbind(clay, clay + silt, clay + silt + sand), c(0.002, 0.02, 2))
}

#' Particle-size conversion for Indian Standard IS 1498 data
#'
#' Converts fractions of the Indian Standard engineering soil classification
#' (clay < 0.002 mm, silt 0.002-0.075 mm, fine sand 0.075-0.425 mm, medium
#' sand 0.425-2 mm; coarser material excluded) to USDA fractions.
#'
#' @param clay,silt,fine_sand,medium_sand Fractions (\%) of the fine earth.
#' @return A data frame with USDA `sand`, `silt` and `clay` (\%).
#' @references Bureau of Indian Standards (1970). *IS 1498: Classification
#'   and Identification of Soils for General Engineering Purposes*.
#' @examples
#' dm_psd_is1498(clay = 20, silt = 40, fine_sand = 25, medium_sand = 15)
#' @export
dm_psd_is1498 <- function(clay, silt, fine_sand, medium_sand) {
  cum <- cbind(clay, clay + silt, clay + silt + fine_sand,
               clay + silt + fine_sand + medium_sand)
  dm_psd_to_usda(cum, c(0.002, 0.075, 0.425, 2))
}

#' Soil Health Card rating limits (India)
#'
#' Soil-test rating limits used for soil health cards in India: organic
#' carbon and available nitrogen, phosphorus and potassium as low, medium or
#' high, and sulphur and micronutrients as deficient or sufficient, with the
#' extraction methods they refer to. Supply your own table with the same
#' columns for state-specific ratings.
#'
#' @return A data frame with columns `parameter`, `unit`, `method`, `low`
#'   (upper limit of the low or deficient class) and `high` (lower limit of
#'   the high class; `NA` for two-class parameters).
#' @references ICAR-Indian Institute of Soil Science, Bhopal. Soil health card
#'   rating limits, *IISS e-Magazine*, 4(1).
#'   <https://iiss.icar.gov.in/eMagazine/v4i1/12.pdf>
#' @examples
#' dm_shc_table()
#' @export
dm_shc_table <- function() {
  data.frame(
    parameter = c("OC", "N", "P", "K", "S", "Zn", "Cu", "Fe", "Mn", "B"),
    unit = c("percent", "kg/ha", "kg/ha", "kg/ha", rep("mg/kg", 6)),
    method = c("Walkley-Black", "Alkaline KMnO4", "Olsen", "1 N NH4OAc", "CaCl2",
               "DTPA", "DTPA", "DTPA", "DTPA", "Hot water"),
    low = c(0.50, 280, 10, 120, 10, 0.6, 0.2, 4.5, 2.0, 0.5),
    high = c(0.75, 560, 25, 280, NA, NA, NA, NA, NA, NA),
    stringsAsFactors = FALSE)
}

#' Rate soil-test values (Soil Health Card)
#'
#' Classifies soil-test values as low, medium or high (organic carbon, N,
#' P, K) or deficient or sufficient (S and micronutrients) using
#' [dm_shc_table()] or a user table. With an uncertain input (a
#' `dm_predictive`, e.g. from kriging of soil-test data), the probability of
#' each class is returned instead, which shows where a rating is reliable.
#'
#' @param x Values (numeric) or a `dm_predictive` object.
#' @param parameter One of the parameters in the table, e.g. `"N"`.
#' @param table Rating table with the columns of [dm_shc_table()].
#' @return A factor of classes, or for a `dm_predictive` a data frame of class
#'   probabilities.
#' @examples
#' dm_shc_rating(c(180, 350, 620), "N")
#' dm_shc_rating(c(0.3, 1.2), "Zn")
#' dm_shc_rating(dm_gaussian(c(9, 18, 30), sd = 3), "P")
#' @export
dm_shc_rating <- function(x, parameter, table = dm_shc_table()) {
  row <- table[table$parameter == parameter, ]
  if (nrow(row) != 1) stop("Unknown parameter '", parameter, "'.", call. = FALSE)
  two <- is.na(row$high)
  labs <- if (two) c("deficient", "sufficient") else c("low", "medium", "high")
  brk <- if (two) c(-Inf, row$low, Inf) else c(-Inf, row$low, row$high, Inf)
  if (inherits(x, "dm_predictive")) {
    F <- if (!is.null(x$draws)) {
      sapply(brk[-c(1, length(brk))], function(b) rowMeans(x$draws < b))
    } else sapply(brk[-c(1, length(brk))], function(b) stats::pnorm((b - x$mean) / pmax(x$sd, 1e-12)))
    F <- cbind(0, matrix(F, nrow = length(x$mean)), 1)
    out <- as.data.frame(t(apply(F, 1, diff)))
    names(out) <- paste0("p_", labs)
    return(out)
  }
  factor(labs[findInterval(x, brk[-1], left.open = FALSE) + 1], levels = labs)
}

#' Rainfall erosivity for India from annual rainfall
#'
#' Estimates the RUSLE rainfall erosivity factor from mean annual rainfall
#' with the relation developed for India by Singh, Chandra and Babu (1981):
#' \eqn{R = 79 + 0.363 P}, with \eqn{P} the mean annual rainfall (mm) and
#' \eqn{R} in MJ mm ha\eqn{^{-1}} h\eqn{^{-1}} yr\eqn{^{-1}}. Inputs can be
#' numbers, rasters (e.g. gridded rainfall) or uncertain values.
#'
#' @param P Mean annual rainfall (mm).
#' @param a,b Coefficients (defaults from Singh et al., 1981).
#' @param n_draws,seed Monte Carlo settings used when an input is uncertain.
#' @return Erosivity R.
#' @references Singh, G., Chandra, S. and Babu, R. (1981). Central Soil and
#'   Water Conservation Research and Training Institute (now ICAR-Indian
#'   Institute of Soil and Water Conservation), Dehradun. As cited in Indian
#'   erosivity studies, e.g. *Indian Journal of Agricultural Sciences*:
#'   <https://epubs.icar.org.in/index.php/IJAgS/article/view/104246>
#' @examples
#' dm_erosivity_india(c(800, 1400, 2200))
#' @export
dm_erosivity_india <- function(P, a = 79, b = 0.363, n_draws = 500, seed = NULL) {
  .check_range(P, "P", 0)
  .apply_formula(function(P) a + b * P, list(P = P), n_draws, seed)
}

#' Slope classes used in Indian soil and land-use surveys
#'
#' Classifies slopes into the classes used by the All India Soil and Land Use
#' Survey (AIS&LUS): A (0-1 \%, nearly level), B (1-3 \%, very gently
#' sloping), C (3-5 \%, gently sloping), D (5-10 \%, moderately sloping),
#' E (10-15 \%, strongly sloping), F (15-25 \%, moderately steep),
#' G (25-33 \%, steep), H (33-50 \%, very steep) and I (above 50 \%, very
#' very steep).
#'
#' @param slope_pct Slope (\%). Convert from degrees with
#'   `100 * tan(deg * pi / 180)`.
#' @param labels Return descriptive labels (`TRUE`) or class codes.
#' @return A factor of slope classes.
#' @examples
#' dm_slope_class_india(c(0.5, 4, 12, 40, 70))
#' @export
dm_slope_class_india <- function(slope_pct, labels = FALSE) {
  if (any(slope_pct < 0, na.rm = TRUE)) stop("Slopes cannot be negative.", call. = FALSE)
  codes <- LETTERS[1:9]
  desc <- c("nearly level", "very gently sloping", "gently sloping", "moderately sloping",
            "strongly sloping", "moderately steep", "steep", "very steep", "very very steep")
  i <- findInterval(slope_pct, c(0, 1, 3, 5, 10, 15, 25, 33, 50), rightmost.closed = TRUE)
  lev <- if (labels) paste0(codes, " (", desc, ")") else codes
  factor(lev[i], levels = lev)
}

#' Soil erosion severity classes
#'
#' Groups soil-loss estimates into severity classes. The default scheme
#' (`"india_usle"`) uses classes applied in USLE-based erosion mapping in
#' India: E0 nil (0-5), E1 slight (5-12), E2 moderate (12-25), E3 severe
#' (25-60) and E4 very severe (above 60 t ha\eqn{^{-1}} yr\eqn{^{-1}}).
#' Schemes differ between studies and agencies, so supply `breaks` and
#' `labels` to use another. With a `dm_predictive` input, the probability of
#' each class is returned.
#'
#' @param x Soil loss (t ha\eqn{^{-1}} yr\eqn{^{-1}}), numeric or
#'   `dm_predictive`.
#' @param breaks Upper limits of all but the last class.
#' @param labels Class labels (one more than `breaks`).
#' @return A factor, or a data frame of class probabilities.
#' @examples
#' dm_erosion_class(c(3, 8, 20, 45, 90))
#' @export
dm_erosion_class <- function(x, breaks = c(5, 12, 25, 60),
                             labels = c("E0 nil", "E1 slight", "E2 moderate", "E3 severe",
                                        "E4 very severe")) {
  if (length(labels) != length(breaks) + 1) stop("Need one more label than breaks.", call. = FALSE)
  if (is.unsorted(breaks, strictly = TRUE)) stop("'breaks' must be increasing.", call. = FALSE)
  if (inherits(x, "dm_predictive")) {
    F <- if (!is.null(x$draws)) sapply(breaks, function(b) rowMeans(x$draws < b)) else
      sapply(breaks, function(b) stats::pnorm((b - x$mean) / pmax(x$sd, 1e-12)))
    F <- cbind(0, matrix(F, nrow = length(x$mean)), 1)
    out <- as.data.frame(t(apply(F, 1, diff)))
    names(out) <- paste0("p_", make.names(labels))
    return(out)
  }
  factor(labels[findInterval(x, breaks) + 1], levels = labels)
}
