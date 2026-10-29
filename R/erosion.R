# ---------------------------------------------------------------------------
# Formula-based runoff, soil loss and nutrient loss
# Every function accepts numbers, SpatRasters or uncertain inputs
# (dm_uncertain / dm_predictive). Uncertain inputs give a dm_predictive result.
# ---------------------------------------------------------------------------

#' Runoff depth by the SCS Curve Number method
#'
#' Computes direct runoff depth from rainfall with the USDA Soil Conservation
#' Service (now NRCS) Curve Number method:
#' \deqn{Q = (P - I_a)^2 / (P - I_a + S) \textrm{ for } P > I_a, \textrm{ else } 0,}
#' with potential retention \eqn{S = 25400 / CN - 254} (mm) and initial
#' abstraction \eqn{I_a = \lambda S}.
#'
#' @param P Rainfall depth (mm).
#' @param CN Curve number (between 1 and 100).
#' @param lambda Initial abstraction ratio (0.2 in the original method; 0.05
#'   is often recommended).
#' @param n_draws,seed Monte Carlo settings used when an input is uncertain.
#' @return Runoff depth (mm), in the form of the inputs: numeric, `SpatRaster`
#'   or `dm_predictive`.
#' @references USDA-SCS (1972). *National Engineering Handbook, Section 4:
#'   Hydrology*. Soil Conservation Service, Washington, DC.
#' @examples
#' dm_runoff_cn(P = c(20, 60, 120), CN = 80)
#' dm_runoff_cn(P = 60, CN = dm_uncertain(80, sd = 5, lower = 1, upper = 100), seed = 1)
#' @export
dm_runoff_cn <- function(P, CN, lambda = 0.2, n_draws = 500, seed = NULL) {
  .check_range(P, "P", 0); .check_range(CN, "CN", 1, 100)
  if (!(lambda >= 0 && lambda <= 1)) stop("'lambda' must lie in [0, 1].", call. = FALSE)
  f <- function(P, CN) {
    CN <- pmin(pmax(CN, 1), 100)
    S <- 25400 / CN - 254
    Ia <- lambda * S
    ifelse(P > Ia, (P - Ia)^2 / (P - Ia + S), 0)
  }
  .apply_formula(f, list(P = P, CN = CN), n_draws, seed)
}

#' Adjust a curve number for antecedent moisture
#'
#' Converts a curve number for average conditions (AMC II) to dry (AMC I) or
#' wet (AMC III) conditions:
#' \eqn{CN_I = 4.2 CN_{II} / (10 - 0.058 CN_{II})},
#' \eqn{CN_{III} = 23 CN_{II} / (10 + 0.13 CN_{II})}.
#'
#' @param CN Curve number for AMC II.
#' @param amc Target condition: `"I"` (dry) or `"III"` (wet).
#' @return Adjusted curve number.
#' @references Chow, V. T., Maidment, D. R. and Mays, L. W. (1988).
#'   *Applied Hydrology*. McGraw-Hill, New York.
#' @examples
#' dm_cn_amc(c(60, 75, 90), "III")
#' @export
dm_cn_amc <- function(CN, amc = c("I", "III")) {
  amc <- match.arg(amc)
  .check_range(CN, "CN", 1, 100)
  f_amc <- function(CN) if (amc == "I") 4.2 * CN / (10 - 0.058 * CN) else 23 * CN / (10 + 0.13 * CN)
  .apply_formula(f_amc, list(CN = CN))
}

#' Peak runoff rate by the rational method
#'
#' \eqn{q_p = C\, i\, A / 360}, with \eqn{q_p} in m3/s, rainfall intensity
#' \eqn{i} in mm/h and area \eqn{A} in hectares.
#'
#' @param C Runoff coefficient (0 to 1).
#' @param i Rainfall intensity (mm/h) for a duration equal to the time of
#'   concentration.
#' @param A Catchment area (ha).
#' @param n_draws,seed Monte Carlo settings used when an input is uncertain.
#' @return Peak discharge (m3/s).
#' @examples
#' dm_peak_rational(C = 0.4, i = 60, A = 25)
#' @export
dm_peak_rational <- function(C, i, A, n_draws = 500, seed = NULL) {
  .check_range(C, "C", 0, 1); .check_range(i, "i", 0); .check_range(A, "A", 0)
  .apply_formula(function(C, i, A) C * i * A / 360, list(C = C, i = i, A = A), n_draws, seed)
}

#' Flow accumulation along a drainage network
#'
#' Accumulates a quantity (cell counts by default) downstream along single
#' flow directions in linear time, without building a dense matrix. Use the
#' `to` vector from [dm_flow_from_dem()].
#'
#' @param to Downstream index for each cell, `NA` for outlets.
#' @param weights Value contributed by each cell (default 1, giving the number
#'   of cells draining through each cell, itself included).
#' @param delivery Fraction passed downstream from each cell (default 1).
#' @return Numeric vector of accumulated values.
#' @examples
#' dm_flow_accumulation(c(3, 3, 4, NA))
#' @export
dm_flow_accumulation <- function(to, weights = 1, delivery = 1) {
  n <- length(to)
  to[is.na(to)] <- 0L
  to <- as.integer(to)
  if (any(to < 0 | to > n) || any(to == seq_len(n))) {
    stop("'to' must index other cells, with NA for outlets.", call. = FALSE)
  }
  w <- rep_len(as.numeric(weights), n)
  d <- rep_len(as.numeric(delivery), n)
  indeg <- tabulate(to[to > 0], nbins = n)
  acc <- w
  queue <- integer(n)                 # preallocated FIFO queue (linear time)
  start <- which(indeg == 0L)
  tail <- length(start)
  queue[seq_len(tail)] <- start
  head <- 0L
  while (head < tail) {
    head <- head + 1L
    u <- queue[head]
    v <- to[u]
    if (v > 0L) {
      acc[v] <- acc[v] + acc[u] * d[u]
      indeg[v] <- indeg[v] - 1L
      if (indeg[v] == 0L) { tail <- tail + 1L; queue[tail] <- v }
    }
  }
  done <- head
  if (done < n) stop("'to' contains a cycle.", call. = FALSE)
  acc
}

.warn_steep <- function(slope_deg, limit = 26.6) {
  sv <- slope_deg[is.finite(slope_deg)]
  if (length(sv) && mean(sv > limit) > 0.05) {
    warning(sprintf("%.0f%% of slopes exceed %.1f degrees (50%%). The (R)USLE slope factors were ",
                    100 * mean(sv > limit), limit),
            "derived from plots on gentler slopes; consider steepness = \"nearing\".",
            call. = FALSE)
  }
  invisible(NULL)
}

#' Slope length and steepness (LS) factor
#'
#' Computes the topographic LS factor of the (Revised) Universal Soil Loss
#' Equation.
#'
#' * `method = "moore"` (Moore and Burch, 1986), from flow accumulation:
#'   \eqn{LS = (A_s / 22.13)^{0.4} (\sin\theta / 0.0896)^{1.3}}, where
#'   \eqn{A_s} is the specific catchment area (upslope area per unit contour
#'   width, m), approximated by flow accumulation (cells) times cell size.
#' * `method = "wischmeier"` (Wischmeier and Smith, 1978), from slope length:
#'   \eqn{LS = (\lambda / 22.13)^m (65.41 \sin^2\theta + 4.56 \sin\theta + 0.065)},
#'   with \eqn{m} = 0.5 for slopes of 5\% or more, 0.4 for 3.5-4.5\%, 0.3 for
#'   1-3\% and 0.2 below 1\% (intermediate gaps take the nearer class).
#'
#' @param slope Slope angle in degrees.
#' @param flow_acc Flow accumulation in cells (method `"moore"`), e.g. from
#'   [dm_flow_accumulation()].
#' @param cell_size Cell size (m) (method `"moore"`).
#' @param length Slope length (m) (method `"wischmeier"`).
#' @param method `"moore"` or `"wischmeier"`.
#' @param steepness Slope-steepness relation: `"standard"` (McCool et al.
#'   1987 form in the chosen method) or `"nearing"`, the single continuous
#'   function of Nearing (1997), \eqn{S = -1.5 + 17 / (1 + e^{2.3 - 6.1 \sin
#'   \theta})}, fitted to data up to about 55\% slope and recommended for steep
#'   terrain such as the Himalaya and the Western Ghats.
#' @return LS factor (dimensionless). A warning is given when more than 5\%
#'   of slopes exceed 50\% (26.6 degrees), beyond the plot data on which the
#'   (R)USLE slope relations were based.
#' @references Moore, I. D. and Burch, G. J. (1986). Physical basis of the
#'   length-slope factor in the Universal Soil Loss Equation. *Soil Science
#'   Society of America Journal*, 50, 1294-1298.
#'
#'   Nearing, M. A. (1997). A single, continuous function for slope
#'   steepness influence on soil loss. *Soil Science Society of America
#'   Journal*, 61, 917-919.
#'
#'   Wischmeier, W. H. and Smith, D. D. (1978). *Predicting Rainfall Erosion
#'   Losses: A Guide to Conservation Planning*. Agriculture Handbook 537,
#'   USDA, Washington, DC.
#' @examples
#' dm_ls_factor(slope = 8, length = 50, method = "wischmeier")
#' dm_ls_factor(slope = 8, flow_acc = 20, cell_size = 30)
#' @export
dm_ls_factor <- function(slope, flow_acc = NULL, cell_size = NULL, length = NULL,
                         method = c("moore", "wischmeier"),
                         steepness = c("standard", "nearing")) {
  method <- match.arg(method); steepness <- match.arg(steepness)
  .check_range(slope, "slope", 0, 90)
  sv <- if (inherits(slope, "SpatRaster")) terra::values(slope[[1]], mat = FALSE) else
    if (.is_uncertain(slope)) slope$mean else slope
  if (steepness == "standard") .warn_steep(sv)
  if (method == "moore") {
    if (is.null(flow_acc) || is.null(cell_size)) {
      stop("method = 'moore' needs 'flow_acc' and 'cell_size'.", call. = FALSE)
    }
    f_moore <- function(slope, flow_acc) {
      As <- flow_acc * cell_size
      S <- if (steepness == "nearing") .s_nearing(slope) else
        (sin(slope * pi / 180) / 0.0896)^1.3
      (As / 22.13)^0.4 * S
    }
    return(.apply_formula(f_moore, list(slope = slope, flow_acc = flow_acc)))
  }
  if (is.null(length)) stop("method = 'wischmeier' needs 'length'.", call. = FALSE)
  f_ws <- function(slope, length) {
    th <- slope * pi / 180
    pct <- 100 * tan(th)
    m <- ifelse(pct >= 4.75, 0.5, ifelse(pct >= 3.25, 0.4, ifelse(pct >= 1, 0.3, 0.2)))
    S <- if (steepness == "nearing") .s_nearing(slope) else
      65.41 * sin(th)^2 + 4.56 * sin(th) + 0.065
    (length / 22.13)^m * S
  }
  .apply_formula(f_ws, list(slope = slope, length = length))
}

# Nearing (1997) continuous slope-steepness factor, valid on steep slopes
.s_nearing <- function(slope_deg) -1.5 + 17 / (1 + exp(2.3 - 6.1 * sin(slope_deg * pi / 180)))

#' Cover-management (C) factor from NDVI
#'
#' \eqn{C = \exp(-\alpha\, NDVI / (\beta - NDVI))}, bounded to \[0, 1\]
#' (van der Knijff et al., 1999).
#'
#' @param ndvi Normalised difference vegetation index.
#' @param alpha,beta Shape parameters (2 and 1 in van der Knijff et al.).
#' @return C factor.
#' @references van der Knijff, J. M., Jones, R. J. A. and Montanarella, L.
#'   (1999). *Soil Erosion Risk Assessment in Italy*. EUR 19022 EN, Joint
#'   Research Centre, European Commission.
#' @examples
#' dm_c_ndvi(c(0.1, 0.4, 0.8))
#' @export
dm_c_ndvi <- function(ndvi, alpha = 2, beta = 1) {
  .check_range(ndvi, "ndvi", -1, 1)
  f <- function(ndvi) {
    nd <- pmin(ndvi, beta - 1e-6)
    pmin(pmax(exp(-alpha * nd / (beta - nd)), 0), 1)
  }
  .apply_formula(f, list(ndvi = ndvi))
}

#' Soil erodibility (K) factor from texture and organic carbon (EPIC)
#'
#' Williams (1995) EPIC equation:
#' \deqn{K = f_{csand}\, f_{cl-si}\, f_{orgc}\, f_{hisand} \times 0.1317,}
#' \deqn{f_{csand} = 0.2 + 0.3 \exp[-0.256\, SAN (1 - SIL/100)],}
#' \deqn{f_{cl-si} = [SIL / (CLA + SIL)]^{0.3},}
#' \deqn{f_{orgc} = 1 - 0.25\, C / [C + \exp(3.72 - 2.95 C)],}
#' \deqn{f_{hisand} = 1 - 0.7\, SN_1 / [SN_1 + \exp(-5.51 + 22.9 SN_1)],
#'       \quad SN_1 = 1 - SAN/100.}
#' The factor 0.1317 converts from US customary units to
#' t ha h ha\eqn{^{-1}} MJ\eqn{^{-1}} mm\eqn{^{-1}}.
#'
#' @param sand,silt,clay Particle-size fractions (\%).
#' @param orgC Organic carbon (\%).
#' @param si Return SI units (`TRUE`) or US customary units.
#' @return K factor.
#' @references Williams, J. R. (1995). The EPIC model. In V. P. Singh (ed.),
#'   *Computer Models of Watershed Hydrology*, 909-1000. Water Resources
#'   Publications, Highlands Ranch, CO.
#' @examples
#' dm_k_epic(sand = 20, silt = 60, clay = 20, orgC = 1.5)
#' @export
dm_k_epic <- function(sand, silt, clay, orgC, si = TRUE) {
  .check_range(sand, "sand", 0, 100); .check_range(silt, "silt", 0, 100)
  .check_range(clay, "clay", 0, 100); .check_range(orgC, "orgC", 0, 100)
  f <- function(sand, silt, clay, orgC) {
    fcsand <- 0.2 + 0.3 * exp(-0.256 * sand * (1 - silt / 100))
    fclsi <- (silt / pmax(clay + silt, 1e-9))^0.3
    forgc <- 1 - 0.25 * orgC / (orgC + exp(3.72 - 2.95 * orgC))
    sn1 <- 1 - sand / 100
    fhisand <- 1 - 0.7 * sn1 / (sn1 + exp(-5.51 + 22.9 * sn1))
    fcsand * fclsi * forgc * fhisand * (if (si) 0.1317 else 1)
  }
  .apply_formula(f, list(sand = sand, silt = silt, clay = clay, orgC = orgC))
}

#' Modified Fournier index and rainfall erosivity
#'
#' Computes the modified Fournier index \eqn{MFI = \sum p_i^2 / P} from
#' monthly rainfall, and optionally rainfall erosivity by a regional linear
#' relation \eqn{R = a + b \cdot MFI}. Regional coefficients `a` and `b` must
#' be supplied by the user; there is no universal default.
#'
#' @param monthly Matrix or data frame of monthly rainfall (mm): one row per
#'   location, 12 columns.
#' @param a,b Optional coefficients of the regional erosivity relation.
#' @return MFI (if `a` and `b` are `NULL`) or erosivity R in the units of the
#'   regional relation.
#' @references Arnoldus, H. M. J. (1980). An approximation of the rainfall
#'   factor in the Universal Soil Loss Equation. In M. De Boodt and D.
#'   Gabriels (eds), *Assessment of Erosion*, 127-132. Wiley, Chichester.
#' @examples
#' rain <- rbind(c(10, 15, 20, 30, 60, 200, 320, 300, 180, 60, 15, 10))
#' dm_fournier(rain)
#' @export
dm_fournier <- function(monthly, a = NULL, b = NULL) {
  m <- as.matrix(monthly)
  if (ncol(m) != 12) stop("'monthly' must have 12 columns.", call. = FALSE)
  if (any(m < 0, na.rm = TRUE)) stop("Rainfall cannot be negative.", call. = FALSE)
  mfi <- rowSums(m^2) / rowSums(m)
  if (is.null(a) || is.null(b)) return(mfi)
  a + b * mfi
}

#' Soil loss by the Revised Universal Soil Loss Equation
#'
#' \eqn{A = R \cdot K \cdot LS \cdot C \cdot P}. Any factor can be a number,
#' a raster or an uncertain input ([dm_uncertain()]); with uncertain factors
#' the result is a predictive distribution of soil loss that can be used
#' directly with [dm_decide()], [dm_voi()] and [dm_design()].
#'
#' @param R Rainfall erosivity (MJ mm ha\eqn{^{-1}} h\eqn{^{-1}} yr\eqn{^{-1}}).
#' @param K Soil erodibility (t ha h ha\eqn{^{-1}} MJ\eqn{^{-1}} mm\eqn{^{-1}}).
#' @param LS Slope length and steepness factor.
#' @param C Cover-management factor.
#' @param P Support-practice factor (1 without conservation practices).
#' @param n_draws,seed Monte Carlo settings used when an input is uncertain.
#' @section Scope: (R)USLE estimates long-term average sheet and rill
#'   erosion. It does not represent gully, ravine, landslide, streambank or
#'   tillage erosion, or deposition within a field; use [dm_sediment_yield()]
#'   for delivery to streams.
#' @return Soil loss (t ha\eqn{^{-1}} yr\eqn{^{-1}}).
#' @references Renard, K. G., Foster, G. R., Weesies, G. A., McCool, D. K. and
#'   Yoder, D. C. (1997). *Predicting Soil Erosion by Water: A Guide to
#'   Conservation Planning with the Revised Universal Soil Loss Equation
#'   (RUSLE)*. Agriculture Handbook 703, USDA, Washington, DC.
#' @examples
#' dm_rusle(R = 800, K = 0.03, LS = 2.5, C = 0.3)
#' A <- dm_rusle(R = dm_uncertain(800, cv = 0.15, lower = 0), K = 0.03,
#'               LS = 2.5, C = dm_uncertain(0.3, sd = 0.05, lower = 0, upper = 1),
#'               seed = 1)
#' dm_decide(A, dm_loss_tolerance())
#' @export
dm_rusle <- function(R, K, LS, C, P = 1, n_draws = 500, seed = NULL) {
  for (nm in c("R", "K", "LS")) .check_range(get(nm), nm, 0)
  .check_range(C, "C", 0, 1); .check_range(P, "P", 0, 1)
  .apply_formula(function(R, K, LS, C, P) R * K * LS * C * P,
                 list(R = R, K = K, LS = LS, C = C, P = P), n_draws, seed)
}

#' Storm sediment yield by the Modified Universal Soil Loss Equation (MUSLE)
#'
#' \eqn{Y = a (Q \cdot q_p)^b K \cdot LS \cdot C \cdot P}, with \eqn{Y} in
#' tonnes per storm, runoff volume \eqn{Q} in m3 and peak discharge
#' \eqn{q_p} in m3/s (Williams, 1975: a = 11.8, b = 0.56).
#'
#' @param Q Runoff volume (m3).
#' @param qp Peak runoff rate (m3/s).
#' @param K,LS,C,P USLE factors.
#' @param a,b Coefficients.
#' @param n_draws,seed Monte Carlo settings used when an input is uncertain.
#' @return Sediment yield (t per storm event).
#' @references Williams, J. R. (1975). Sediment-yield prediction with
#'   universal equation using runoff energy factor. In *Present and
#'   Prospective Technology for Predicting Sediment Yields and Sources*,
#'   ARS-S-40, 244-252. USDA Agricultural Research Service.
#' @examples
#' dm_musle(Q = 1500, qp = 0.8, K = 0.03, LS = 2, C = 0.3, P = 1)
#' @export
dm_musle <- function(Q, qp, K, LS, C, P = 1, a = 11.8, b = 0.56, n_draws = 500, seed = NULL) {
  .check_range(Q, "Q", 0); .check_range(qp, "qp", 0)
  .apply_formula(function(Q, qp, K, LS, C, P) a * (Q * qp)^b * K * LS * C * P,
                 list(Q = Q, qp = qp, K = K, LS = LS, C = C, P = P), n_draws, seed)
}

#' Sediment yield from gross soil loss
#'
#' Sediment yield = gross soil loss x sediment delivery ratio (SDR).
#'
#' @param soil_loss Gross soil loss (e.g. from [dm_rusle()]).
#' @param sdr Sediment delivery ratio (0 to 1), a number, raster or
#'   uncertain input.
#' @param n_draws,seed Monte Carlo settings used when an input is uncertain.
#' @return Sediment yield, in the units of `soil_loss`.
#' @examples
#' dm_sediment_yield(20, sdr = 0.35)
#' @export
dm_sediment_yield <- function(soil_loss, sdr, n_draws = 500, seed = NULL) {
  .check_range(sdr, "sdr", 0, 1)
  .apply_formula(function(soil_loss, sdr) soil_loss * sdr,
                 list(soil_loss = soil_loss, sdr = sdr), n_draws, seed)
}

#' Nutrient loss attached to eroded sediment (mass balance)
#'
#' Sediment-bound nutrient loss (kg ha\eqn{^{-1}}) = sediment
#' (t ha\eqn{^{-1}}) x nutrient content of the soil x enrichment ratio.
#' With content in percent this is \eqn{10 \cdot S \cdot c \cdot ER}; with
#' content in mg kg\eqn{^{-1}} it is \eqn{0.001 \cdot S \cdot c \cdot ER}.
#'
#' @param sediment Sediment or soil loss (t ha\eqn{^{-1}}).
#' @param content Nutrient content of the source soil.
#' @param unit Unit of `content`: `"percent"` or `"mg_kg"`.
#' @param enrichment Enrichment ratio (nutrient content of sediment / of the
#'   source soil; usually above 1), a number, raster, uncertain input, or the
#'   result of [dm_enrichment()].
#' @param n_draws,seed Monte Carlo settings used when an input is uncertain.
#' @return Nutrient loss (kg ha\eqn{^{-1}}).
#' @examples
#' # 15 t/ha of sediment, 0.08 % total N, enrichment ratio 1.5
#' dm_nutrient_sediment(15, 0.08, "percent", enrichment = 1.5)
#' @export
dm_nutrient_sediment <- function(sediment, content, unit = c("percent", "mg_kg"),
                                 enrichment = 1, n_draws = 500, seed = NULL) {
  unit <- match.arg(unit)
  .check_range(sediment, "sediment", 0); .check_range(content, "content", 0)
  .check_range(enrichment, "enrichment", 0)
  k <- if (unit == "percent") 10 else 0.001
  .apply_formula(function(sediment, content, enrichment) k * sediment * content * enrichment,
                 list(sediment = sediment, content = content, enrichment = enrichment),
                 n_draws, seed)
}

#' Dissolved nutrient loss in runoff
#'
#' Dissolved loss (kg ha\eqn{^{-1}}) = runoff depth (mm) x concentration
#' (mg L\eqn{^{-1}}) x 0.01.
#'
#' @param runoff Runoff depth (mm), e.g. from [dm_runoff_cn()].
#' @param conc Nutrient concentration in runoff water (mg/L).
#' @param n_draws,seed Monte Carlo settings used when an input is uncertain.
#' @return Dissolved nutrient loss (kg ha\eqn{^{-1}}).
#' @examples
#' dm_nutrient_runoff(runoff = 40, conc = 2.5)
#' @export
dm_nutrient_runoff <- function(runoff, conc, n_draws = 500, seed = NULL) {
  .check_range(runoff, "runoff", 0); .check_range(conc, "conc", 0)
  .apply_formula(function(runoff, conc) 0.01 * runoff * conc,
                 list(runoff = runoff, conc = conc), n_draws, seed)
}

#' Enrichment ratio from sediment load (power law)
#'
#' Enrichment ratios usually fall as sediment load rises, and are commonly
#' described by a power law \eqn{ER = a \cdot S^b} with \eqn{b < 0} (Menzel,
#' 1980). Coefficients depend on the nutrient, soil and region, so they must
#' be supplied by the user. Results are bounded below by `min_er`.
#'
#' @param sediment Sediment load (in the units the coefficients were derived
#'   for).
#' @param a,b Power-law coefficients.
#' @param min_er Lower bound for the ratio (default 1).
#' @return Enrichment ratio.
#' @references Menzel, R. G. (1980). Enrichment ratios for water quality
#'   modeling. In W. G. Knisel (ed.), *CREAMS: A Field-Scale Model for
#'   Chemicals, Runoff, and Erosion from Agricultural Management Systems*,
#'   Conservation Research Report 26, 486-492. USDA.
#' @examples
#' dm_enrichment(c(100, 1000, 5000), a = 7.4, b = -0.2)
#' @export
dm_enrichment <- function(sediment, a, b, min_er = 1) {
  if (missing(a) || missing(b)) stop("Supply regional coefficients 'a' and 'b'.", call. = FALSE)
  .check_range(sediment, "sediment", 0)
  .apply_formula(function(sediment) pmax(a * pmax(sediment, 1e-12)^b, min_er),
                 list(sediment = sediment))
}

#' Soil-loss tolerance decision
#'
#' Ready-made [dm_loss()] for the question "is soil loss above the
#' tolerable limit, so that conservation treatment is needed?". The default
#' limit of 11.2 t ha\eqn{^{-1}} yr\eqn{^{-1}} is the upper tolerance commonly
#' used in conservation planning; supply site-specific limits where
#' available (for shallow or degraded soils they are lower).
#'
#' Tolerable limits differ between soils and regions. In India, permissible
#' limits of 2.5 to 12.5 t ha\eqn{^{-1}} yr\eqn{^{-1}} have been estimated
#' for soil mapping units from soil quality indicators (Lenka, Mandal and
#' Sudhishri, 2014, *Current Science* 107(4), for West Bengal). For a map of
#' site-specific limits, decide on the ratio of soil loss to the local limit
#' with `T_limit = 1` (see the example).
#'
#' @param T_limit Tolerable soil loss (t ha\eqn{^{-1}} yr\eqn{^{-1}}).
#' @param ratio Cost of leaving an eroding area untreated relative to the cost
#'   of treating an area that did not need it.
#' @param type `"fixed"` or `"linear"` loss (see [dm_loss()]).
#' @return A `dm_loss` object.
#' @examples
#' dm_loss_tolerance()
#' # site-specific limits: decide on soil loss / local limit
#' loss <- c(6, 9, 14); local_T <- c(5, 12.5, 10)
#' dm_decide(dm_gaussian(loss / local_T, sd = 0.2), dm_loss_tolerance(T_limit = 1))
#' @export
dm_loss_tolerance <- function(T_limit = 11.2, ratio = 3, type = c("fixed", "linear")) {
  dm_loss(T_limit, "above", ratio = ratio, type = match.arg(type),
          action = "conservation treatment")
}

#' Nutrient-loss limit decision
#'
#' Ready-made [dm_loss()] for "does nutrient loss exceed an acceptable
#' limit?".
#'
#' @param limit Acceptable nutrient loss (kg ha\eqn{^{-1}}).
#' @param nutrient Label used in outputs, e.g. `"N"`, `"P"`, `"K"`.
#' @param ratio,type See [dm_loss()].
#' @return A `dm_loss` object.
#' @examples
#' dm_loss_nutrient(10, "N")
#' @export
dm_loss_nutrient <- function(limit, nutrient = "nutrient", ratio = 2, type = c("fixed", "linear")) {
  dm_loss(limit, "above", ratio = ratio, type = match.arg(type),
          action = paste("reduce", nutrient, "loss"))
}

#' Load delivered to the watershed outlet(s)
#'
#' Routes a per-cell quantity (soil loss, sediment or nutrient loss, in mass
#' per cell) down the drainage network and returns the load at each outlet,
#' with uncertainty carried through when `x` is a predictive distribution.
#'
#' @param x Per-cell values: numeric, or a `dm_predictive`.
#' @param to Flow directions from [dm_flow_from_dem()].
#' @param delivery Fraction passed from each cell to the next (one value or
#'   one per cell).
#' @param cell_area Area of a cell (ha), to convert per-hectare rates to
#'   per-cell mass.
#' @return For numeric `x`, a named vector of outlet loads; for a
#'   `dm_predictive`, a `dm_predictive` of outlet loads.
#' @examples
#' elev <- outer(1:6, 1:6, "+")
#' to <- dm_flow_from_dem(elev)
#' dm_outlet_load(rep(10, 36), to, delivery = 0.9, cell_area = 0.09)
#' @export
dm_outlet_load <- function(x, to, delivery = 1, cell_area = 1) {
  outlets <- which(is.na(to))
  if (inherits(x, "dm_predictive")) {
    if (!is.null(x$draws)) {
      # route every draw in linear time; no dense matrix needed
      L <- vapply(seq_len(ncol(x$draws)), function(j) {
        dm_flow_accumulation(to, weights = x$draws[, j] * cell_area, delivery = delivery)[outlets]
      }, numeric(length(outlets)))
      L <- matrix(L, nrow = length(outlets))
      return(dm_draws(L))
    }
    if (length(to) > 5000) {
      stop("For more than 5000 cells, pass draws (dm_draws) so loads can be routed per draw.",
           call. = FALSE)
    }
    H <- dm_flow_matrix(to, delivery)[outlets, , drop = FALSE] * cell_area
    return(dm_aggregate(x, H))
  }
  acc <- dm_flow_accumulation(to, weights = as.numeric(x) * cell_area, delivery = delivery)
  stats::setNames(acc[outlets], paste0("outlet_", outlets))
}

#' Daily runoff by the curve number method with antecedent moisture
#'
#' Applies the curve number method to each day of a rainfall series, adjusting
#' the curve number for antecedent moisture condition (AMC) from the rainfall
#' of the previous five days (USDA SCS National Engineering Handbook, Section
#' 4). The curve number method is an event method; summing daily runoff is the
#' valid way to estimate seasonal or annual runoff, rather than applying the
#' method to annual rainfall.
#'
#' @param P Daily rainfall (mm): a numeric vector (one location) or a matrix
#'   with one row per location and one column per day.
#' @param CN Curve number for AMC II: one value or one per location. With a
#'   single rainfall series and several curve numbers, the series is applied
#'   to each.
#' @param growing Logical, one per day (or one value): whether the day falls
#'   in the growing season, which sets the AMC limits (growing season: below
#'   35.6 mm is AMC I, above 53.3 mm AMC III; dormant season: 12.7 and
#'   27.9 mm).
#' @param lambda Initial abstraction ratio.
#' @return A list with `daily` (runoff matrix, mm, locations by days),
#'   `total` (runoff per location, mm), `amc` (AMC class per location and day)
#'   and `coefficient` (total runoff / total rainfall).
#' @examples
#' set.seed(1)
#' rain <- ifelse(runif(120) < 0.3, rexp(120, 1 / 25), 0)
#' out <- dm_runoff_cn_daily(rain, CN = 78)
#' out$total; out$coefficient
#' table(out$amc)
#' @export
dm_runoff_cn_daily <- function(P, CN, growing = TRUE, lambda = 0.2) {
  P <- if (is.matrix(P)) P else matrix(P, nrow = 1)
  .check_range(P, "P", 0); .check_range(CN, "CN", 1, 100)
  if (nrow(P) == 1L && length(CN) > 1L) P <- P[rep(1L, length(CN)), , drop = FALSE]
  nl <- nrow(P); nd <- ncol(P)
  if (!length(CN) %in% c(1L, nl)) stop("'CN' must have length 1 or ", nl, ".", call. = FALSE)
  if (!length(growing) %in% c(1L, nd)) stop("'growing' must have length 1 or ", nd, ".", call. = FALSE)
  CN <- rep_len(CN, nl); growing <- rep_len(as.logical(growing), nd)
  cn1 <- dm_cn_amc(CN, "I"); cn3 <- pmin(dm_cn_amc(CN, "III"), 100)
  Q <- A <- matrix(NA, nl, nd)
  for (d in seq_len(nd)) {
    ant <- if (d == 1L) rep(0, nl) else rowSums(P[, max(1L, d - 5L):(d - 1L), drop = FALSE])
    lo <- if (growing[d]) 35.6 else 12.7; hi <- if (growing[d]) 53.3 else 27.9
    amc <- ifelse(ant < lo, 1L, ifelse(ant > hi, 3L, 2L))
    cn <- ifelse(amc == 1L, cn1, ifelse(amc == 3L, cn3, CN))
    S <- 25400 / cn - 254; Ia <- lambda * S
    Q[, d] <- ifelse(P[, d] > Ia, (P[, d] - Ia)^2 / (P[, d] - Ia + S), 0)
    A[, d] <- amc
  }
  tot <- rowSums(Q)
  list(daily = Q, total = tot, amc = factor(A, 1:3, c("I", "II", "III")),
       coefficient = tot / pmax(rowSums(P), .Machine$double.eps))
}


#' Terrain indices for gully and ravine erosion risk
#'
#' (R)USLE does not predict gully or ravine erosion. Where concentrated flow
#' is strong enough, gullies form; this is screened with the stream power
#' index \eqn{SPI = A_s \tan\beta}, the topographic wetness index
#' \eqn{TWI = \ln(A_s / \tan\beta)} (Moore et al., 1991) and the topographic
#' threshold \eqn{\tan\beta \, A^b > k} for gully initiation (Patton and
#' Schumm, 1975; Vandaele et al., 1996).
#'
#' @param slope Slope (degrees).
#' @param area Specific contributing area (m2 per m of contour width), for
#'   example flow accumulation times cell size.
#' @param threshold Gully-initiation threshold \eqn{k}. Thresholds are
#'   site-specific: calibrate it by mapping gully heads on imagery and
#'   taking a low quantile of \eqn{\tan\beta \, A^b} at those heads. If
#'   `NULL`, no cells are flagged.
#' @param b Exponent of contributing area in the threshold (commonly 0.2 to
#'   0.5).
#' @return A data frame with `spi`, `twi`, `topo` (\eqn{\tan\beta\,A^b}) and,
#'   when `threshold` is given, logical `at_risk`.
#' @references Moore, I. D., Grayson, R. B. and Ladson, A. R. (1991). Digital
#'   terrain modelling: a review of hydrological, geomorphological, and
#'   biological applications. *Hydrological Processes*, 5, 3-30.
#'
#'   Patton, P. C. and Schumm, S. A. (1975). Gully erosion, northwestern
#'   Colorado: a threshold phenomenon. *Geology*, 3, 88-90.
#'
#'   Vandaele, K., Poesen, J., Govers, G. and van Wesemael, B. (1996).
#'   Geomorphic threshold conditions for ephemeral gully incision.
#'   *Geomorphology*, 16, 161-173.
#' @examples
#' dm_gully_risk(slope = c(2, 8, 15), area = c(50, 400, 2000), threshold = 5)
#' @export
dm_gully_risk <- function(slope, area, threshold = NULL, b = 0.4) {
  .check_range(slope, "slope", 0, 90); .check_range(area, "area", 0)
  n <- max(length(slope), length(area))
  if (!all(c(length(slope), length(area)) %in% c(1L, n))) {
    stop("'slope' and 'area' must have length 1 or the same length.", call. = FALSE)
  }
  tb <- pmax(tan(slope * pi / 180), 1e-6)
  out <- data.frame(spi = area * tb, twi = log(pmax(area, 1e-6) / tb), topo = tb * area^b)
  if (!is.null(threshold)) out$at_risk <- out$topo > threshold
  out
}
