# ---------------------------------------------------------------------------
# Simulated micro-watershed with known runoff, soil loss and nutrient loss
# ---------------------------------------------------------------------------

.smooth_field <- function(grid, range, n) {
  C <- dm_cov(grid[, c("x", "y")], model = dm_covmodel("exponential", 1, range))
  L <- chol(C + diag(1e-8, n))
  as.numeric(crossprod(L, stats::rnorm(n)))
}

#' Simulate a micro-watershed with known runoff, soil and nutrient loss
#'
#' Creates a synthetic micro-watershed on a regular grid (cell size in metres)
#' with a valley draining to one outlet, and the layers needed by the formula
#' functions: elevation, slope, flow directions and accumulation, LS factor,
#' NDVI and C factor, soil texture and K factor, rainfall, erosivity and curve
#' numbers. Annual runoff is approximated as 25 equal storms of the
#' curve-number event model. "True" annual runoff, soil loss and nitrogen, phosphorus and
#' potassium losses are generated from the formulas plus spatially correlated
#' deviations, so that the truth differs from any single model, as in real
#' landscapes. A set of runoff plots with noisy measurements is included.
#'
#' The erosivity relation used here (R = 0.5 x annual rainfall) is a
#' convenience for simulation only, not a recommended regional equation.
#'
#' @param n_side Grid cells along each side.
#' @param cell_size Cell size (m).
#' @param n_plots Number of runoff plots with measurements.
#' @param seed Optional random seed; the global random state is restored.
#' @return A list of class `dm_watershed` with `grid` (one row per cell),
#'   `to` (flow directions), `outlet` (outlet cell index), `plots` (measured
#'   soil loss, runoff and nutrient losses at plot cells) and `cell_size`.
#' @examples
#' ws <- dm_simulate_watershed(n_side = 20, seed = 1)
#' ws
#' @export
dm_simulate_watershed <- function(n_side = 30, cell_size = 30, n_plots = 40, seed = NULL) {
  if (n_side < 5) stop("'n_side' must be at least 5.", call. = FALSE)
  N <- n_side^2
  if (n_plots < 1 || n_plots > N) stop("'n_plots' must be between 1 and n_side^2.", call. = FALSE)
  .with_seed(seed, {
    g <- (seq_len(n_side) - 0.5) / n_side
    grid <- expand.grid(x = g, y = g)
    # valley along x = 0.5 sloping to the south edge; outlet at the south centre
    z <- 60 * grid$y + 90 * abs(grid$x - 0.5)^1.3 + 1.0 * .smooth_field(grid, 0.3, N)
    elev <- matrix(z, n_side, n_side)            # rows = x index, cols = y index
    to <- dm_flow_from_dem(elev)
    acc <- dm_flow_accumulation(to)
    # slope from central differences (degrees)
    h <- cell_size
    zx <- (elev[c(2:n_side, n_side), ] - elev[c(1, 1:(n_side - 1)), ]) /
      (h * c(1, rep(2, n_side - 2), 1))
    zy <- t((t(elev)[c(2:n_side, n_side), ] - t(elev)[c(1, 1:(n_side - 1)), ]) /
              (h * c(1, rep(2, n_side - 2), 1)))
    slope <- as.numeric(atan(sqrt(zx^2 + zy^2)) * 180 / pi)
    ls <- dm_ls_factor(pmax(slope, 0.1), flow_acc = acc, cell_size = cell_size)
    ndvi <- pmin(pmax(0.45 + 0.15 * .smooth_field(grid, 0.2, N) - 0.1 * (grid$y - 0.5), 0.05), 0.85)
    Cf <- dm_c_ndvi(ndvi)
    sand <- pmin(pmax(35 + 8 * .smooth_field(grid, 0.25, N), 5), 80)
    clay <- pmin(pmax(22 + 5 * .smooth_field(grid, 0.25, N), 5), 50)
    silt <- pmax(100 - sand - clay, 5)
    orgC <- pmin(pmax(0.9 + 0.25 * .smooth_field(grid, 0.2, N), 0.2), 3)
    K <- dm_k_epic(sand, silt, clay, orgC)
    rain <- 1400 + 120 * grid$x + 40 * .smooth_field(grid, 0.4, N)
    R <- 0.5 * rain
    cn <- pmin(pmax(78 - 25 * (ndvi - 0.45) + 3 * .smooth_field(grid, 0.2, N), 40), 95)
    # truth: formula x spatially correlated multiplicative deviation
    dev <- exp(0.35 * .smooth_field(grid, 0.12, N) - 0.35^2 / 2)
    soil_loss <- dm_rusle(R, K, ls, Cf) * dev
    # annual runoff as 25 equal storms of the curve-number event model
    runoff_mm <- 25 * dm_runoff_cn(rain / 25, cn) * exp(0.15 * .smooth_field(grid, 0.15, N))
    runoff_pct <- pmin(100 * runoff_mm / rain, 100)
    N_loss <- dm_nutrient_sediment(soil_loss, 0.08, "percent", enrichment = 1.6) +
      dm_nutrient_runoff(runoff_mm, 1.5)
    P_loss <- dm_nutrient_sediment(soil_loss, 0.03, "percent", enrichment = 2.0) +
      dm_nutrient_runoff(runoff_mm, 0.15)
    K_loss <- dm_nutrient_sediment(soil_loss, 0.6, "percent", enrichment = 1.2) +
      dm_nutrient_runoff(runoff_mm, 2.0)
    grid <- cbind(grid, elev = z, slope = slope, flow_acc = acc, ls = ls, ndvi = ndvi, C = Cf,
                  sand = sand, silt = silt, clay = clay, orgC = orgC, K = K, rain = rain,
                  R = R, cn = cn, runoff_mm = runoff_mm, runoff_pct = runoff_pct,
                  soil_loss = soil_loss, N_loss = N_loss, P_loss = P_loss, K_loss = K_loss)
    cells <- sort(sample.int(N, n_plots))
    noise <- function(v, cv) v * exp(stats::rnorm(length(v), 0, cv))
    plots <- data.frame(cell = cells, x = grid$x[cells], y = grid$y[cells],
                        soil_loss = noise(soil_loss[cells], 0.15),
                        runoff_pct = noise(runoff_pct[cells], 0.1),
                        N_loss = noise(N_loss[cells], 0.15),
                        P_loss = noise(P_loss[cells], 0.15),
                        K_loss = noise(K_loss[cells], 0.15))
    outlet <- which(is.na(to))
    structure(list(grid = grid, to = to, outlet = outlet, plots = plots,
                   cell_size = cell_size), class = "dm_watershed")
  })
}

#' @export
print.dm_watershed <- function(x, ...) {
  g <- x$grid
  ha <- nrow(g) * x$cell_size^2 / 1e4
  cat("<dm_watershed> ", nrow(g), " cells of ", x$cell_size, " m (", format(round(ha, 1)),
      " ha), ", length(x$outlet), " outlet(s), ", nrow(x$plots), " runoff plots\n", sep = "")
  q <- function(v) paste0(format(signif(stats::median(v), 3)), " (",
                          format(signif(min(v), 3)), "-", format(signif(max(v), 3)), ")")
  cat("  soil loss (t/ha/yr), median (range): ", q(g$soil_loss), "\n", sep = "")
  cat("  runoff (% of rainfall)             : ", q(g$runoff_pct), "\n", sep = "")
  cat("  N / P / K loss (kg/ha), median     : ",
      paste(format(signif(c(stats::median(g$N_loss), stats::median(g$P_loss),
                            stats::median(g$K_loss)), 3)), collapse = " / "), "\n", sep = "")
  cat("  share above 11.2 t/ha/yr           : ",
      sprintf("%.0f%%", 100 * mean(g$soil_loss > 11.2)), "\n", sep = "")
  invisible(x)
}
