# ---------------------------------------------------------------------------
# Global look-up tables and a raster pipeline for use anywhere
# ---------------------------------------------------------------------------

#' Runoff curve numbers (USDA TR-55)
#'
#' Curve numbers for average antecedent moisture (AMC II, \eqn{I_a = 0.2S})
#' by cover type, treatment and hydrologic condition, for hydrologic soil
#' groups A to D, from USDA-NRCS Technical Release 55. Use with
#' [dm_cn_lookup()]; replace or extend with local tables where available.
#'
#' @return A data frame with columns `cover`, `A`, `B`, `C`, `D`.
#' @references USDA-NRCS (1986). *Urban Hydrology for Small Watersheds*.
#'   Technical Release 55, 2nd edn. Natural Resources Conservation Service,
#'   Washington, DC.
#' @examples
#' dm_cn_table()
#' @export
dm_cn_table <- function() {
  data.frame(
    cover = c("fallow_bare", "row_crops_straight_poor", "row_crops_straight_good",
              "row_crops_contoured_poor", "row_crops_contoured_good",
              "contoured_terraced_poor", "contoured_terraced_good",
              "small_grain_straight_poor", "small_grain_straight_good",
              "pasture_poor", "pasture_fair", "pasture_good", "meadow",
              "woods_grass_poor", "woods_grass_fair", "woods_grass_good",
              "woods_poor", "woods_fair", "woods_good", "farmsteads", "impervious"),
    A = c(77, 72, 67, 70, 65, 66, 62, 65, 63, 68, 49, 39, 30, 57, 43, 32, 45, 36, 30, 59, 98),
    B = c(86, 81, 78, 79, 75, 74, 71, 76, 75, 79, 69, 61, 58, 73, 65, 58, 66, 60, 55, 74, 98),
    C = c(91, 88, 85, 84, 82, 80, 78, 84, 83, 86, 79, 74, 71, 82, 76, 72, 77, 73, 70, 82, 98),
    D = c(94, 91, 89, 88, 86, 82, 81, 88, 87, 89, 84, 80, 78, 86, 82, 79, 83, 79, 77, 86, 98),
    stringsAsFactors = FALSE)
}

#' Look up curve numbers by cover and hydrologic soil group
#'
#' @param cover Cover names as in [dm_cn_table()] (vector, one per location).
#' @param hsg Hydrologic soil group(s): `"A"`, `"B"`, `"C"` or `"D"` (for
#'   example from [dm_hsg_texture()]).
#' @param table A curve-number table with the same columns as
#'   [dm_cn_table()], e.g. a regional table.
#' @return Numeric vector of curve numbers.
#' @examples
#' dm_cn_lookup(c("row_crops_contoured_good", "woods_fair"), c("B", "C"))
#' @export
dm_cn_lookup <- function(cover, hsg, table = dm_cn_table()) {
  n <- max(length(cover), length(hsg))
  cover <- rep_len(as.character(cover), n); hsg <- rep_len(toupper(as.character(hsg)), n)
  bad <- setdiff(stats::na.omit(unique(cover)), table$cover)
  if (length(bad)) stop("Unknown cover: ", paste(bad, collapse = ", "),
                        ". See dm_cn_table().", call. = FALSE)
  if (any(!stats::na.omit(hsg) %in% c("A", "B", "C", "D"))) {
    stop("'hsg' must be A, B, C or D.", call. = FALSE)
  }
  out <- rep(NA_real_, n)
  ok <- !is.na(cover) & !is.na(hsg)
  out[ok] <- as.matrix(table[, c("A", "B", "C", "D")])[
    cbind(match(cover[ok], table$cover), match(hsg[ok], c("A", "B", "C", "D")))]
  out
}

#' USDA soil texture class
#'
#' Classifies soils into the 12 USDA texture classes from sand, silt and clay
#' percentages (which should sum to about 100).
#'
#' @param sand,silt,clay Particle-size fractions (\%).
#' @return Character vector of texture classes.
#' @references Soil Science Division Staff (2017). *Soil Survey Manual*.
#'   USDA Handbook 18. Government Printing Office, Washington, DC.
#' @examples
#' dm_texture_class(c(90, 40, 20, 10), c(5, 40, 60, 30), c(5, 20, 20, 60))
#' @export
dm_texture_class <- function(sand, silt, clay) {
  n <- max(length(sand), length(silt), length(clay))
  sa <- rep_len(sand, n); si <- rep_len(silt, n); cl <- rep_len(clay, n)
  if (any(c(sa, si, cl) < 0, na.rm = TRUE) || any(c(sa, si, cl) > 100, na.rm = TRUE)) {
    stop("Fractions must lie between 0 and 100.", call. = FALSE)
  }
  tot <- sa + si + cl
  if (any(abs(tot - 100) > 5, na.rm = TRUE)) warning("Some fractions do not sum to about 100.",
                                                      call. = FALSE)
  out <- rep(NA_character_, n)
  set <- function(cond, cls) { i <- is.na(out) & cond & !is.na(cond); out[i] <<- cls }
  set(si + 1.5 * cl < 15, "sand")
  set(si + 1.5 * cl >= 15 & si + 2 * cl < 30, "loamy sand")
  set((cl >= 7 & cl < 20 & sa > 52 & si + 2 * cl >= 30) | (cl < 7 & si < 50 & si + 2 * cl >= 30),
      "sandy loam")
  set(cl >= 7 & cl < 27 & si >= 28 & si < 50 & sa <= 52, "loam")
  set((si >= 50 & cl >= 12 & cl < 27) | (si >= 50 & si < 80 & cl < 12), "silt loam")
  set(si >= 80 & cl < 12, "silt")
  set(cl >= 20 & cl < 35 & si < 28 & sa > 45, "sandy clay loam")
  set(cl >= 27 & cl < 40 & sa > 20 & sa <= 45, "clay loam")
  set(cl >= 27 & cl < 40 & sa <= 20, "silty clay loam")
  set(cl >= 35 & sa > 45, "sandy clay")
  set(cl >= 40 & si >= 40, "silty clay")
  set(cl >= 40 & sa <= 45 & si < 40, "clay")
  out
}

#' Approximate hydrologic soil group from texture
#'
#' Assigns hydrologic soil groups from USDA texture classes when measured
#' infiltration or drainage data are not available: sand, loamy sand and
#' sandy loam to A; loam, silt loam and silt to B; sandy clay loam to C; and
#' clay loam, silty clay loam, sandy clay, silty clay and clay to D. Shallow
#' soils, high water tables or impeding layers should be moved to group D
#' using local knowledge.
#'
#' @param sand,silt,clay Particle-size fractions (\%), or supply `texture`.
#' @param texture Texture classes from [dm_texture_class()].
#' @return Character vector of hydrologic soil groups.
#' @examples
#' dm_hsg_texture(sand = c(85, 40, 30), silt = c(10, 40, 30), clay = c(5, 20, 40))
#' @export
dm_hsg_texture <- function(sand = NULL, silt = NULL, clay = NULL,
                           texture = dm_texture_class(sand, silt, clay)) {
  map <- c("sand" = "A", "loamy sand" = "A", "sandy loam" = "A",
           "loam" = "B", "silt loam" = "B", "silt" = "B",
           "sandy clay loam" = "C",
           "clay loam" = "D", "silty clay loam" = "D", "sandy clay" = "D",
           "silty clay" = "D", "clay" = "D")
  unname(map[texture])
}

#' Support-practice (P) factor by slope
#'
#' P factors for contouring and contour strip cropping by land slope, from
#' Wischmeier and Smith (1978). Slopes outside the tabulated range (1-25\%)
#' return `NA`. Terracing and other practices vary widely and should be
#' supplied from local tables.
#'
#' @param practice `"none"`, `"contouring"` or `"strip_cropping"` (one value
#'   or one per location).
#' @param slope_pct Land slope (\%).
#' @return Numeric vector of P factors.
#' @references Wischmeier, W. H. and Smith, D. D. (1978). *Predicting Rainfall
#'   Erosion Losses: A Guide to Conservation Planning*. Agriculture Handbook
#'   537, USDA, Washington, DC.
#' @examples
#' dm_p_factor("contouring", c(3, 10, 15, 22))
#' @export
dm_p_factor <- function(practice = c("none", "contouring", "strip_cropping"), slope_pct) {
  n <- max(length(practice), length(slope_pct))
  pr <- rep_len(practice, n); s <- rep_len(slope_pct, n)
  if (any(!pr %in% c("none", "contouring", "strip_cropping"))) {
    stop("'practice' must be 'none', 'contouring' or 'strip_cropping'.", call. = FALSE)
  }
  cont <- ifelse(s < 1 | s > 25, NA_real_,
          ifelse(s <= 8, 0.5, ifelse(s <= 12, 0.6, ifelse(s <= 16, 0.7,
          ifelse(s <= 20, 0.8, 0.9)))))
  strip <- ifelse(s < 1 | s > 25, NA_real_,
           ifelse(s <= 8, 0.25, ifelse(s <= 16, 0.30, 0.40)))
  ifelse(pr == "none", 1, ifelse(pr == "contouring", cont, strip))
}

#' Soil loss maps from global or national raster data
#'
#' Runs the RUSLE chain on rasters for any area in the world: slope and D8
#' flow routing from a digital elevation model, the LS factor (Moore and
#' Burch), the C factor from NDVI (or a supplied C raster), the K factor from
#' soil texture and organic carbon (EPIC, or a supplied K raster), and soil
#' loss \eqn{A = R K LS C P}. Geographic (longitude/latitude) grids are
#' handled by using each cell's true size in metres. Typical inputs are a
#' DEM (e.g. SRTM or Copernicus), NDVI (e.g. Sentinel-2 or MODIS), soil
#' texture and organic carbon (e.g. SoilGrids or national soil maps) and
#' rainfall erosivity (a national map or a global erosivity dataset).
#'
#' @param dem Elevation raster (m).
#' @param R Erosivity: a raster or a single value
#'   (MJ mm ha\eqn{^{-1}} h\eqn{^{-1}} yr\eqn{^{-1}}).
#' @param ndvi,C NDVI raster (to derive C) or a C raster/value.
#' @param sand,silt,clay,soc Texture and organic carbon rasters (to derive K),
#'   or supply `K`.
#' @param K Erodibility raster or value.
#' @param P Support-practice raster or value (default 1).
#' @param texture_unit Unit of the texture rasters: `"percent"` or `"g_kg"`
#'   (as in SoilGrids).
#' @param soc_unit Unit of organic carbon: `"percent"`, `"g_kg"` or `"dg_kg"`
#'   (as in SoilGrids).
#' @param max_slope_length Cap on the specific catchment area used in the LS
#'   factor (m), to avoid treating channels as hillslopes.
#' @param steepness Slope-steepness relation for the LS factor; see
#'   [dm_ls_factor()]. `"nearing"` is recommended for steep terrain.
#' @param fill_sinks Fill pits and flats before flow routing
#'   ([dm_fill_sinks()]); strongly recommended for real elevation models.
#' @return A `SpatRaster` with layers `slope` (degrees), `flow_acc` (cells),
#'   `LS`, `C`, `K`, `P`, `A` (soil loss, t ha\eqn{^{-1}} yr\eqn{^{-1}}),
#'   and the gully-risk indices `SPI` and `TWI` (see [dm_gully_risk()]),
#'   with the D8 flow directions attached as attribute `to` (for
#'   [dm_outlet_load()]).
#' @examples
#' if (requireNamespace("terra", quietly = TRUE)) {
#'   dem <- terra::rast(nrows = 30, ncols = 30, xmin = 77, xmax = 77.03,
#'                      ymin = 30, ymax = 30.03, crs = "EPSG:4326")
#'   terra::values(dem) <- outer(1:30, 1:30, function(i, j) 600 + 4 * i + 2 * j)
#'   ndvi <- terra::init(dem, 0.4)
#'   out <- dm_erosion_raster(dem, R = 900, ndvi = ndvi, K = 0.03)
#'   out
#' }
#' @export
dm_erosion_raster <- function(dem, R, ndvi = NULL, C = NULL, sand = NULL, silt = NULL,
                              clay = NULL, soc = NULL, K = NULL, P = 1,
                              texture_unit = c("percent", "g_kg"),
                              soc_unit = c("percent", "g_kg", "dg_kg"),
                              max_slope_length = 300, fill_sinks = TRUE,
                              steepness = c("standard", "nearing")) {
  if (!requireNamespace("terra", quietly = TRUE)) stop("Package 'terra' is needed.", call. = FALSE)
  if (!inherits(dem, "SpatRaster")) stop("'dem' must be a SpatRaster.", call. = FALSE)
  texture_unit <- match.arg(texture_unit); soc_unit <- match.arg(soc_unit)
  steepness <- match.arg(steepness)
  dem <- dem[[1]]
  align <- function(r, method = "bilinear") {
    if (is.null(r) || !inherits(r, "SpatRaster")) return(r)
    r <- r[[1]]
    if (!terra::compareGeom(r, dem, stopOnError = FALSE)) {
      if (terra::crs(r) != terra::crs(dem)) r <- terra::project(r, dem, method = method)
      else r <- terra::resample(r, dem, method = method)
    }
    terra::values(r, mat = FALSE)
  }
  val <- function(r, nm) {
    if (is.null(r)) return(NULL)
    if (inherits(r, "SpatRaster")) return(align(r))
    if (length(r) == 1L) return(rep(r, terra::ncell(dem)))
    stop("'", nm, "' must be a raster or a single value.", call. = FALSE)
  }
  z <- terra::values(dem, mat = FALSE)
  slope <- terra::values(terra::terrain(dem, v = "slope", unit = "degrees"), mat = FALSE)
  slope[is.na(slope) & !is.na(z)] <- 0
  to <- dm_flow_from_dem(if (fill_sinks) dm_fill_sinks(dem) else dem)
  acc <- dm_flow_accumulation(to)
  acc[is.na(z)] <- NA
  cell_m <- if (terra::is.lonlat(dem)) sqrt(terra::values(terra::cellSize(dem, unit = "m"),
                                                          mat = FALSE)) else
    rep(mean(terra::res(dem)), terra::ncell(dem))
  if (steepness == "standard") .warn_steep(slope)
  As_full <- acc * cell_m
  As <- pmin(As_full, max_slope_length)
  S <- if (steepness == "nearing") .s_nearing(pmax(slope, 0)) else
    (sin(pmax(slope, 0) * pi / 180) / 0.0896)^1.3
  LS <- (As / 22.13)^0.4 * S
  gi <- dm_gully_risk(pmax(slope, 0), ifelse(is.na(As_full), 0, As_full))
  Cv <- if (!is.null(C)) val(C, "C") else if (!is.null(ndvi)) {
    nd <- val(ndvi, "ndvi"); dm_c_ndvi(pmin(pmax(nd, -1), 1))
  } else stop("Supply 'ndvi' or 'C'.", call. = FALSE)
  Kv <- if (!is.null(K)) val(K, "K") else {
    if (any(vapply(list(sand, silt, clay, soc), is.null, logical(1)))) {
      stop("Supply 'K', or all of 'sand', 'silt', 'clay' and 'soc'.", call. = FALSE)
    }
    tf <- if (texture_unit == "g_kg") 0.1 else 1
    sf <- switch(soc_unit, percent = 1, g_kg = 0.1, dg_kg = 0.01)
    sa <- val(sand, "sand") * tf; si <- val(silt, "silt") * tf; cl <- val(clay, "clay") * tf
    oc <- val(soc, "soc") * sf
    ok <- stats::complete.cases(sa, si, cl, oc)
    kk <- rep(NA_real_, length(sa))
    kk[ok] <- dm_k_epic(pmin(pmax(sa[ok], 0), 100), pmin(pmax(si[ok], 0), 100),
                        pmin(pmax(cl[ok], 0), 100), pmax(oc[ok], 0))
    kk
  }
  Rv <- val(R, "R"); Pv <- val(P, "P")
  A <- Rv * Kv * LS * Cv * Pv
  mk <- function(v, nm) { r <- terra::rast(dem); terra::values(r) <- v; names(r) <- nm; r }
  out <- c(mk(slope, "slope"), mk(acc, "flow_acc"), mk(LS, "LS"), mk(Cv, "C"), mk(Kv, "K"),
           mk(Pv, "P"), mk(A, "A"), mk(ifelse(is.na(z), NA, gi$spi), "SPI"),
           mk(ifelse(is.na(z), NA, gi$twi), "TWI"))
  attr(out, "to") <- to
  out
}
