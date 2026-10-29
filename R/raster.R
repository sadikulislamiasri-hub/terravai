# ---------------------------------------------------------------------------
# Raster output (optional, via terra)
# ---------------------------------------------------------------------------

#' Convert decision or value-of-information results to a raster
#'
#' Writes chosen columns of a `dm_decision` or `dm_voi` result into a
#' `SpatRaster` for use in GIS software. If the predictions were created from
#' rasters, the original grid is reused; otherwise a grid is built from
#' regularly spaced coordinates.
#'
#' @param x A `dm_decision` or `dm_voi` object.
#' @param layers Column names to write. Defaults to the main outputs.
#' @param crs Coordinate reference system for a newly built grid.
#' @return A `SpatRaster` (package `terra`).
#' @examples
#' \donttest{
#' if (requireNamespace("terra", quietly = TRUE)) {
#'   sim <- dm_simulate(n_side = 10, n_obs = 20, seed = 8)
#'   fit <- dm_krige(sim$obs[, c("x", "y")], sim$obs$z, sim$grid[, c("x", "y")],
#'                   model = sim$model)
#'   r <- dm_as_raster(dm_decide(fit, dm_loss(0, "above")))
#'   r
#' }
#' }
#' @export
dm_as_raster <- function(x, layers = NULL, crs = "") {
  if (!requireNamespace("terra", quietly = TRUE)) {
    stop("Package 'terra' is needed for raster output.", call. = FALSE)
  }
  if (!inherits(x, c("dm_decision", "dm_voi"))) {
    stop("'x' must be a dm_decision or dm_voi object.", call. = FALSE)
  }
  layers <- layers %||% if (inherits(x, "dm_voi")) "voi" else
    c("act", "p_need", "p_wrong", "exp_loss", "breakeven")
  df <- as.data.frame(unclass(x))
  miss <- setdiff(c("x", "y", layers), names(df))
  if (length(miss)) stop("Missing column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  info <- attr(x, "pred")
  if (!is.null(info$template) && inherits(x, "dm_decision")) {
    out <- terra::rast(info$template, nlyrs = length(layers))
    for (i in seq_along(layers)) {
      vals <- rep(NA_real_, terra::ncell(out))
      vals[info$cells] <- as.numeric(df[[layers[i]]])
      terra::values(out[[i]]) <- vals
    }
    names(out) <- layers
    return(out)
  }
  xyz <- data.frame(x = df$x, y = df$y, lapply(df[layers], as.numeric))
  terra::rast(xyz, type = "xyz", crs = crs)
}
