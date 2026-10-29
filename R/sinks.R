# ---------------------------------------------------------------------------
# Depression (sink) filling for real elevation models
# ---------------------------------------------------------------------------

#' Fill sinks in an elevation model
#'
#' Real elevation models (SRTM, Copernicus, CartoDEM) contain many small
#' artificial pits and flat areas. Without filling, single-direction flow
#' routing stops at every pit, so drainage areas, LS factors and outlet loads
#' are badly underestimated. This function fills depressions so that every
#' cell drains to the edge of the grid (or to a missing-value area), with a
#' small gradient `epsilon` imposed across filled areas so that flats also
#' drain. It uses the priority-flood algorithm (Barnes et al., 2014) in
#' compiled code, which gives the same surface as Planchon and Darboux (2002)
#' in \eqn{O(n \log n)} time: about a million cells per second.
#'
#' @param elev Numeric matrix or single-layer `SpatRaster` of elevations.
#' @param epsilon Minimum elevation rise imposed from the outlet side across
#'   filled areas (in elevation units). Use 0 for a pure fill (flats remain).
#' @return The filled elevations, in the form of the input.
#' @references Barnes, R., Lehman, C. and Mulla, D. (2014). Priority-flood:
#'   an optimal depression-filling and watershed-labeling algorithm for
#'   digital elevation models. *Computers & Geosciences*, 62, 117-127.
#'   \doi{10.1016/j.cageo.2013.04.024}
#'
#'   Planchon, O. and Darboux, F. (2002). A fast, simple and versatile
#'   algorithm to fill the depressions of digital elevation models.
#'   *Catena*, 46, 159-176. \doi{10.1016/S0341-8162(01)00164-3}
#' @examples
#' z <- outer(1:20, 1:20, "+"); z[10, 10] <- 0       # an artificial pit
#' sum(is.na(dm_flow_from_dem(z)))                    # 2 sinks
#' sum(is.na(dm_flow_from_dem(dm_fill_sinks(z))))     # 1 outlet
#' @export
dm_fill_sinks <- function(elev, epsilon = 1e-4) {
  if (inherits(elev, "SpatRaster")) {
    if (!requireNamespace("terra", quietly = TRUE)) stop("Package 'terra' is needed.", call. = FALSE)
    m <- terra::as.matrix(elev[[1]], wide = TRUE)
    f <- dm_fill_sinks(m, epsilon)
    out <- terra::rast(elev[[1]])
    terra::values(out) <- as.numeric(t(f))
    return(out)
  }
  Z <- as.matrix(elev)
  if (!is.numeric(Z)) stop("'elev' must be numeric.", call. = FALSE)
  if (!(is.numeric(epsilon) && length(epsilon) == 1 && epsilon >= 0)) {
    stop("'epsilon' must be a single non-negative number.", call. = FALSE)
  }
  storage.mode(Z) <- "double"
  out <- .Call(C_fill_sinks, Z, nrow(Z), ncol(Z), as.double(epsilon))
  dim(out) <- dim(Z)
  out
}
