# ---------------------------------------------------------------------------
# Depression (sink) filling for real elevation models
# ---------------------------------------------------------------------------

#' Fill sinks in an elevation model
#'
#' Real elevation models (SRTM, Copernicus, CartoDEM) contain many small
#' artificial pits and flat areas. Without filling, single-direction flow
#' routing stops at every pit, so drainage areas, LS factors and outlet loads
#' are badly underestimated. This function fills depressions so that every
#' cell drains to the edge of the grid (or to a missing-value area), using the
#' algorithm of Planchon and Darboux (2002) with a small gradient
#' `epsilon` imposed across filled areas so that flats also drain.
#'
#' @param elev Numeric matrix or single-layer `SpatRaster` of elevations.
#' @param epsilon Minimum elevation drop imposed between neighbouring cells
#'   in filled areas (in elevation units).
#' @param max_iter Maximum number of sweeps.
#' @return The filled elevations, in the form of the input.
#' @references Planchon, O. and Darboux, F. (2002). A fast, simple and
#'   versatile algorithm to fill the depressions of digital elevation models.
#'   *Catena*, 46, 159-176. \doi{10.1016/S0341-8162(01)00164-3}
#' @examples
#' z <- outer(1:20, 1:20, "+"); z[10, 10] <- 0       # an artificial pit
#' sum(is.na(dm_flow_from_dem(z)))                    # 2 sinks
#' sum(is.na(dm_flow_from_dem(dm_fill_sinks(z))))     # 1 outlet
#' @export
dm_fill_sinks <- function(elev, epsilon = 1e-4, max_iter = 10000) {
  if (inherits(elev, "SpatRaster")) {
    if (!requireNamespace("terra", quietly = TRUE)) stop("Package 'terra' is needed.", call. = FALSE)
    m <- terra::as.matrix(elev[[1]], wide = TRUE)
    f <- dm_fill_sinks(m, epsilon, max_iter)
    out <- terra::rast(elev[[1]])
    terra::values(out) <- as.numeric(t(f))
    return(out)
  }
  Z <- as.matrix(elev)
  nr <- nrow(Z); nc <- ncol(Z)
  na <- is.na(Z)
  pad <- function(M, v) { P <- matrix(v, nr + 2, nc + 2); P[2:(nr + 1), 2:(nc + 1)] <- M; P }
  # boundary cells: grid edge or next to a missing value
  Pna <- pad(na, TRUE)
  border <- matrix(FALSE, nr, nc)
  dr <- c(-1, -1, -1, 0, 0, 1, 1, 1); dc <- c(-1, 0, 1, -1, 1, -1, 0, 1)
  for (k in 1:8) border <- border | Pna[(2:(nr + 1)) + dr[k], (2:(nc + 1)) + dc[k]]
  W <- ifelse(border | na, Z, Inf)
  for (it in seq_len(max_iter)) {
    W_old <- W
    for (k in 1:8) {
      Wn <- pad(W, Inf)[(2:(nr + 1)) + dr[k], (2:(nc + 1)) + dc[k]]
      Wn[is.na(Wn)] <- Inf
      cand <- pmax(Z, Wn + epsilon)
      upd <- !border & !na & cand < W
      W[upd] <- cand[upd]
    }
    if (identical(W, W_old)) break
    if (it == max_iter) warning("Sink filling did not converge in ", max_iter, " sweeps.",
                                call. = FALSE)
  }
  W[na] <- NA
  W
}
