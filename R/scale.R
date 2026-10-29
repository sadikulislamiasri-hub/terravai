# ---------------------------------------------------------------------------
# Scale operators: points -> blocks -> watersheds
# ---------------------------------------------------------------------------

#' Block-aggregation matrix (districts, fields, management zones)
#'
#' Builds the scale operator `H` that turns location-level predictions into
#' unit-level means (or sums), for decisions made per district, field or
#' management zone.
#'
#' @param group A vector assigning each location to a unit (e.g. district
#'   names). `NA` locations are left out of every unit.
#' @param fun `"mean"` (area-weighted average) or `"sum"` (total).
#' @param weights Optional non-negative weights per location (e.g. cell area or
#'   cropland fraction).
#' @return A numeric matrix with one row per unit (row names = unit labels)
#'   and one column per location.
#' @examples
#' dm_block_matrix(c("A", "A", "B", "B", "B"))
#' @export
dm_block_matrix <- function(group, fun = c("mean", "sum"), weights = NULL) {
  fun <- match.arg(fun)
  g <- as.factor(group)
  n <- length(g)
  w <- weights %||% rep(1, n)
  if (length(w) != n || any(w < 0, na.rm = TRUE)) {
    stop("'weights' must be non-negative with one value per location.", call. = FALSE)
  }
  H <- matrix(0, nlevels(g), n, dimnames = list(levels(g), NULL))
  ok <- !is.na(g)
  H[cbind(as.integer(g)[ok], which(ok))] <- w[ok]
  if (fun == "mean") H <- H / pmax(rowSums(H), .Machine$double.eps)
  H
}

#' Flow-accumulation matrix for watershed decisions
#'
#' Builds the scale operator `H` such that `H %*% z` is the quantity
#' delivered to every unit from all units upstream of it (including itself),
#' for example sediment or runoff accumulated along a drainage network.
#'
#' Each unit passes a fraction `delivery` of what reaches it to the next unit
#' downstream, so the contribution of unit `u` to unit `v` is the product of
#' the delivery ratios along the path from `u` to `v` (excluding `v`).
#'
#' @param to Integer vector: for each unit, the index of the unit it drains
#'   to, or `NA`/`0` for an outlet.
#' @param delivery Delivery ratio(s) between 0 and 1 (one value, or one per
#'   unit).
#' @return An `n x n` matrix. Select rows for the units where decisions are
#'   made, e.g. `H[outlet, , drop = FALSE]`.
#' @seealso [dm_flow_from_dem()]
#' @examples
#' # 1 -> 3, 2 -> 3, 3 -> 4 (outlet)
#' dm_flow_matrix(c(3, 3, 4, NA), delivery = 0.9)
#' @export
dm_flow_matrix <- function(to, delivery = 1) {
  n <- length(to)
  to[is.na(to)] <- 0L
  to <- as.integer(to)
  if (any(to < 0 | to > n) || any(to == seq_len(n))) {
    stop("'to' must index other units, with NA or 0 for outlets.", call. = FALSE)
  }
  if (length(delivery) == 1L) delivery <- rep(delivery, n)
  if (length(delivery) != n || any(delivery < 0 | delivery > 1)) {
    stop("'delivery' must lie in [0, 1], one value or one per unit.", call. = FALSE)
  }
  H <- diag(1, n)
  for (u in seq_len(n)) {
    w <- 1; v <- u; steps <- 0L
    while (to[v] != 0L) {
      w <- w * delivery[v]
      v <- to[v]
      H[v, u] <- H[v, u] + w
      steps <- steps + 1L
      if (steps > n) stop("'to' contains a cycle.", call. = FALSE)
    }
  }
  H
}

#' Flow directions from a gridded elevation model (D8)
#'
#' Derives single-direction (D8, steepest descent) flow routing from a matrix
#' of elevations, giving the `to` vector needed by [dm_flow_matrix()]. Cells
#' with no lower neighbour drain nowhere (outlets or pits).
#'
#' @param elev Numeric matrix of elevations (rows = y, columns = x), or a
#'   single-layer `SpatRaster`.
#' @return An integer vector `to` of length `length(elev)`, indexing cells in
#'   column-major order of the matrix (or raster cell order for a
#'   `SpatRaster`), with `NA` for cells that drain nowhere.
#' @examples
#' elev <- outer(1:5, 1:5, function(i, j) i + j)   # slopes towards cell [1, 1]
#' to <- dm_flow_from_dem(elev)
#' H <- dm_flow_matrix(to)
#' sum(H[1, ])   # every cell drains to the corner
#' @export
dm_flow_from_dem <- function(elev) {
  if (inherits(elev, "SpatRaster")) {
    if (!requireNamespace("terra", quietly = TRUE)) stop("Package 'terra' is needed.", call. = FALSE)
    m <- terra::as.matrix(elev, wide = TRUE)
    to_m <- dm_flow_from_dem(m)
    # matrix is row-major in raster cell order: convert index
    nr <- nrow(m); nc <- ncol(m)
    rc <- function(k) ((k - 1) %% nr) * nc + (k - 1) %/% nr + 1
    out <- integer(nr * nc)
    out[rc(seq_len(nr * nc))] <- ifelse(is.na(to_m), NA_integer_, rc(to_m))
    return(out)
  }
  elev <- as.matrix(elev)
  nr <- nrow(elev); nc <- ncol(elev)
  # vectorised D8: compare every cell with its 8 neighbours at once
  pad <- matrix(NA_real_, nr + 2, nc + 2)
  pad[2:(nr + 1), 2:(nc + 1)] <- elev
  dr <- c(-1, -1, -1, 0, 0, 1, 1, 1); dc <- c(-1, 0, 1, -1, 1, -1, 0, 1)
  dd <- sqrt(dr^2 + dc^2)
  S <- matrix(-Inf, nr * nc, 8)
  for (k in 1:8) {
    nb <- pad[(2:(nr + 1)) + dr[k], (2:(nc + 1)) + dc[k]]
    sl <- (elev - nb) / dd[k]
    sl[is.na(sl)] <- -Inf
    S[, k] <- as.numeric(sl)
  }
  best <- max.col(S, ties.method = "first")
  smax <- S[cbind(seq_len(nr * nc), best)]
  ii <- rep(seq_len(nr), nc) + dr[best]
  jj <- rep(seq_len(nc), each = nr) + dc[best]
  to <- as.integer((jj - 1L) * nr + ii)
  to[!(smax > 0) | is.na(as.numeric(elev))] <- NA_integer_
  to
}

#' Aggregate a predictive distribution to decision units
#'
#' Applies a scale operator `H` (from [dm_block_matrix()],
#' [dm_flow_matrix()] or your own) to a predictive distribution, giving the
#' predictive distribution of the decision quantities `Q = H Z` with their
#' correlation fully propagated. Ignoring correlation (aggregating
#' pixel-by-pixel standard deviations) usually misstates the uncertainty of
#' an aggregate.
#'
#' @param pred A `dm_predictive` object with a covariance matrix or draws.
#' @param H Scale operator, one row per decision unit and one column per
#'   location.
#' @param coords Optional coordinates of the decision units. By default the
#'   weighted centroid of each unit's locations.
#' @return A `dm_predictive` object for the decision units.
#' @examples
#' sim <- dm_simulate(n_side = 12, n_obs = 25, seed = 3)
#' fit <- dm_krige(sim$obs[, c("x", "y")], sim$obs$z, sim$grid[, c("x", "y")],
#'                 model = sim$model)
#' quad <- paste0(sim$grid$x > 0.5, sim$grid$y > 0.5)
#' dm_aggregate(fit, dm_block_matrix(quad))
#' @export
dm_aggregate <- function(pred, H, coords = NULL) {
  .check_pred(pred)
  H <- as.matrix(H)
  if (ncol(H) != length(pred$mean)) stop("'H' must have one column per location.", call. = FALSE)
  if (is.null(coords) && !is.null(pred$coords)) {
    W <- abs(H) / pmax(rowSums(abs(H)), .Machine$double.eps)
    coords <- data.frame(x = as.numeric(W %*% pred$coords$x),
                         y = as.numeric(W %*% pred$coords$y))
  }
  if (!is.null(pred$draws)) {
    out <- dm_draws(H %*% pred$draws, coords = coords)
  } else {
    Sigma <- .Sigma(pred)
    out <- dm_gaussian(as.numeric(H %*% pred$mean), cov = H %*% Sigma %*% t(H),
                       coords = coords)
  }
  out
}
