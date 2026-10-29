# ---------------------------------------------------------------------------
# Base-graphics map helpers (no plotting dependencies)
# ---------------------------------------------------------------------------

.grid_index <- function(coords) {
  ux <- sort(unique(round(coords$x, 10)))
  uy <- sort(unique(round(coords$y, 10)))
  if (length(ux) * length(uy) < nrow(coords) || length(ux) < 2 || length(uy) < 2) return(NULL)
  if (length(ux) * length(uy) > 4 * nrow(coords)) return(NULL)
  dx <- diff(ux); dy <- diff(uy)
  if (max(dx) - min(dx) > 1e-6 * max(dx) || max(dy) - min(dy) > 1e-6 * max(dy)) return(NULL)
  list(ux = ux, uy = uy, ix = match(round(coords$x, 10), ux),
       iy = match(round(coords$y, 10), uy))
}

.map_plot <- function(coords, z, main = "", palette = "Viridis", rev = FALSE,
                      discrete = FALSE, labels = NULL, cols = NULL, zlim = NULL,
                      pts = NULL, pts_labels = NULL, legend_note = NULL, cex_pts = 0.9) {
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  graphics::layout(matrix(1:2, 1), widths = c(5, 1.25))
  graphics::par(mar = c(2.5, 2.5, 3.2, 0.5), mgp = c(1.5, 0.5, 0), cex.main = 1)
  ok <- is.finite(z)
  if (discrete) {
    k <- length(labels)
    cols <- cols %||% grDevices::hcl.colors(k, "Dark 3")
    zlim <- c(0, k - 1)
    pal <- cols
    brk <- seq(-0.5, k - 0.5, by = 1)
  } else {
    zlim <- zlim %||% range(z[ok])
    if (diff(zlim) == 0) zlim <- zlim + c(-0.5, 0.5)
    pal <- grDevices::hcl.colors(64, palette, rev = rev)
    brk <- seq(zlim[1], zlim[2], length.out = 65)
  }
  zc <- pmin(pmax(z, zlim[1]), zlim[2])
  gi <- .grid_index(coords)
  if (!is.null(gi)) {
    m <- matrix(NA_real_, length(gi$ux), length(gi$uy))
    m[cbind(gi$ix, gi$iy)] <- zc
    graphics::image(gi$ux, gi$uy, m, col = pal, breaks = brk, asp = 1,
                    xlab = "x", ylab = "y", main = main, useRaster = TRUE)
  } else {
    idx <- findInterval(zc, brk, all.inside = TRUE)
    graphics::plot.new()
    graphics::plot.window(range(coords$x), range(coords$y), asp = 1)
    graphics::axis(1); graphics::axis(2); graphics::box()
    graphics::title(main = main, xlab = "x", ylab = "y")
    graphics::points(coords$x, coords$y, pch = 16, col = pal[idx], cex = 1.1)
  }
  if (!is.null(legend_note)) graphics::mtext(legend_note, side = 3, line = 0.2, cex = 0.7)
  if (!is.null(pts)) {
    graphics::points(pts$x, pts$y, pch = 21, bg = "white", col = "black", cex = cex_pts * 1.6)
    if (!is.null(pts_labels)) graphics::text(pts$x, pts$y, pts_labels, cex = 0.55)
  }
  # legend panel
  graphics::par(mar = c(2.5, 0.2, 3, 0.5))
  graphics::plot.new()
  if (discrete) {
    graphics::legend("center", legend = labels, fill = cols, bty = "n", cex = 0.85)
  } else {
    graphics::plot.window(c(0, 1), c(0, 1))
    yb <- seq(0.1, 0.85, length.out = 65)
    graphics::rect(0.05, yb[-65], 0.35, yb[-1], col = pal, border = NA)
    tk <- pretty(zlim, 5); tk <- tk[tk >= zlim[1] & tk <= zlim[2]]
    ty <- 0.1 + (tk - zlim[1]) / diff(zlim) * 0.75
    graphics::text(0.42, ty, format(tk), adj = 0, cex = 0.75, xpd = NA)
  }
  invisible(NULL)
}
