# Draws the TerravAI hex sticker with base R graphics: Rscript tools/make_logo.R
make_logo <- function(file = "man/figures/logo.png", size = 600) {
  grDevices::png(file, width = size, height = round(size * 1.1547), bg = "transparent", res = 144)
  op <- graphics::par(mar = c(0, 0, 0, 0)); on.exit({ graphics::par(op); grDevices::dev.off() })
  graphics::plot.new(); graphics::plot.window(c(-1, 1), c(-1.1547, 1.1547), asp = 1)
  a <- seq(90, 450, by = 60) * pi / 180
  hx <- cos(a) * 1.1547 * 0.98; hy <- sin(a) * 1.1547 * 0.98
  graphics::polygon(hx, hy, col = "#F7F5F0", border = "#2F6690", lwd = 9)
  # terrain contours, clipped to the hexagon interior by drawing within its band
  g <- seq(-0.85, 0.85, length.out = 120)
  z <- outer(g, g, function(x, y) exp(-((x + 0.25)^2 + (y - 0.05)^2) / 0.12) +
               0.7 * exp(-((x - 0.35)^2 + (y + 0.15)^2) / 0.08) + 0.15 * y)
  inside <- outer(g, g, function(x, y) abs(x) < 0.84 & abs(y) < 0.97 - abs(x) * 0.577 & y > -0.5 & y < 0.55)
  z[!inside] <- NA
  graphics::contour(g, g, z, nlevels = 9, add = TRUE, drawlabels = FALSE,
                    col = grDevices::colorRampPalette(c("#6B4F3A", "#2F6690"))(9), lwd = 2.2)
  # decision boundary (threshold contour) highlighted
  graphics::contour(g, g, z, levels = 0.55, add = TRUE, drawlabels = FALSE, col = "#14B8A6", lwd = 5)
  # value-of-information pin
  px <- 0.12; py <- 0.13
  graphics::polygon(px + c(-0.09, 0, 0.09), py + c(0.05, -0.16, 0.05), col = "#D9822B", border = NA)
  graphics::symbols(px, py + 0.07, circles = 0.095, inches = FALSE, add = TRUE, bg = "#D9822B", fg = NA)
  graphics::symbols(px, py + 0.07, circles = 0.04, inches = FALSE, add = TRUE, bg = "#F7F5F0", fg = NA)
  # wordmark
  cx <- 2.5
  w1 <- graphics::strwidth("Terrav", cex = cx, font = 2)
  w2 <- graphics::strwidth("AI", cex = cx, font = 2)
  x0 <- -(w1 + w2) / 2
  graphics::text(x0, -0.68, "Terrav", adj = c(0, 0.5), cex = cx, font = 2, col = "#6B4F3A")
  graphics::text(x0 + w1, -0.68, "AI", adj = c(0, 0.5), cex = cx, font = 2, col = "#14B8A6")
  graphics::text(0, 0.78, "VALUE  OF  INFORMATION", cex = 0.95, font = 2, col = "#2F6690")
  invisible(file)
}
make_logo()
