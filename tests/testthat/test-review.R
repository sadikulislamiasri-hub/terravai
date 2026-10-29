test_that("sink filling removes pits and never lowers the surface", {
  z <- outer(1:20, 1:20, function(i, j) i + j); z[10, 10] <- 0; z[5, 15] <- -3
  f <- dm_fill_sinks(z)
  expect_true(all(f >= z))
  to <- dm_flow_from_dem(f)
  expect_equal(sum(is.na(to)), sum(is.na(dm_flow_from_dem(outer(1:20, 1:20, `+`)))))
})

test_that("mismatched input lengths are an error, not recycled", {
  expect_error(dm_rusle(R = 1:3, K = 0.03, LS = 1:2, C = 0.2, P = 1), "length")
})

test_that("Gaussian predictive covariance is kept in formulas", {
  S <- matrix(0.95, 2, 2); diag(S) <- 1
  p <- dm_gaussian(mean = c(10, 10), cov = S)
  out <- dm_rusle(R = p, K = 1, LS = 1, C = 1, P = 1, n_draws = 4000, seed = 1)
  expect_gt(stats::cor(out$draws[1, ], out$draws[2, ]), 0.9)
})

test_that("outlet load routes draws without a dense matrix", {
  ws <- dm_simulate_watershed(n_side = 15, seed = 1)
  d <- dm_draws(matrix(rexp(225 * 50), 225, 50))
  out <- dm_outlet_load(d, ws$to)
  expect_equal(ncol(out$draws), 50)
  expect_equal(sum(out$mean), sum(d$mean), tolerance = 1e-8)
})

test_that("ensemble VOI warns with few draws", {
  d <- dm_draws(matrix(rnorm(20 * 200, 10, 3), 20, 200))
  expect_warning(dm_voi(d, dm_loss_tolerance(), method = "ensemble", seed = 1), "draws")
})

test_that("steep slopes are flagged", {
  expect_warning(dm_ls_factor(slope = rep(40, 10), flow_acc = 10, cell_size = 30), "exceed")
})

test_that("daily curve number runoff uses antecedent moisture", {
  rain <- c(0, 60, 60, 60, 0, 0, 0, 0, 0, 0, 40)
  out <- dm_runoff_cn_daily(rain, CN = 75)
  expect_equal(as.character(out$amc[3]), "III")
  expect_true(out$total <= sum(rain))
  expect_true(all(out$daily >= 0))
})

test_that("compiled sink fill matches Planchon-Darboux and handles NA", {
  set.seed(2)
  z <- outer(1:30, 1:30, "+") + matrix(rnorm(900, 0, 2), 30)
  z[10:12, 10:12] <- NA
  f <- dm_fill_sinks(z)
  expect_true(all(is.na(f) == is.na(z)))
  expect_true(all(f >= z, na.rm = TRUE))
  # every non-NA interior cell now has a downslope neighbour
  to <- dm_flow_from_dem(f)
  border <- row(z) %in% c(1, 30) | col(z) %in% c(1, 30)
  nearna <- as.vector(Reduce(`|`, lapply(list(c(-1,0),c(1,0),c(0,-1),c(0,1),c(-1,-1),c(1,1),c(-1,1),c(1,-1)),
    function(d) { P <- matrix(FALSE, 32, 32); P[2:31, 2:31] <- is.na(z); P[(2:31) + d[1], (2:31) + d[2]] })))
  expect_true(all(!is.na(to[!border & !nearna & !is.na(as.vector(z))])))
  expect_error(dm_fill_sinks(z, epsilon = -1))
})

test_that("Nearing steepness factor is continuous and bounded on steep slopes", {
  s <- suppressWarnings(dm_ls_factor(c(1, 30, 60, 80), flow_acc = 1, cell_size = 22.13,
                                     steepness = "nearing"))
  expect_true(all(diff(s) > 0))
  expect_lt(max(s), 15.5)
  expect_silent(dm_ls_factor(rep(40, 5), flow_acc = 1, cell_size = 22.13, steepness = "nearing"))
})

test_that("gully risk indices", {
  g <- dm_gully_risk(c(2, 15), c(50, 2000), threshold = 5)
  expect_equal(g$at_risk, c(FALSE, TRUE))
  expect_equal(g$spi[2], 2000 * tan(15 * pi / 180))
})

test_that("local kriging agrees with global kriging", {
  sim <- dm_simulate(n_side = 15, n_obs = 120, seed = 5)
  xy <- sim$grid[, c("x", "y")]; o <- sim$obs
  f1 <- dm_krige(o[, c("x", "y")], o$z, xy, sim$model, full_cov = FALSE)
  f2 <- dm_krige(o[, c("x", "y")], o$z, xy, sim$model, full_cov = FALSE, nmax = 60)
  expect_gt(cor(f1$mean, f2$mean), 0.999)
  expect_lt(max(abs(f1$sd - f2$sd)), 0.05)
})

test_that("validation statistics are correct", {
  v <- dm_validate(c(1, 2, 3, 4), c(1, 2, 3, 4))
  expect_equal(unname(v$stats[c("nse", "kge", "rmse", "pbias")]), c(1, 1, 0, 0))
  p <- dm_gaussian(mean = rep(0, 2000), sd = rep(1, 2000))
  set.seed(1); v2 <- dm_validate(rnorm(2000), p)
  expect_equal(v2$coverage$observed, c(0.5, 0.9), tolerance = 0.03)
  expect_equal(unname(v2$stats["crps"]), 1 / sqrt(pi), tolerance = 0.08)
})

test_that("regional erosivity equations", {
  tab <- dm_erosivity_table()
  expect_equal(nrow(tab), 14)
  expect_equal(dm_erosivity_india(1000), 79 + 363)
  expect_equal(dm_erosivity_india(1000, equation = "bundelkhand_jhansi"), -1670.38 + 6190)
  expect_warning(dm_erosivity_india(300, equation = "bundelkhand_jhansi"), "range")
  expect_error(dm_erosivity_india(1000, equation = "nowhere"))
})

test_that("nutrient transfer model recovers parameters and narrows with local plots", {
  set.seed(3)
  site <- rep(letters[1:6], each = 40)
  A <- rlnorm(240, 1.5, 0.8)
  mu <- setNames(rnorm(6, 0, 0.6), letters[1:6])[site]
  N <- exp(mu + 0.5 * log(A) + rnorm(240, 0, 0.4))
  p0 <- dm_nutrient_transfer(A, N, site, new_soil_loss = c(2, 5, 10), seed = 1)
  expect_equal(unname(attr(p0, "params")["b"]), 0.5, tolerance = 0.1)
  expect_equal(unname(attr(p0, "params")["sigma"]), 0.4, tolerance = 0.15)
  p5 <- dm_nutrient_transfer(A, N, site, new_soil_loss = c(2, 5, 10),
                             local_soil_loss = c(2, 4, 6, 8, 10), local_nutrient = c(3, 4, 5, 6, 6),
                             seed = 1)
  expect_lt(attr(p5, "site_effect")["var"], attr(p0, "site_effect")["var"])
  expect_error(dm_nutrient_transfer(A[1:80], N[1:80], site[1:80], 5), "3 sites")
})
