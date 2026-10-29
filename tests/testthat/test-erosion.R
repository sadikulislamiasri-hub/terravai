test_that("curve number runoff matches hand calculation", {
  # P = 60 mm, CN = 80: S = 63.5, Ia = 12.7, Q = 47.3^2 / 110.8
  expect_equal(dm_runoff_cn(60, 80), 47.3^2 / (47.3 + 63.5), tolerance = 1e-6)
  expect_equal(dm_runoff_cn(10, 80), 0)                      # below initial abstraction
  expect_true(dm_runoff_cn(60, 80, lambda = 0.05) > dm_runoff_cn(60, 80))
  expect_error(dm_runoff_cn(-1, 80), "P")
  expect_error(dm_runoff_cn(50, 120), "CN")
  expect_equal(dm_cn_amc(75, "III"), 23 * 75 / (10 + 0.13 * 75))
  expect_equal(dm_cn_amc(75, "I"), 4.2 * 75 / (10 - 0.058 * 75))
})

test_that("RUSLE factors behave as published", {
  expect_equal(dm_rusle(800, 0.03, 2.5, 0.3, 0.5), 800 * 0.03 * 2.5 * 0.3 * 0.5)
  expect_equal(dm_c_ndvi(0), 1)
  expect_true(all(diff(dm_c_ndvi(c(0.1, 0.3, 0.6))) < 0))
  # Wischmeier-Smith LS = 1 on the unit plot (22.13 m, 9 % slope)
  expect_equal(dm_ls_factor(atan(0.09) * 180 / pi, length = 22.13, method = "wischmeier"),
               1, tolerance = 0.01)
  expect_true(dm_ls_factor(10, flow_acc = 50, cell_size = 30) >
                dm_ls_factor(10, flow_acc = 5, cell_size = 30))
  k <- dm_k_epic(20, 60, 20, 1.5)
  expect_true(k > 0.01 && k < 0.08)
  expect_true(dm_k_epic(20, 60, 20, 3) < k)                   # more organic carbon, less erodible
  expect_equal(dm_fournier(rbind(rep(100, 12))), 100)
  expect_equal(dm_fournier(rbind(rep(100, 12)), a = 10, b = 2), 210)
  expect_equal(dm_musle(1000, 1, 0.03, 2, 0.3), 11.8 * 1000^0.56 * 0.03 * 2 * 0.3)
  expect_equal(dm_peak_rational(0.5, 36, 10), 0.5)
})

test_that("nutrient mass balance and enrichment", {
  expect_equal(dm_nutrient_sediment(10, 0.1, "percent", 2), 20)
  expect_equal(dm_nutrient_sediment(10, 500, "mg_kg"), 5)
  expect_equal(dm_nutrient_runoff(50, 2), 1)
  er <- dm_enrichment(c(10, 1000), a = 7.4, b = -0.2)
  expect_true(er[1] > er[2] && all(er >= 1))
  expect_error(dm_enrichment(10), "coefficients")
})

test_that("empirical N, P, K equations reproduce the published coefficients", {
  expect_equal(dm_nutrient_empirical(3, "N", "soil"), 2.3647 * 3)
  expect_equal(dm_nutrient_empirical(14, "N", "soil"), 0.7909 * 14 - 0.2909)
  expect_equal(dm_nutrient_empirical(25, "P", "runoff"), 0.027 * 25 - 0.0066)
  expect_equal(dm_nutrient_empirical(10, "K", "runoff"), 0.4252 * 10 - 0.1895)
  expect_equal(dm_nutrient_empirical(c(3, 14), "N", "soil"),
               c(2.3647 * 3, 0.7909 * 14 - 0.2909))
  expect_error(dm_nutrient_empirical(-1, "N"), "x")
  expect_error(dm_nutrient_empirical(120, "N", "runoff"), "x")
  cf <- dm_nutrient_coef()
  expect_equal(nrow(cf), 24)
  expect_equal(nrow(attr(cf, "jumps")), 18)
})

test_that("uncertainty propagates through formulas into decisions", {
  A <- dm_rusle(R = dm_uncertain(800, cv = 0.15, lower = 0), K = 0.03, LS = 2.5,
                C = dm_uncertain(0.3, sd = 0.05, lower = 0, upper = 1), n_draws = 2000, seed = 1)
  expect_s3_class(A, "dm_predictive")
  expect_equal(A$mean, 18, tolerance = 0.05)
  d <- dm_decide(A, dm_loss_tolerance())
  expect_true(d$act)
  ln <- dm_uncertain(5, sd = 2, dist = "lognormal")
  expect_true(all(dm_runoff_cn(ln, 80, seed = 1)$draws >= 0))
  tn <- dm_uncertain(c(0.9, 0.2), sd = 0.5, lower = 0, upper = 1)
  D <- terravai:::.draw_uncertain(tn, 500)
  expect_true(all(D >= 0 & D <= 1))
  N <- dm_nutrient_empirical(dm_uncertain(14, sd = 3, lower = 0), "N", seed = 1)
  expect_s3_class(N, "dm_predictive")
  expect_s3_class(dm_loss_nutrient(10, "N"), "dm_loss")
})

test_that("simulated watershed routes loads to one outlet", {
  ws <- dm_simulate_watershed(n_side = 15, n_plots = 10, seed = 2)
  expect_s3_class(ws, "dm_watershed")
  expect_equal(length(ws$outlet), 1)
  expect_equal(nrow(ws$plots), 10)
  expect_true(all(ws$grid$soil_loss >= 0))
  acc <- dm_flow_accumulation(ws$to)
  expect_equal(acc[ws$outlet], nrow(ws$grid))
  ld <- dm_outlet_load(rep(1, nrow(ws$grid)), ws$to)
  expect_equal(unname(ld), nrow(ws$grid))
  p <- dm_gaussian(ws$grid$soil_loss, sd = 1, coords = ws$grid[, c("x", "y")])
  expect_warning(ol <- dm_outlet_load(p, ws$to, delivery = 0.9), "independent")
  expect_s3_class(ol, "dm_predictive")
  expect_output(print(ws), "outlet")
})

test_that("raster inputs return rasters", {
  skip_if_not_installed("terra")
  r <- terra::rast(nrows = 3, ncols = 3, vals = seq(10, 90, by = 10))
  out <- dm_runoff_cn(r, 80)
  expect_s4_class(out, "SpatRaster")
  expect_equal(terra::values(out, mat = FALSE), dm_runoff_cn(seq(10, 90, by = 10), 80))
})
