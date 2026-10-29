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
