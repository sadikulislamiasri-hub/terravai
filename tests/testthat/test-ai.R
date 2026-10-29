test_that("spatial copula correlates quantile draws", {
  q <- matrix(rep(c(-1.64, 0, 1.64), each = 3), 3)
  xy <- cbind(c(0, 0.05, 5), 0)
  ind <- dm_quantiles(q, c(0.05, 0.5, 0.95), coords = xy, n_draws = 3000, seed = 1)
  cop <- dm_quantiles(q, c(0.05, 0.5, 0.95), coords = xy, n_draws = 3000, seed = 1,
                      copula = dm_covmodel("exponential", 1, 0.5))
  r_ind <- cor(t(ind$draws)); r_cop <- cor(t(cop$draws))
  expect_lt(abs(r_ind[1, 2]), 0.1)
  expect_gt(r_cop[1, 2], 0.8)
  expect_lt(abs(r_cop[1, 3]), 0.1)
  # marginals are unchanged by the copula
  expect_equal(quantile(cop$draws[1, ], 0.5, names = FALSE), 0, tolerance = 0.1)
  expect_error(dm_quantiles(q, c(0.05, 0.5, 0.95), copula = dm_covmodel()), "coords")
})

test_that("conformal distributions cover at the nominal rate", {
  set.seed(5)
  x <- runif(600); y <- sin(4 * x) + rnorm(600, 0, 0.3)
  pred <- sin(4 * x)
  cal <- 1:300; new <- 301:600
  p <- dm_conformal(pred[new], y[cal], pred[cal], n_draws = 400, seed = 2)
  lo <- apply(p$draws, 1, quantile, 0.05); hi <- apply(p$draws, 1, quantile, 0.95)
  cover <- mean(y[new] >= lo & y[new] <= hi)
  expect_gt(cover, 0.84); expect_lt(cover, 0.96)
  s <- dm_conformal(pred[new], y[cal], pred[cal], scale = rep(2, 300),
                    calib_scale = rep(1, 300), seed = 2)
  expect_gt(mean(s$sd), 1.5 * mean(p$sd))
  expect_error(dm_conformal(1, 1:5, 1:5), "10 calibration")
})

test_that("ensemble VOI agrees with the Gaussian closed form", {
  sim <- dm_simulate(n_side = 10, n_obs = 20, seed = 1)
  fit <- dm_krige(sim$obs[, c("x", "y")], sim$obs$z, sim$grid[, c("x", "y")], model = sim$model)
  l <- dm_loss(0, "above", ratio = 3)
  cand <- seq(1, 100, by = 9)
  vg <- dm_voi(fit, l, candidates = cand, noise = 0.05)
  set.seed(4)
  L <- chol(fit$cov + diag(1e-8, 100))
  D <- fit$mean + crossprod(L, matrix(rnorm(100 * 1500), 100))
  ve <- dm_voi(dm_draws(D, coords = sim$grid[, c("x", "y")]), l, candidates = cand,
               noise = 0.05, method = "ensemble", seed = 1)
  expect_equal(attr(ve, "method"), "ensemble")
  expect_gt(cor(vg$voi, ve$voi), 0.9)
  expect_equal(sum(ve$voi) / sum(vg$voi), 1, tolerance = 0.2)
  expect_error(dm_voi(fit, l, method = "ensemble"), "draws")
  # auto picks the ensemble method for draws
  expect_equal(attr(suppressWarnings(dm_voi(dm_draws(D[, 1:200]), l, candidates = 1:5,
                                             seed = 1)), "method"), "ensemble")
})

test_that("plain-language narration reports the decision", {
  sim <- dm_simulate(n_side = 10, n_obs = 20, seed = 2)
  fit <- dm_krige(sim$obs[, c("x", "y")], sim$obs$z, sim$grid[, c("x", "y")], model = sim$model)
  dec <- dm_decide(fit, dm_loss(0, "above", ratio = 2, action = "apply lime"))
  txt <- dm_narrate(dec, voi = dm_voi(fit, dm_loss(0, "above", ratio = 2)), area = "the block")
  expect_s3_class(txt, "dm_narrative")
  expect_match(txt, "apply lime")
  expect_match(txt, "the block")
  expect_match(txt, "most valuable")
  expect_output(print(txt), "Across the block")
})

test_that("quantile regression forest adapter returns spatial draws", {
  skip_if_not_installed("ranger")
  sim <- dm_simulate(n_side = 10, n_obs = 60, seed = 3)
  sim$grid$cov1 <- sim$grid$z + stats::rnorm(100, 0, 0.3)
  train <- cbind(sim$obs, cov1 = sim$grid$cov1[sim$obs$cell])
  p <- dm_fit_qrf(z ~ cov1, train, sim$grid, num.trees = 50, n_draws = 100, seed = 1)
  expect_s3_class(p, "dm_predictive")
  expect_equal(dim(p$draws), c(100, 100))
  expect_true(inherits(attr(p, "forest"), "ranger"))
  p2 <- dm_fit_qrf(z ~ cov1, train, sim$grid, num.trees = 50, n_draws = 100, seed = 1,
                   copula = FALSE)
  expect_null(attr(p2, "copula"))
  expect_error(dm_fit_qrf(z ~ cov1, train, sim$grid, coords = c("a", "b")), "Coordinate")
})
