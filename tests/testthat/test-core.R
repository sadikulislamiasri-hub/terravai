test_that("predictive constructors validate and summarise", {
  p <- dm_gaussian(c(1, 2, 3), sd = 0.5)
  expect_s3_class(p, "dm_predictive")
  expect_equal(length(p), 3)
  expect_error(dm_gaussian(c(1, 2), sd = c(1, 2, 3)), "same length")
  expect_error(dm_gaussian(1, sd = -1), "sd")
  S <- matrix(c(1, 0.5, 0.5, 2), 2)
  expect_equal(dm_gaussian(c(0, 0), cov = S)$sd, sqrt(c(1, 2)))
  d <- dm_draws(matrix(rnorm(40), 2))
  expect_equal(d$type, "draws")
  q <- dm_quantiles(rbind(c(1, 2, 3)), c(0.1, 0.5, 0.9), n_draws = 4000, seed = 1)
  expect_equal(median(q$draws[1, ]), 2, tolerance = 0.03)
})

test_that("loss and cut-off agree", {
  l <- dm_loss(5, "below", ratio = 4)
  expect_equal(l$cost_fn, 4)
  expect_equal(dm_cutoff(4), 0.2)
  expect_equal(dm_cutoff(cutoff = 0.2), 4)
  expect_error(dm_loss(5, ratio = -1), "positive")
})

test_that("Gaussian expected losses match Monte Carlo", {
  set.seed(10)
  D <- matrix(rnorm(4e5, 0.2, 0.7), nrow = 1)
  for (type in c("fixed", "linear")) for (dir in c("above", "below")) {
    l <- dm_loss(0.5, dir, ratio = 3, type = type)
    g <- dm_decide(dm_gaussian(0.2, sd = 0.7), l)
    m <- dm_decide(dm_draws(D), l)
    expect_equal(g$loss_act, m$loss_act, tolerance = 0.01)
    expect_equal(g$loss_no, m$loss_no, tolerance = 0.01)
  }
})

test_that("decision rule follows the probability cut-off", {
  l <- dm_loss(0, "above", ratio = 4)            # act if P(Z > 0) > 0.2
  p <- dm_gaussian(qnorm(c(0.1, 0.3)), sd = 1)   # P = 0.1, 0.3
  d <- dm_decide(p, l)
  expect_equal(d$p_need, c(0.1, 0.3), tolerance = 1e-8)
  expect_equal(d$act, c(FALSE, TRUE))
  expect_equal(d$breakeven, c(9, 0.7 / 0.3), tolerance = 1e-8)
  expect_true(all(d$exp_loss <= pmin(d$loss_act, d$loss_no) + 1e-12))
})

test_that("VOI of a perfect measurement equals EVPI and is non-negative", {
  l <- dm_loss(0.5, "above", ratio = 2, type = "linear")
  p <- dm_gaussian(0.3, cov = matrix(0.16))
  expect_equal(dm_voi(p, l)$voi, dm_decide(p, l)$exp_loss, tolerance = 1e-6)
  noisy <- dm_voi(p, l, noise = 0.5)$voi
  expect_true(noisy > 0 && noisy < dm_decide(p, l)$exp_loss)
  far <- dm_gaussian(10, cov = matrix(0.01))
  expect_equal(dm_voi(far, l)$voi, 0, tolerance = 1e-10)
})

test_that("kriging, VOI and design run end to end", {
  sim <- dm_simulate(n_side = 10, n_obs = 25, seed = 1)
  expect_s3_class(sim, "dm_sim")
  cm <- dm_fit_cov(sim$obs[, c("x", "y")], sim$obs$z)
  expect_s3_class(cm, "dm_covmodel")
  fit <- dm_krige(sim$obs[, c("x", "y")], sim$obs$z, sim$grid[, c("x", "y")], model = cm)
  expect_equal(dim(fit$cov), c(100, 100))
  # kriging interpolates observed cells up to the nugget
  expect_true(all(fit$sd[sim$obs$cell] < sqrt(cm$psill)))
  l <- dm_loss(0, "above", ratio = 3)
  v <- dm_voi(fit, l, noise = cm$nugget)
  expect_equal(nrow(v), 100)
  expect_true(all(v$voi >= 0))
  des <- dm_design(fit, l, n = 5, noise = cm$nugget)
  des2 <- dm_design(fit, l, n = 5, noise = cm$nugget, lazy = FALSE)
  rl <- attr(des, "remaining_loss")
  expect_true(all(diff(rl) <= 1e-10))
  # value of information is usually, but not strictly, diminishing
  expect_true(all(diff(des2$voi) <= 0.02 * des2$voi[-length(des2$voi)]))
  expect_gt(sum(des$voi) / sum(des2$voi), 0.99)
  stop_des <- dm_design(fit, l, n = 50, noise = cm$nugget, stop_below = max(v$voi) * 2)
  expect_equal(nrow(stop_des), 0)
  b <- dm_design(fit, l, n = 50, cost = rep(2, 100), budget = 6, noise = cm$nugget)
  expect_equal(nrow(b), 3)
})

test_that("random state is restored by dm_simulate", {
  set.seed(42); a <- runif(1)
  set.seed(42); invisible(dm_simulate(5, 5, seed = 1)); b <- runif(1)
  expect_equal(a, b)
})

test_that("scale operators aggregate correctly", {
  H <- dm_block_matrix(c("a", "a", "b"))
  expect_equal(unname(H), rbind(c(0.5, 0.5, 0), c(0, 0, 1)))
  Hs <- dm_block_matrix(c("a", "a", "b"), fun = "sum")
  expect_equal(unname(rowSums(Hs)), c(2, 1))
  F <- dm_flow_matrix(c(3, 3, 4, NA), delivery = 0.5)
  expect_equal(F[4, ], c(0.25, 0.25, 0.5, 1))
  expect_error(dm_flow_matrix(c(2, 1)), "cycle")
  elev <- outer(1:4, 1:4, "+")
  to <- dm_flow_from_dem(elev)
  expect_equal(sum(is.na(to)), 1)
  expect_equal(sum(dm_flow_matrix(to)[which(is.na(to)), ]), 16)
  S <- diag(3); S[1, 2] <- S[2, 1] <- 0.5
  p <- dm_gaussian(c(1, 2, 3), cov = S, coords = cbind(1:3, 0))
  a <- dm_aggregate(p, H)
  expect_equal(a$mean, c(1.5, 3))
  expect_equal(a$sd[1], sqrt((1 + 1 + 2 * 0.5) / 4))
  # VOI for an aggregate decision uses the operator
  v <- dm_voi(p, dm_loss(1.6, "above"), target = H[1, , drop = FALSE])
  expect_equal(v$voi[3], 0)
  expect_true(v$voi[1] > 0)
})

test_that("evaluation detects ranking reversals", {
  set.seed(3)
  obs <- rnorm(500)
  near <- obs + ifelse(abs(obs) < 0.5, -sign(obs) * 0.6, 0)  # errors only near 0
  far <- obs + ifelse(abs(obs) > 1.5, sign(obs) * 1.2, 0)    # larger errors far from 0
  ev <- dm_evaluate(obs, list(near = near, far = far), dm_loss(0, "above"))
  s <- summary(ev)
  expect_equal(s$best_by_rmse, "near")
  expect_equal(s$best_by_regret, "far")
  expect_true(s$reversal)
  cur <- dm_threshold_curve(obs, list(near = near, far = far), c(-1, 0, 1))
  expect_s3_class(cur, "dm_curve")
  expect_equal(nrow(cur), 6)
})

test_that("plots and raster output work", {
  sim <- dm_simulate(n_side = 8, n_obs = 15, seed = 2)
  fit <- dm_krige(sim$obs[, c("x", "y")], sim$obs$z, sim$grid[, c("x", "y")], model = sim$model)
  dec <- dm_decide(fit, dm_loss(0, "above", ratio = 2))
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  on.exit({ grDevices::dev.off(); unlink(f) }, add = TRUE)
  for (w in c("act", "p_need", "p_wrong", "exp_loss", "breakeven", "mean")) {
    expect_invisible(plot(dec, w))
  }
  expect_invisible(plot(dm_voi(fit, dm_loss(0, "above"))))
  expect_invisible(plot(dm_design(fit, dm_loss(0, "above"), n = 3)))
  ev <- dm_evaluate(sim$grid$z, list(k = fit, m = fit$mean), dm_loss(0, "above"),
                    ratios = c(1, 4))
  expect_invisible(plot(ev))
  expect_output(print(summary(dec)), "Decision summary")
  expect_output(print(summary(ev)), "Accuracy versus decision quality")
  skip_if_not_installed("terra")
  r <- dm_as_raster(dec)
  expect_s4_class(r, "SpatRaster")
  expect_equal(terra::nlyr(r), 5)
})
