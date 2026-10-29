test_that("Monte Carlo error shrinks at the 1/sqrt(M) rate", {
  f <- function(M, s) dm_rusle(R = dm_uncertain(800, cv = 0.2, lower = 0), K = 0.03, LS = 2,
                               C = 0.25, n_draws = M, seed = s)$mean
  e100 <- sd(sapply(1:60, function(s) f(100, s)))
  e1600 <- sd(sapply(1:60, function(s) f(1600, s)))
  expect_equal(e100 / e1600, 4, tolerance = 0.35)
})

test_that("closed-form VOI quadrature converges and matches Monte Carlo", {
  p <- dm_gaussian(0.3, cov = matrix(0.16))
  for (type in c("fixed", "linear")) for (dir in c("above", "below")) {
    l <- dm_loss(0.5, dir, ratio = 2, type = type)
    ref <- dm_voi(p, l, noise = 0.1, n_quad = 200)$voi
    expect_equal(dm_voi(p, l, noise = 0.1, n_quad = 20)$voi, ref, tolerance = 1e-3)
    set.seed(1)
    vr <- 0.16^2 / 0.26
    mp <- 0.3 + sqrt(vr) * rnorm(4e5)
    el <- terravai:::.exp_loss_gauss(mp, sqrt(0.16 - vr), l)
    mc <- dm_decide(p, l)$exp_loss - mean(pmin(el$act, el$no))
    expect_equal(ref, mc, tolerance = 0.01)
  }
})

test_that("covariance fitting reports convergence and boundary solutions", {
  sm <- dm_simulate(10, 80, model = dm_covmodel("exponential", 1, 0.15, 0.1), seed = 1)
  cm <- dm_fit_cov(sm$obs[, 1:2], sm$obs$z)
  expect_equal(cm$convergence, 0)
  expect_false(any(cm$at_bound))
  sm4 <- dm_simulate(10, 80, model = dm_covmodel("exponential", 1, 0.15, 0.1), seed = 4)
  expect_warning(dm_fit_cov(sm4$obs[, 1:2], sm4$obs$z), "search limit")
})

test_that("kriging handles one prediction location and duplicate points", {
  xy <- rbind(c(0, 0), c(0, 0), c(1, 1), c(0.5, 0.2))
  expect_warning(k <- dm_krige(xy, c(1, 1.1, 2, 1.5), rbind(c(0.2, 0.2)),
                               dm_covmodel("exponential", 1, 0.3, 0)), "positive definite")
  expect_length(k$mean, 1)
})

test_that("input uncertainty can be independent, common or spatial", {
  xy <- cbind(c(0, 0.05, 2), 0)
  draw <- function(cr) terravai:::.draw_uncertain(
    dm_uncertain(rep(800, 3), cv = 0.2, lower = 0, coords = xy, correlation = cr), 5000)
  set.seed(1)
  r_ind <- cor(t(draw("independent"))); r_com <- cor(t(draw("common")))
  r_sp <- cor(t(draw(dm_covmodel("exponential", 1, 0.3))))
  expect_lt(abs(r_ind[1, 2]), 0.05)
  expect_gt(r_com[1, 3], 0.99)
  expect_gt(r_sp[1, 2], 0.7); expect_lt(abs(r_sp[1, 3]), 0.05)
  # common errors do not average out over an area
  tot <- function(cr) {
    u <- dm_uncertain(rep(10, 50), sd = 2, correlation = cr)
    sd(colSums(terravai:::.draw_uncertain(u, 4000)))
  }
  set.seed(2)
  expect_gt(tot("common") / tot("independent"), 5)
  expect_error(dm_uncertain(1:3, sd = 1, correlation = dm_covmodel()), "coords")
})

test_that("convergence diagnostic recommends enough draws", {
  mk <- function(M) dm_rusle(R = dm_uncertain(800, cv = 0.2, lower = 0), K = 0.03, LS = 2,
                             C = dm_uncertain(0.25, sd = 0.07, lower = 0, upper = 1),
                             n_draws = M, seed = 1)
  small <- dm_convergence(mk(300), dm_loss_tolerance())
  big <- dm_convergence(mk(5000), dm_loss_tolerance())
  expect_s3_class(small, "dm_convergence")
  expect_gt(small$recommended, 300)
  expect_lte(big$recommended, 5000)
  expect_equal(big$mcse_p, sqrt(big$p * (1 - big$p) / 5000))
  expect_output(print(small), "increase n_draws")
  expect_output(print(big), "sufficient")
  lin <- dm_convergence(mk(500), dm_loss_tolerance(type = "linear"))
  expect_true(is.logical(lin$unstable))
  expect_error(dm_convergence(dm_gaussian(1, sd = 1)), "draws")
  f <- tempfile(fileext = ".pdf"); grDevices::pdf(f); on.exit({ grDevices::dev.off(); unlink(f) })
  expect_invisible(plot(small))
})
