test_that("particle-size conversion to USDA", {
  # USDA data passed through unchanged
  u <- dm_psd_to_usda(rbind(c(20, 60, 100)), c(0.002, 0.05, 2))
  expect_equal(unlist(u), c(sand = 40, silt = 40, clay = 20))
  # ISSS: log-linear interpolation between 0.02 and 0.2 mm
  r <- dm_psd_isss(clay = 25, silt = 15, fine_sand = 35, coarse_sand = 25)
  f <- (log10(0.05) - log10(0.02)) / (log10(0.2) - log10(0.02))
  expect_equal(r$silt, 15 + f * 35)
  expect_equal(r$clay, 25)
  expect_equal(r$sand + r$silt + r$clay, 100)
  expect_warning(dm_psd_isss(25, 15, sand = 60), "less accurate")
  # IS 1498: USDA silt is part of IS silt
  is <- dm_psd_is1498(20, 40, 25, 15)
  expect_true(is$silt < 40 && is$silt > 30)
  expect_error(dm_psd_to_usda(rbind(c(60, 20, 100)), c(0.002, 0.05, 2)), "decrease")
  expect_error(dm_psd_to_usda(rbind(c(20, 100)), c(0.01, 2)), "0.002")
})

test_that("Soil Health Card ratings", {
  expect_equal(as.character(dm_shc_rating(c(180, 280, 560, 620), "N")),
               c("low", "medium", "high", "high"))
  expect_equal(as.character(dm_shc_rating(c(0.3, 0.6, 1), "OC")), c("low", "medium", "high"))
  expect_equal(as.character(dm_shc_rating(c(0.3, 1.2), "Zn")), c("deficient", "sufficient"))
  pr <- dm_shc_rating(dm_gaussian(c(5, 17.5, 40), sd = 2), "P")
  expect_equal(rowSums(pr), rep(1, 3))
  expect_gt(pr$p_low[1], 0.9); expect_gt(pr$p_high[3], 0.9)
  expect_error(dm_shc_rating(1, "Q"), "Unknown")
  expect_equal(nrow(dm_shc_table()), 10)
})

test_that("Indian erosivity, slope and erosion classes", {
  expect_equal(dm_erosivity_india(1000), 79 + 363)
  expect_s3_class(dm_erosivity_india(dm_uncertain(1200, cv = 0.1, lower = 0), seed = 1),
                  "dm_predictive")
  expect_equal(as.character(dm_slope_class_india(c(0.5, 2, 4, 7, 12, 20, 30, 40, 70))),
               LETTERS[1:9])
  expect_match(as.character(dm_slope_class_india(4, labels = TRUE)), "gently sloping")
  expect_equal(as.character(dm_erosion_class(c(3, 8, 20, 45, 90))),
               c("E0 nil", "E1 slight", "E2 moderate", "E3 severe", "E4 very severe"))
  pc <- dm_erosion_class(dm_gaussian(c(4, 30), sd = 1))
  expect_equal(rowSums(pc), c(1, 1))
  expect_error(dm_erosion_class(1, breaks = c(5, 2), labels = c("a", "b", "c")), "increasing")
})
