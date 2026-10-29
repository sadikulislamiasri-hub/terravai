test_that("curve number tables and lookups", {
  tab <- dm_cn_table()
  expect_equal(tab$B[tab$cover == "row_crops_straight_good"], 78)
  expect_equal(tab$D[tab$cover == "fallow_bare"], 94)
  expect_true(all(tab$A <= tab$B & tab$B <= tab$C & tab$C <= tab$D))
  expect_equal(dm_cn_lookup(c("woods_good", "meadow"), c("C", "A")), c(70, 30))
  expect_error(dm_cn_lookup("rice", "A"), "Unknown cover")
  expect_error(dm_cn_lookup("meadow", "E"), "hsg")
})

test_that("USDA texture classes and hydrologic groups", {
  cls <- dm_texture_class(c(92, 82, 65, 40, 20, 5, 60, 30, 10, 50, 5, 20),
                          c(4, 12, 25, 40, 65, 90, 12, 37, 58, 5, 50, 20),
                          c(4, 6, 10, 20, 15, 5, 28, 33, 32, 45, 45, 60))
  expect_equal(cls, c("sand", "loamy sand", "sandy loam", "loam", "silt loam", "silt",
                      "sandy clay loam", "clay loam", "silty clay loam", "sandy clay",
                      "silty clay", "clay"))
  expect_equal(dm_hsg_texture(texture = c("sand", "loam", "sandy clay loam", "clay")),
               c("A", "B", "C", "D"))
  expect_error(dm_texture_class(120, 0, 0), "between")
})

test_that("support-practice factors follow Wischmeier and Smith", {
  expect_equal(dm_p_factor("contouring", c(5, 10, 15, 18, 22)), c(0.5, 0.6, 0.7, 0.8, 0.9))
  expect_equal(dm_p_factor("strip_cropping", c(5, 12, 20)), c(0.25, 0.30, 0.40))
  expect_equal(dm_p_factor("none", 30), 1)
  expect_true(is.na(dm_p_factor("contouring", 30)))
})

test_that("raster erosion pipeline works on projected and geographic grids", {
  skip_if_not_installed("terra")
  mk_dem <- function(crs, ext) {
    dem <- terra::rast(nrows = 25, ncols = 25, ext = ext, crs = crs)
    terra::values(dem) <- as.numeric(outer(1:25, 1:25, function(i, j) 500 + 6 * i + 3 * j))
    dem
  }
  for (cfg in list(list("EPSG:4326", terra::ext(77, 77.025, 30, 30.025)),
                   list("EPSG:32643", terra::ext(500000, 500750, 3300000, 3300750)))) {
    dem <- mk_dem(cfg[[1]], cfg[[2]])
    sand <- terra::init(dem, 300); silt <- terra::init(dem, 450); clay <- terra::init(dem, 250)
    soc <- terra::init(dem, 120)
    out <- dm_erosion_raster(dem, R = 800, ndvi = terra::init(dem, 0.4), sand = sand,
                             silt = silt, clay = clay, soc = soc, texture_unit = "g_kg",
                             soc_unit = "dg_kg")
    expect_s4_class(out, "SpatRaster")
    expect_equal(names(out), c("slope", "flow_acc", "LS", "C", "K", "P", "A"))
    v <- terra::values(out)
    expect_true(all(v[, "A"] >= 0, na.rm = TRUE))
    expect_equal(max(v[, "flow_acc"], na.rm = TRUE), 625)
    expect_equal(unname(v[1, "K"]), dm_k_epic(30, 45, 25, 1.2), tolerance = 1e-6)
    to <- attr(out, "to")
    expect_equal(length(to), 625)
  }
  expect_error(dm_erosion_raster(dem, R = 800), "ndvi")
})
