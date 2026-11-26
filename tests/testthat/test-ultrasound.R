test_that("ultrasound regions parse into a calibrated table", {
  reg <- dcm_us_regions(make_us_dataset())
  expect_equal(nrow(reg), 2)
  expect_equal(reg$data_type_label, c("tissue", "PW spectral Doppler"))
  expect_equal(reg$units_x_label, c("cm", "seconds"))
  expect_equal(reg$units_y_label, c("cm", "cm/sec"))
  expect_equal(reg$delta_x, c(0.05, 0.01))
})

test_that("pixel positions convert to physical units and distances", {
  ds <- make_us_dataset()
  p <- dcm_us_pixel_to_physical(ds, row = c(4, 14), col = c(3, 3))
  expect_equal(p$region, c(1L, 1L))
  expect_equal(diff(p$y), 0.5)          # 10 rows x 0.05 cm
  d <- dcm_us_distance(ds, c(4, 3), c(10, 11))
  expect_equal(as.numeric(d), sqrt((6 * 0.05)^2 + (8 * 0.05)^2))
  expect_equal(attr(d, "units"), "cm")
  dop <- dcm_us_pixel_to_physical(ds, row = 20, col = 11)
  expect_equal(dop$region, 2L)
  expect_equal(dop$x, 0.1)
  expect_equal(dop$units_y, "cm/sec")
  expect_error(dcm_us_distance(ds, c(4, 3), c(20, 11)), "same calibrated region")
})

test_that("frame times come from Frame Time", {
  ft <- dcm_frame_times(make_us_dataset(frames = 4))
  expect_equal(as.numeric(ft), c(0, 33.3, 66.6, 99.9))
  expect_equal(attr(ft, "frame_rate"), 1000 / 33.3)
})
