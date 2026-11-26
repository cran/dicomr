# Uses the synthetic example files shipped in inst/extdata, so these run
# everywhere, including on CRAN, and exercise the libjpeg decoder.

echo_truth <- function() {
  frames <- 4L; rows <- 48L; cols <- 64L
  truth <- array(0L, c(frames, rows, cols, 3L))
  for (f in seq_len(frames)) {
    truth[f, , , 1] <- outer(seq_len(rows), seq_len(cols), function(r, c) 40L + 3L * r + 10L * f)
    truth[f, , , 2] <- outer(seq_len(rows), seq_len(cols), function(r, c) 30L + 2L * c)
    truth[f, , , 3] <- 120L + 20L * f
  }
  truth
}

test_that("JPEG Baseline colour cine decodes close to the source image", {
  f <- system.file("extdata", "echo_synthetic.dcm", package = "dicomr")
  ds <- read_dicom(f)
  expect_equal(as.character(dcm_transfer_syntax(ds)), "1.2.840.10008.1.2.4.50")
  truth <- echo_truth()
  for (mode in c("pydicom", "accurate")) {
    a <- dcm_pixel_array(ds, jpeg_mode = mode)
    expect_equal(dim(a), c(4L, 48L, 64L, 3L))
    expect_equal(attr(a, "photometric_interpretation"), "RGB")
    d <- abs(a - truth)
    expect_lt(max(d), 12)
    expect_lt(mean(d), 2)
  }
  one <- dcm_pixel_array(ds, frames = 3)
  expect_equal(dim(one), c(48L, 64L, 3L))
  expect_identical(as.vector(one), as.vector(dcm_pixel_array(ds)[3, , , ]))
  raw_ybr <- dcm_pixel_array(ds, frames = 1, as_rgb = FALSE)
  expect_equal(attr(raw_ybr, "photometric_interpretation"), "YBR_FULL_422")
})

test_that("ultrasound calibration and frame timing are read", {
  ds <- read_dicom(system.file("extdata", "echo_synthetic.dcm", package = "dicomr"))
  reg <- dcm_us_regions(ds)
  expect_equal(nrow(reg), 2)
  expect_equal(reg$delta_x[1], 0.05)
  expect_equal(as.vector(dcm_us_distance(ds, c(1, 1), c(1, 21))), 1)
  expect_equal(as.vector(dcm_frame_times(ds)), c(0, 33.3, 66.6, 99.9))
})

test_that("the CT example rescales to Hounsfield units", {
  ct <- read_dicom(system.file("extdata", "ct_synthetic.dcm", package = "dicomr"))
  expect_equal(as.character(dcm_transfer_syntax(ct)), "1.2.840.10008.1.2")
  hu <- dcm_apply_rescale(dcm_pixel_array(ct), ct)
  expect_equal(dim(hu), c(32L, 32L))
  expect_equal(range(hu), c(-1000, 0))
})

test_that("corrupt JPEG data gives an R error, not a crash", {
  ds <- read_dicom(system.file("extdata", "echo_synthetic.dcm", package = "dicomr"))
  px <- ds$PixelData
  px$fragments[[1]] <- px$fragments[[1]][1:40]
  px$fragments[[2]] <- as.raw(c(0xFF, 0xD8, 0x00, 0x01, 0x02))
  ds <- dcm_set(ds, "PixelData", px, vr = "OB")
  expect_error(dcm_pixel_array(ds, frames = 2))
})
