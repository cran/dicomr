# Cross-checks against pydicom on its bundled sample files. Skipped unless
# Python 3 with pydicom, numpy and Pillow is available.
pydicom_dir <- function() {
  py <- Sys.which("python3")
  if (!nzchar(py)) return("")
  out <- suppressWarnings(tryCatch(system2(py, c("-c", shQuote(paste(
    "import os, pydicom.data, PIL, numpy;",
    "print(os.path.dirname(pydicom.data.get_testdata_file('CT_small.dcm')))"))),
    stdout = TRUE, stderr = FALSE), error = function(e) ""))
  if (length(out) == 0) "" else out[length(out)]
}

test_that("tags and pixels match pydicom", {
  skip_on_cran()
  skip_if_not_installed("jsonlite")
  dir <- pydicom_dir()
  skip_if(!nzchar(dir) || !dir.exists(dir), "pydicom is not available")
  files <- c("examples_ybr_color.dcm", "examples_rgb_color.dcm", "examples_palette.dcm",
             "ExplVR_BigEnd.dcm", "CT_small.dcm", "MR_small_implicit.dcm", "MR_small_RLE.dcm",
             "MR_small_bigendian.dcm", "SC_rgb_rle_32bit_2frame.dcm", "SC_rgb_jpeg_dcmtk.dcm",
             "SC_rgb_dcmtk_+eb+cy+np.dcm", "SC_jpeg_no_color_transform.dcm",
             "SC_ybr_full_422_uncompressed.dcm", "image_dfl.dcm", "liver_1frame.dcm",
             "rtdose_rle.dcm", "test-SR.dcm", "waveform_ecg.dcm")
  for (f in files) {
    path <- file.path(dir, f)
    if (!file.exists(path)) next
    r <- dcm_compare_pydicom(path)
    expect_equal(nrow(r$issues), 0, label = paste(f, "tag differences"))
    if (r$pixels$status != "none") expect_equal(r$pixels$status, "equal", label = paste(f, "pixels"))
  }
})
