uncompressed <- c("1.2.840.10008.1.2", "1.2.840.10008.1.2.1", "1.2.840.10008.1.2.2", "1.2.840.10008.1.2.1.99")

test_that("a data set survives writing and reading in every uncompressed syntax", {
  ds <- make_us_dataset()
  ds$PatientName <- "M\u00fcller^J\u00fcrgen"   # non-ASCII, written as Latin-1
  ds$SpecificCharacterSet <- "ISO_IR 100"
  ds$WindowCenter <- c(40.5, 400)
  ds$FrameIncrementPointer <- "FrameTime"
  arr <- dcm_pixel_array(ds)
  for (ts in uncompressed) {
    f <- tempfile(fileext = ".dcm")
    write_dicom(ds, f, transfer_syntax = ts)
    back <- read_dicom(f)
    expect_equal(as.character(dcm_transfer_syntax(back)), ts)
    expect_equal(enc2utf8(back$PatientName), "M\u00fcller^J\u00fcrgen")
    expect_equal(back$WindowCenter, c(40.5, 400))
    expect_equal(back$FrameIncrementPointer, "00181063")
    expect_equal(back$ImageType, c("ORIGINAL", "PRIMARY"))
    expect_equal(length(back$SequenceOfUltrasoundRegions), 2)
    expect_equal(back$SequenceOfUltrasoundRegions[[2]]$PhysicalDeltaY, -2.5)
    expect_identical(unclass(dcm_pixel_array(back)), unclass(arr))
    unlink(f)
  }
})

test_that("pixel data layout matches the frames x rows x cols x samples convention", {
  ds <- make_us_dataset(frames = 3, rows = 20, cols = 30)
  a <- dcm_pixel_array(ds)
  expect_equal(dim(a), c(3L, 20L, 30L, 3L))
  expect_equal(a[2, 1, 1, 1], 2L)
  expect_equal(a[1, 5, 7, 2], 200L)
  expect_equal(a[1, 4, 9, 3], 13L)
  one <- dcm_pixel_array(ds, frames = 2)
  expect_equal(dim(one), c(20L, 30L, 3L))
  expect_identical(as.vector(one), as.vector(a[2, , , ]))
  expect_error(dcm_pixel_array(ds, frames = 4), "between 1 and 3")
})

test_that("greyscale 16-bit signed data round trips and rescales", {
  m <- matrix(as.integer(seq(-1000, 1000, length.out = 64)), 8, 8)
  ds <- dcm_set_pixel_data(dicom_dataset(Modality = "CT", RescaleSlope = 1, RescaleIntercept = -1024), m)
  expect_equal(ds$BitsAllocated, 16L)
  expect_equal(ds$PixelRepresentation, 1L)
  f <- tempfile(fileext = ".dcm")
  write_dicom(ds, f)
  back <- read_dicom(f)
  a <- dcm_pixel_array(back)
  expect_identical(as.vector(a), as.vector(m))
  expect_equal(as.vector(dcm_apply_rescale(a, back)), as.vector(m) - 1024)
})

# Minimal PackBits encoder (literal runs only) to build RLE test data
rle_frame <- function(segments) {
  enc <- lapply(segments, function(s) {
    out <- raw(0)
    for (st in seq(1, length(s), by = 128)) {
      chunk <- s[st:min(length(s), st + 127)]
      out <- c(out, as.raw(length(chunk) - 1), chunk)
    }
    if (length(out) %% 2) out <- c(out, as.raw(0x80))
    out
  })
  offs <- 64 + c(0, cumsum(lengths(enc)))[seq_along(enc)]
  hdr <- writeBin(as.integer(c(length(enc), offs, rep(0, 15 - length(enc)))), raw(), size = 4, endian = "little")
  c(hdr, unlist(enc))
}

test_that("RLE Lossless frames decode", {
  rows <- 4; cols <- 5
  v <- matrix(as.integer(1:20 * 300), rows, cols, byrow = FALSE)
  px <- as.vector(t(v))                      # row-major pixel order
  hi <- as.raw(px %/% 256); lo <- as.raw(px %% 256)
  frame <- rle_frame(list(hi, lo))
  ds <- dicom_dataset(Rows = rows, Columns = cols, BitsAllocated = 16L, BitsStored = 16L,
                      HighBit = 15L, PixelRepresentation = 0L, SamplesPerPixel = 1L,
                      PhotometricInterpretation = "MONOCHROME2")
  enc <- structure(list(offsets = numeric(0), fragments = list(frame)), class = "dicom_encapsulated")
  ds <- dcm_set(ds, "PixelData", enc, vr = "OB")
  attr(ds, "encoding") <- list(transfer_syntax = "1.2.840.10008.1.2.5", little_endian = TRUE,
                               implicit_vr = FALSE, encapsulated = TRUE, deflated = FALSE, has_meta = TRUE)
  expect_identical(as.vector(dcm_pixel_array(ds)), as.vector(v))
  f <- tempfile(fileext = ".dcm")
  write_dicom(ds, f)
  back <- read_dicom(f)
  expect_identical(as.vector(dcm_pixel_array(back)), as.vector(v))
})
