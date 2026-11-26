# Pixel data ------------------------------------------------------------------

.ts_rle <- "1.2.840.10008.1.2.5"
.ts_jpeg8 <- c("1.2.840.10008.1.2.4.50", "1.2.840.10008.1.2.4.51")

.int1 <- function(ds, kw, default) {
  v <- dcm_get(ds, kw)
  if (is.null(v) || length(v) == 0 || is.na(v[1])) return(default)
  as.integer(v[1])
}

#' Image geometry and pixel description
#'
#' @param ds A `dicom_dataset`.
#' @return A list with `rows`, `columns`, `samples_per_pixel`, `frames`,
#'   `bits_allocated`, `bits_stored`, `pixel_representation`,
#'   `planar_configuration`, `photometric_interpretation` and
#'   `transfer_syntax`.
#' @examples
#' echo <- read_dicom(system.file("extdata", "echo_synthetic.dcm", package = "dicomr"))
#' str(dcm_pixel_info(echo))
#' @export
dcm_pixel_info <- function(ds) {
  pi <- dcm_get(ds, "PhotometricInterpretation")
  nf <- dcm_get(ds, "NumberOfFrames")
  nf <- if (is.null(nf) || length(nf) == 0 || is.na(suppressWarnings(as.integer(nf[1])))) 1L else max(1L, as.integer(nf[1]))
  ba <- .int1(ds, "BitsAllocated", NA_integer_)
  list(
    rows = .int1(ds, "Rows", NA_integer_),
    columns = .int1(ds, "Columns", NA_integer_),
    samples_per_pixel = .int1(ds, "SamplesPerPixel", 1L),
    frames = nf,
    bits_allocated = ba,
    bits_stored = .int1(ds, "BitsStored", ba),
    pixel_representation = .int1(ds, "PixelRepresentation", 0L),
    planar_configuration = .int1(ds, "PlanarConfiguration", 0L),
    photometric_interpretation = if (is.null(pi)) NA_character_ else trimws(pi[1]),
    transfer_syntax = as.character(dcm_transfer_syntax(ds))
  )
}

# Split encapsulated fragments into one raw vector per frame (0-based index).
.encapsulated_frames <- function(ds, enc, nframes, ts) {
  frags <- enc$fragments
  nfrag <- length(frags)
  if (nfrag == 0) stop("Encapsulated pixel data has no fragments", call. = FALSE)
  if (nframes == 1) return(list(do.call(c, frags)))

  offsets <- enc$offsets
  if (length(offsets) == 0) {
    eot <- dcm_get(ds, "ExtendedOffsetTable")
    if (is.raw(eot) && length(eot) >= 8 * nframes) {
      v <- matrix(as.numeric(eot), nrow = 8)
      offsets <- colSums(v * 256^(0:7))
    }
  }
  if (length(offsets) >= nframes) {
    starts <- c(0, cumsum(8 + lengths(frags)))[seq_len(nfrag)]
    bounds <- c(offsets[seq_len(nframes)], Inf)
    out <- vector("list", nframes)
    for (f in seq_len(nframes)) {
      idx <- which(starts >= bounds[f] & starts < bounds[f + 1])
      out[[f]] <- do.call(c, frags[idx])
    }
    return(out)
  }
  if (nfrag == nframes) return(frags)
  if (ts %in% .ts_jpeg8 || grepl("^1\\.2\\.840\\.10008\\.1\\.2\\.4\\.", ts)) {
    soi <- .fragment_starts_soi(frags)
    soi[1] <- TRUE
    grp <- cumsum(soi)
    if (max(grp) == nframes) return(lapply(split(frags, grp), function(x) do.call(c, x)))
  }
  stop(sprintf("Can't match %d fragments to %d frames without an offset table", nfrag, nframes),
       call. = FALSE)
}

#' Decode pixel data into an R array
#'
#' Returns pixel data with the same axis order as pydicom's `pixel_array`:
#' `frames x rows x columns x samples`. With `drop = TRUE` (the default) the
#' frame axis is dropped for single-frame images and the sample axis is
#' dropped for monochrome images, so a single-frame greyscale image is a
#' `rows x columns` matrix and a 30-frame colour echo clip is a
#' `30 x rows x columns x 3` array. `arr[f, r, c, s]` in R is
#' `arr[f-1, r-1, c-1, s-1]` in Python.
#'
#' Supported encodings: uncompressed (implicit/explicit VR, little/big
#' endian, deflated), RLE Lossless, and 8-bit JPEG Baseline / Extended
#' through libjpeg-turbo. YBR_FULL and YBR_FULL_422 data are converted to RGB
#' when `as_rgb = TRUE`, using the same float formula as pydicom.
#'
#' The values are the stored values. Apply [dcm_apply_rescale()] for
#' Hounsfield units and similar.
#'
#' @param ds A `dicom_dataset` from [read_dicom()].
#' @param frames Optional 1-based frame numbers to decode, for example
#'   `1:10`. Only these frames are decoded, which saves time and memory on
#'   long cine loops.
#' @param as_rgb Convert YBR colour data to RGB.
#' @param drop Drop length-one frame and sample axes as pydicom does.
#' @param correct_unused_bits Clear (or sign-extend) the bits above
#'   Bits Stored. pydicom's default is `FALSE`.
#' @param jpeg_mode How to decode colour JPEG. `"pydicom"` (the default)
#'   reproduces pydicom with its Pillow plugin exactly: libjpeg-turbo decodes
#'   without a colour transform, then pydicom's YBR to RGB formula is applied.
#'   Use it when results must match a Python pipeline, such as a model
#'   trained on pydicom output. `"accurate"` lets libjpeg-turbo do the
#'   colour conversion in integer arithmetic; values can differ from
#'   `"pydicom"` by a grey level or two.
#' @return An integer array (double for 32-bit unsigned or float data), with
#'   attribute `photometric_interpretation` giving the colour space of the
#'   returned values.
#' @examples
#' echo <- read_dicom(system.file("extdata", "echo_synthetic.dcm", package = "dicomr"))
#' clip <- dcm_pixel_array(echo)               # frames x rows x cols x 3
#' dim(clip)
#' first <- dcm_pixel_array(echo, frames = 1)  # rows x cols x 3
#' dim(first)
#' ct <- read_dicom(system.file("extdata", "ct_synthetic.dcm", package = "dicomr"))
#' dim(dcm_pixel_array(ct))                    # rows x cols
#' @export
dcm_pixel_array <- function(ds, frames = NULL, as_rgb = TRUE, drop = TRUE,
                            correct_unused_bits = FALSE,
                            jpeg_mode = c("pydicom", "accurate")) {
  jpeg_mode <- match.arg(jpeg_mode)
  info <- dcm_pixel_info(ds)
  pe <- dcm_element(ds, "PixelData")
  is_float <- FALSE
  if (is.null(pe)) {
    pe <- dcm_element(ds, "FloatPixelData")
    if (is.null(pe)) pe <- dcm_element(ds, "DoubleFloatPixelData")
    if (is.null(pe)) stop("The data set has no pixel data", call. = FALSE)
    is_float <- TRUE
  }
  if (is.na(info$rows) || is.na(info$columns) || is.na(info$bits_allocated))
    stop("Rows, Columns and Bits Allocated are required to decode pixel data", call. = FALSE)

  nf <- info$frames
  if (is.null(frames)) frames <- seq_len(nf)
  frames <- as.integer(frames)
  if (any(is.na(frames)) || any(frames < 1) || any(frames > nf))
    stop(sprintf("frames must be between 1 and %d", nf), call. = FALSE)
  f0 <- frames - 1L

  ts <- info$transfer_syntax
  pi <- info$photometric_interpretation
  spp <- info$samples_per_pixel
  enc <- attr(ds, "encoding")
  little <- if (is.null(enc)) TRUE else isTRUE(enc$little_endian)
  value <- pe$value

  if (inherits(value, "dicom_encapsulated")) {
    fr <- .encapsulated_frames(ds, value, nf, ts)[frames]
    if (ts == .ts_rle) {
      arr <- .decode_rle_frames(fr, info$rows, info$columns, spp, info$bits_allocated,
                                info$bits_stored, info$pixel_representation, correct_unused_bits)
    } else if (ts %in% .ts_jpeg8) {
      arr <- .decode_jpeg_frames(fr, info$rows, info$columns, spp, pi, ts == .ts_jpeg8[1],
                                 as_rgb, jpeg_mode == "accurate")
      pi <- attr(arr, "photometric_interpretation")
    } else {
      nm <- .ts_names[ts]
      stop(sprintf("Decoding '%s' (%s) is not supported yet", if (is.na(nm)) ts else nm, ts),
           call. = FALSE)
    }
  } else {
    if (!is.raw(value)) stop("Pixel data is not a byte array", call. = FALSE)
    if (identical(pi, "YBR_FULL_422") && info$bits_allocated == 8 && spp == 3) {
      arr <- .decode_ybr422(value, info$rows, info$columns, f0)
      pi <- "YBR_FULL"
    } else {
      swap8 <- !little && info$bits_allocated == 8 && identical(pe$vr, "OW")
      arr <- .decode_native(value, info$rows, info$columns, spp, info$bits_allocated,
                            info$bits_stored, info$pixel_representation,
                            info$planar_configuration, nf, f0, little, swap8,
                            is_float, correct_unused_bits)
    }
  }

  if (as_rgb && spp == 3 && pi %in% c("YBR_FULL", "YBR_FULL_422") && is.integer(arr)) {
    arr <- .ybr_to_rgb(arr)
    pi <- "RGB"
  }

  if (drop) {
    d <- dim(arr)
    keep <- c(length(frames) > 1, TRUE, TRUE, spp > 1)
    dim(arr) <- d[keep]
  }
  attr(arr, "photometric_interpretation") <- pi
  arr
}

#' Apply the modality rescale (slope and intercept)
#'
#' Converts stored values to output units such as Hounsfield units, using
#' Rescale Slope (0028,1053) and Rescale Intercept (0028,1052).
#'
#' @param arr An array from [dcm_pixel_array()].
#' @param ds The data set it came from.
#' @return A double array.
#' @examples
#' ct <- read_dicom(system.file("extdata", "ct_synthetic.dcm", package = "dicomr"))
#' hu <- dcm_apply_rescale(dcm_pixel_array(ct), ct)
#' range(hu)
#' @export
dcm_apply_rescale <- function(arr, ds) {
  slope <- dcm_get(ds, "RescaleSlope", 1)
  intercept <- dcm_get(ds, "RescaleIntercept", 0)
  out <- arr * as.numeric(slope[1]) + as.numeric(intercept[1])
  attr(out, "photometric_interpretation") <- attr(arr, "photometric_interpretation")
  out
}

#' Show a frame
#'
#' A quick look at one frame with base graphics.
#'
#' @param ds A `dicom_dataset`.
#' @param frame Frame number (1-based).
#' @param ... Passed to [graphics::plot()] (colour) or [graphics::image()] (greyscale).
#' @return `NULL`, invisibly.
#' @examples
#' echo <- read_dicom(system.file("extdata", "echo_synthetic.dcm", package = "dicomr"))
#' dcm_show(echo, frame = 2)
#' ct <- read_dicom(system.file("extdata", "ct_synthetic.dcm", package = "dicomr"))
#' dcm_show(ct)
#' @export
dcm_show <- function(ds, frame = 1, ...) {
  a <- dcm_pixel_array(ds, frames = frame, drop = FALSE)
  d <- dim(a)
  op <- graphics::par(mar = c(0, 0, 0, 0))
  on.exit(graphics::par(op))
  if (d[4] == 3) {
    img <- a[1, , , ] / 255
    img[img < 0] <- 0; img[img > 1] <- 1
    r <- grDevices::as.raster(img)
    graphics::plot(r, ...)
  } else {
    m <- a[1, , , 1]
    if (identical(attr(a, "photometric_interpretation"), "MONOCHROME1")) m <- max(m) - m
    graphics::image(t(m[nrow(m):1, , drop = FALSE]), col = grDevices::gray.colors(256, 0, 1),
                    axes = FALSE, useRaster = TRUE, ...)
  }
  invisible(NULL)
}

#' Store an R array as uncompressed pixel data
#'
#' The inverse of [dcm_pixel_array()]: writes `arr` into Pixel Data and sets
#' Rows, Columns, Samples per Pixel, Number of Frames, Bits Allocated/Stored,
#' High Bit, Pixel Representation, Photometric Interpretation and Planar
#' Configuration to match. The data set's transfer syntax becomes Explicit
#' VR Little Endian if it was compressed.
#'
#' @param ds A `dicom_dataset`.
#' @param arr An integer array: `rows x cols` (one greyscale frame),
#'   `frames x rows x cols` (greyscale cine), `rows x cols x 3` (one colour
#'   frame) or `frames x rows x cols x 3`.
#' @param photometric_interpretation `"MONOCHROME2"`, `"MONOCHROME1"` or
#'   `"RGB"`. Guessed from the array when missing.
#' @param bits_allocated 8, 16 or 32. Chosen from the value range when missing.
#' @return The modified data set.
#' @examples
#' img <- matrix(as.integer(seq(0, 4095, length.out = 64 * 64)), 64, 64)
#' ds <- dicom_dataset(Modality = "OT")
#' ds <- dcm_set_pixel_data(ds, img)
#' ds$BitsAllocated
#' identical(as.vector(dcm_pixel_array(ds)), as.vector(img))
#' @export
dcm_set_pixel_data <- function(ds, arr, photometric_interpretation = NULL, bits_allocated = NULL) {
  d <- dim(arr)
  if (is.null(d)) stop("arr must be a matrix or array", call. = FALSE)
  if (is.null(photometric_interpretation)) {
    photometric_interpretation <- if (length(d) == 4 || (length(d) == 3 && d[3] == 3)) "RGB" else "MONOCHROME2"
  }
  spp <- if (photometric_interpretation == "RGB") 3L else 1L
  dims <- switch(as.character(length(d)),
    "2" = c(1L, d, 1L),
    "3" = if (spp == 3) c(1L, d) else c(d, 1L),
    "4" = d,
    stop("arr has too many dimensions", call. = FALSE))
  if (dims[4] != spp) stop("The last axis must have 3 samples for RGB", call. = FALSE)
  v <- as.numeric(arr)
  if (anyNA(v)) stop("arr contains NA", call. = FALSE)
  signed <- any(v < 0)
  if (is.null(bits_allocated)) {
    rng <- range(v)
    bits_allocated <- if (!signed && rng[2] <= 255) 8L
                      else if ((signed && rng[1] >= -32768 && rng[2] <= 32767) || (!signed && rng[2] <= 65535)) 16L
                      else 32L
  }
  dim(arr) <- dims
  # DICOM order: samples fastest, then columns, rows, frames
  bytes <- writeBin(as.integer(aperm(arr, c(4, 3, 2, 1))), raw(), size = bits_allocated / 8,
                    endian = "little")
  if (length(bytes) %% 2) bytes <- c(bytes, as.raw(0))
  ds <- dcm_set(ds, "SamplesPerPixel", spp)
  ds <- dcm_set(ds, "PhotometricInterpretation", photometric_interpretation)
  if (spp > 1) ds <- dcm_set(ds, "PlanarConfiguration", 0L) else ds <- dcm_set(ds, "PlanarConfiguration", NULL)
  ds <- dcm_set(ds, "NumberOfFrames", if (dims[1] > 1) as.integer(dims[1]) else NULL)
  ds <- dcm_set(ds, "Rows", dims[2])
  ds <- dcm_set(ds, "Columns", dims[3])
  ds <- dcm_set(ds, "BitsAllocated", bits_allocated)
  ds <- dcm_set(ds, "BitsStored", bits_allocated)
  ds <- dcm_set(ds, "HighBit", bits_allocated - 1L)
  ds <- dcm_set(ds, "PixelRepresentation", as.integer(signed))
  ds <- dcm_set(ds, "PixelData", bytes, vr = if (bits_allocated > 8) "OW" else "OB")
  enc <- attr(ds, "encoding")
  if (is.null(enc) || isTRUE(enc$encapsulated)) {
    attr(ds, "encoding") <- list(transfer_syntax = "1.2.840.10008.1.2.1", little_endian = TRUE,
                                 implicit_vr = FALSE, encapsulated = FALSE, deflated = FALSE,
                                 has_meta = TRUE)
  } else if (!isTRUE(enc$little_endian)) {
    # bytes above are little endian; keep the data set consistent
    enc$little_endian <- TRUE
    enc$transfer_syntax <- "1.2.840.10008.1.2.1"
    attr(ds, "encoding") <- enc
  }
  ds
}
