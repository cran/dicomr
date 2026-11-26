# Ultrasound and cine helpers --------------------------------------------------

.us_units <- c("0" = "none", "1" = "percent", "2" = "dB", "3" = "cm", "4" = "seconds",
               "5" = "hertz", "6" = "dB/seconds", "7" = "cm/sec", "8" = "cm2",
               "9" = "cm2/sec", "10" = "cm3", "11" = "cm3/sec", "12" = "degrees")

.us_spatial <- c("0" = "none", "1" = "2D", "2" = "M-mode", "3" = "spectral",
                 "4" = "waveform", "5" = "graphics")

.us_datatype <- c("0" = "none", "1" = "tissue", "2" = "color flow",
                  "3" = "PW spectral Doppler", "4" = "CW spectral Doppler",
                  "5" = "Doppler mean trace", "6" = "Doppler mode trace",
                  "7" = "Doppler max trace", "8" = "volume trace",
                  "10" = "ECG trace", "11" = "pulse trace", "12" = "phonocardiogram trace",
                  "13" = "gray bar", "14" = "color bar", "15" = "integrated backscatter",
                  "16" = "area trace", "17" = "d(area)/dt", "18" = "other physiological input")

.lab <- function(map, x) {
  out <- unname(map[as.character(x)])
  out[is.na(out) & !is.na(x)] <- "unknown"
  out
}

#' Ultrasound region calibration
#'
#' Parses the Sequence of Ultrasound Regions (0018,6011) into a data frame,
#' one row per region. Each region says where it sits in the image (in
#' 0-based pixel coordinates, as in the standard), what it shows (2D tissue,
#' colour flow, spectral Doppler, M-mode, ECG trace...), and how pixels map to
#' physical units: `delta_x` and `delta_y` are the physical size of one pixel
#' in `units_x` and `units_y` (for 2D echo usually cm per pixel; for spectral
#' Doppler usually seconds per pixel across and cm/sec per pixel up).
#'
#' @param ds A `dicom_dataset`.
#' @return A data frame (zero rows if the file has no regions).
#' @examples
#' echo <- read_dicom(system.file("extdata", "echo_synthetic.dcm", package = "dicomr"))
#' reg <- dcm_us_regions(echo)
#' reg[, c("region", "data_type_label", "delta_x", "units_x_label")]
#' @export
dcm_us_regions <- function(ds) {
  seq <- dcm_get(ds, "SequenceOfUltrasoundRegions")
  if (is.null(seq) || length(seq) == 0) {
    return(data.frame(region = integer(0), spatial_format = integer(0), data_type = integer(0),
                      delta_x = numeric(0), delta_y = numeric(0)))
  }
  num <- function(item, kw) {
    v <- dcm_get(item, kw)
    if (is.null(v) || length(v) == 0) NA_real_ else as.numeric(v[1])
  }
  rows <- lapply(seq_along(seq), function(i) {
    it <- seq[[i]]
    data.frame(
      region = i,
      spatial_format = num(it, "RegionSpatialFormat"),
      data_type = num(it, "RegionDataType"),
      flags = num(it, "RegionFlags"),
      x0 = num(it, "RegionLocationMinX0"),
      y0 = num(it, "RegionLocationMinY0"),
      x1 = num(it, "RegionLocationMaxX1"),
      y1 = num(it, "RegionLocationMaxY1"),
      ref_pixel_x0 = num(it, "ReferencePixelX0"),
      ref_pixel_y0 = num(it, "ReferencePixelY0"),
      units_x = num(it, "PhysicalUnitsXDirection"),
      units_y = num(it, "PhysicalUnitsYDirection"),
      ref_physical_x = num(it, "ReferencePixelPhysicalValueX"),
      ref_physical_y = num(it, "ReferencePixelPhysicalValueY"),
      delta_x = num(it, "PhysicalDeltaX"),
      delta_y = num(it, "PhysicalDeltaY"),
      transducer_frequency = num(it, "TransducerFrequency"),
      pulse_repetition_frequency = num(it, "PulseRepetitionFrequency"),
      doppler_correction_angle = num(it, "DopplerCorrectionAngle"),
      steering_angle = num(it, "SteeringAngle"),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$spatial_format_label <- .lab(.us_spatial, out$spatial_format)
  out$data_type_label <- .lab(.us_datatype, out$data_type)
  out$units_x_label <- .lab(.us_units, out$units_x)
  out$units_y_label <- .lab(.us_units, out$units_y)
  out
}

.pick_region <- function(reg, x0, y0, region) {
  if (!is.null(region)) {
    if (!region %in% reg$region) stop("No region ", region, call. = FALSE)
    return(reg[reg$region == region, , drop = FALSE])
  }
  inside <- which(x0 >= reg$x0 & x0 <= reg$x1 & y0 >= reg$y0 & y0 <= reg$y1)
  if (length(inside) == 0) return(NULL)
  # prefer calibrated regions, then the smallest one containing the point
  cal <- inside[!is.na(reg$delta_x[inside]) & reg$units_x[inside] != 0]
  if (length(cal)) inside <- cal
  area <- (reg$x1[inside] - reg$x0[inside]) * (reg$y1[inside] - reg$y0[inside])
  reg[inside[which.min(area)], , drop = FALSE]
}

#' Convert pixel positions to physical coordinates
#'
#' Uses the ultrasound region that contains each point (or the one you name)
#' to turn pixel positions into physical values:
#' `physical = reference_physical + (pixel - (region_origin + reference_pixel)) * delta`.
#'
#' @param ds A `dicom_dataset`, or a data frame from [dcm_us_regions()].
#' @param row,col Pixel positions as R array indices (1-based), matching
#'   `arr[row, col]` on a frame from [dcm_pixel_array()].
#' @param region Optional region number to use for every point.
#' @return A data frame with `row`, `col`, `region`, `x`, `y`, `units_x` and
#'   `units_y`; `x` runs across the image and `y` down it.
#' @examples
#' echo <- read_dicom(system.file("extdata", "echo_synthetic.dcm", package = "dicomr"))
#' dcm_us_pixel_to_physical(echo, row = c(1, 11), col = c(33, 33))
#' @export
dcm_us_pixel_to_physical <- function(ds, row, col, region = NULL) {
  reg <- if (is.data.frame(ds)) ds else dcm_us_regions(ds)
  if (nrow(reg) == 0) stop("No ultrasound regions in this data set", call. = FALSE)
  n <- max(length(row), length(col))
  row <- rep_len(row, n); col <- rep_len(col, n)
  out <- data.frame(row = row, col = col, region = NA_integer_, x = NA_real_, y = NA_real_,
                    units_x = NA_character_, units_y = NA_character_, stringsAsFactors = FALSE)
  for (i in seq_len(n)) {
    px <- col[i] - 1; py <- row[i] - 1
    r <- .pick_region(reg, px, py, region)
    if (is.null(r)) next
    ox <- r$x0 + ifelse(is.na(r$ref_pixel_x0), 0, r$ref_pixel_x0)
    oy <- r$y0 + ifelse(is.na(r$ref_pixel_y0), 0, r$ref_pixel_y0)
    rx <- ifelse(is.na(r$ref_physical_x), 0, r$ref_physical_x)
    ry <- ifelse(is.na(r$ref_physical_y), 0, r$ref_physical_y)
    out$region[i] <- r$region
    out$x[i] <- rx + (px - ox) * r$delta_x
    out$y[i] <- ry + (py - oy) * r$delta_y
    out$units_x[i] <- r$units_x_label
    out$units_y[i] <- r$units_y_label
  }
  out
}

#' Distance between two points in physical units
#'
#' For 2D regions where both axes share units (cm in echo), returns the
#' straight-line distance, for example a caliper measurement of an LV
#' dimension. Both points must fall in the same region.
#'
#' @param ds A `dicom_dataset` or a data frame from [dcm_us_regions()].
#' @param p1,p2 Points as `c(row, col)`, 1-based array indices.
#' @param region Optional region number.
#' @return A number with a `units` attribute.
#' @examples
#' echo <- read_dicom(system.file("extdata", "echo_synthetic.dcm", package = "dicomr"))
#' dcm_us_distance(echo, p1 = c(row = 5, col = 10), p2 = c(row = 25, col = 40))
#' @export
dcm_us_distance <- function(ds, p1, p2, region = NULL) {
  pts <- dcm_us_pixel_to_physical(ds, c(p1[1], p2[1]), c(p1[2], p2[2]), region)
  if (anyNA(pts$region) || pts$region[1] != pts$region[2])
    stop("Both points must lie in the same calibrated region", call. = FALSE)
  if (pts$units_x[1] != pts$units_y[1])
    stop("The region's axes have different units (", pts$units_x[1], ", ", pts$units_y[1],
         "); use dcm_us_pixel_to_physical() instead", call. = FALSE)
  structure(sqrt(diff(pts$x)^2 + diff(pts$y)^2), units = pts$units_x[1])
}

#' Frame times for a cine loop
#'
#' Returns the time of each frame in milliseconds from the first frame, from
#' Frame Time Vector (0018,1065) if present, otherwise Frame Time (0018,1063),
#' otherwise Cine Rate (0018,0040) or Recommended Display Frame Rate.
#'
#' @param ds A `dicom_dataset`.
#' @return A numeric vector with one time per frame (ms), or `NULL` if the
#'   file has no timing information. The frame rate (frames per second) is
#'   in the `frame_rate` attribute.
#' @examples
#' echo <- read_dicom(system.file("extdata", "echo_synthetic.dcm", package = "dicomr"))
#' dcm_frame_times(echo)
#' @export
dcm_frame_times <- function(ds) {
  nf <- dcm_pixel_info(ds)$frames
  ftv <- dcm_get(ds, "FrameTimeVector")
  ft <- dcm_get(ds, "FrameTime")
  rate <- dcm_get(ds, "CineRate")
  if (is.null(rate)) rate <- dcm_get(ds, "RecommendedDisplayFrameRate")
  if (!is.null(ftv) && length(ftv) == nf && is.numeric(ftv)) {
    times <- cumsum(as.numeric(ftv))
  } else if (!is.null(ft) && is.numeric(ft) && ft[1] > 0) {
    times <- (seq_len(nf) - 1) * as.numeric(ft[1])
  } else if (!is.null(rate) && is.numeric(rate) && rate[1] > 0) {
    times <- (seq_len(nf) - 1) * 1000 / as.numeric(rate[1])
  } else {
    return(NULL)
  }
  step <- if (nf > 1) mean(diff(times)) else NA_real_
  structure(times, frame_rate = if (is.na(step) || step <= 0) NA_real_ else 1000 / step)
}
