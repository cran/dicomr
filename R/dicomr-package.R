#' dicomr: read, write and decode DICOM files natively in R
#'
#' A native R take on 'pydicom'. Start with [read_dicom()], access elements
#' with `ds$Keyword`, `ds[["Keyword"]]` or `ds[["0010,0010"]]`, decode images
#' with [dcm_pixel_array()], and write files with [write_dicom()].
#' For echo and other ultrasound, see [dcm_us_regions()] and [dcm_frame_times()].
#'
#' @useDynLib dicomr, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @keywords internal
"_PACKAGE"
