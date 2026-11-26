# Cross-checking against pydicom ------------------------------------------------

.md5_raw <- function(b) {
  f <- tempfile()
  on.exit(unlink(f))
  writeBin(b, f)
  unname(tools::md5sum(f))
}

.same_value <- function(r, p) {
  # r: R value; p: value parsed from pydicom JSON (list or NULL)
  if (is.null(p)) return(is.null(r) || length(r) == 0 || (is.character(r) && all(r == "")))
  if (is.list(p) && !is.null(p$md5)) return(is.raw(r) && length(r) == p$len && .md5_raw(r) == p$md5)
  if (is.list(p) && !is.null(p$fragments)) {
    if (!inherits(r, "dicom_encapsulated")) return(FALSE)
    return(identical(vapply(r$fragments, .md5_raw, ""), unlist(p$fragments)))
  }
  if (is.null(r) || length(r) != length(p)) return(FALSE)
  for (i in seq_along(p)) {
    pv <- p[[i]]
    rv <- r[[i]]
    if (is.null(pv)) { if (!(is.na(rv) || identical(rv, ""))) return(FALSE); next }
    if (is.character(pv)) {
      if (is.numeric(rv)) { if (!isTRUE(abs(as.numeric(pv) - rv) < 1e-9)) return(FALSE); next }
      if (!identical(enc2utf8(as.character(rv)), enc2utf8(pv))) return(FALSE)
    } else {
      if (is.character(rv)) return(FALSE)
      if (is.na(rv) || !isTRUE(abs(pv - rv) <= 1e-9 * max(1, abs(pv)))) return(FALSE)
    }
  }
  TRUE
}

.compare_ds <- function(rds, pds, path = "") {
  issues <- list()
  add <- function(tag, issue) issues[[length(issues) + 1]] <<- data.frame(
    path = path, tag = dcm_tag_string(tag), issue = issue, stringsAsFactors = FALSE)
  rnames <- names(rds)
  pnames <- names(pds)
  for (t in setdiff(pnames, rnames)) add(t, "missing in dicomr")
  for (t in setdiff(rnames, pnames)) add(t, "only in dicomr")
  n <- 0L
  for (t in intersect(pnames, rnames)) {
    n <- n + 1L
    rel <- .subset2(rds, t)
    pel <- pds[[t]]
    if (!identical(rel$vr, pel$vr)) add(t, sprintf("VR %s vs pydicom %s", rel$vr, pel$vr))
    if (identical(pel$vr, "SQ") && identical(rel$vr, "SQ")) {
      if (length(rel$value) != length(pel$value)) { add(t, "different number of items"); next }
      for (i in seq_along(pel$value)) {
        sub <- .compare_ds(rel$value[[i]], pel$value[[i]], sprintf("%s%s[%d]/", path, dcm_tag_string(t), i))
        n <- n + sub$n
        issues <- c(issues, list(sub$issues))
      }
    } else if (!.same_value(rel$value, pel$value)) {
      add(t, "value differs")
    }
  }
  iss <- if (length(issues)) do.call(rbind, issues) else NULL
  if (is.null(iss)) iss <- data.frame(path = character(0), tag = character(0), issue = character(0))
  list(n = n, issues = iss)
}

#' Compare dicomr's reading of a file with pydicom's
#'
#' Runs pydicom (through Python) on the same file and compares every data
#' element, including nested sequence items, and the decoded pixel array.
#' Useful for checking your own files. Requires Python 3 with `pydicom` and
#' `numpy` installed, plus the 'jsonlite' package.
#'
#' @param file Path to a DICOM file.
#' @param python The Python executable.
#' @param pixels Also compare the decoded pixel arrays.
#' @param jpeg_plugin pydicom decoding plugin for JPEG Baseline/Extended data:
#'   `"pillow"` (the default, which `jpeg_mode = "pydicom"` reproduces) or
#'   `"pylibjpeg"`, or `"default"` to let pydicom choose.
#' @param ... Passed to [dcm_pixel_array()], for example `jpeg_mode`.
#' @return A list with `elements` (number compared), `issues` (a data frame
#'   of differences), and `pixels`: a list with `status` ("equal",
#'   "close", "different", "pydicom_failed", "dicomr_failed", "both_failed"
#'   or "none"), `max_abs_diff`, `fraction_equal` and any error messages.
#' @examples
#' # Needs Python 3 with pydicom and numpy, so it is not run here
#' \dontrun{
#' r <- dcm_compare_pydicom(system.file("extdata", "ct_synthetic.dcm", package = "dicomr"))
#' r$issues
#' r$pixels$status
#' }
#' @export
dcm_compare_pydicom <- function(file, python = "python3", pixels = TRUE,
                                jpeg_plugin = "pillow", ...) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Install 'jsonlite' to use this", call. = FALSE)
  script <- system.file("python", "pydicom_reference.py", package = "dicomr")
  out <- tempfile("pydicom_ref")
  on.exit(unlink(out, recursive = TRUE))
  status <- system2(python, c("-I", shQuote(script), shQuote(file), shQuote(out), jpeg_plugin), stdout = TRUE, stderr = TRUE)
  ref_file <- file.path(out, "tags.json")
  if (!file.exists(ref_file)) stop("pydicom could not read the file: ", paste(status, collapse = "\n"), call. = FALSE)
  ref <- jsonlite::fromJSON(ref_file, simplifyVector = FALSE)

  ds <- suppressWarnings(read_dicom(file, force = TRUE))
  meta <- .compare_ds(dcm_file_meta(ds), ref$meta, "file_meta/")
  main <- .compare_ds(ds, ref$dataset)
  res <- list(elements = meta$n + main$n, issues = rbind(meta$issues, main$issues))

  px <- list(status = "none", max_abs_diff = NA_real_, fraction_equal = NA_real_,
             pydicom_error = ref$pixel_error, dicomr_error = NULL)
  if (pixels && dcm_has(ds, "PixelData")) {
    arr <- tryCatch(dcm_pixel_array(ds, ...), error = function(e) e)
    rerr <- inherits(arr, "error")
    if (rerr) px$dicomr_error <- conditionMessage(arr)
    perr <- is.null(ref$pixel_shape)
    if (rerr && perr) px$status <- "both_failed"
    else if (perr) px$status <- "pydicom_failed"
    else if (rerr) px$status <- "dicomr_failed"
    else {
      shape <- unlist(ref$pixel_shape)
      p <- readBin(file.path(out, "pixels.bin"), "double", n = prod(shape), endian = "little")
      if (!identical(as.integer(dim(arr) %||% length(arr)), as.integer(shape))) {
        px$status <- "different"
        px$dicomr_error <- sprintf("shape %s vs pydicom %s", paste(dim(arr), collapse = "x"),
                                   paste(shape, collapse = "x"))
      } else {
        d <- abs(as.numeric(arr) - p)
        px$max_abs_diff <- max(d)
        px$fraction_equal <- mean(d == 0)
        px$status <- if (px$max_abs_diff == 0) "equal" else if (px$max_abs_diff <= 2) "close" else "different"
      }
    }
  }
  res$pixels <- px
  res
}

`%||%` <- function(a, b) if (is.null(a)) b else a
