# Reading files and working with data sets ------------------------------------

.ts_names <- c(
  "1.2.840.10008.1.2" = "Implicit VR Little Endian",
  "1.2.840.10008.1.2.1" = "Explicit VR Little Endian",
  "1.2.840.10008.1.2.1.99" = "Deflated Explicit VR Little Endian",
  "1.2.840.10008.1.2.2" = "Explicit VR Big Endian",
  "1.2.840.10008.1.2.4.50" = "JPEG Baseline (Process 1)",
  "1.2.840.10008.1.2.4.51" = "JPEG Extended (Process 2 & 4)",
  "1.2.840.10008.1.2.4.57" = "JPEG Lossless, Non-Hierarchical (Process 14)",
  "1.2.840.10008.1.2.4.70" = "JPEG Lossless, Non-Hierarchical, First-Order Prediction",
  "1.2.840.10008.1.2.4.80" = "JPEG-LS Lossless",
  "1.2.840.10008.1.2.4.81" = "JPEG-LS Lossy (Near-Lossless)",
  "1.2.840.10008.1.2.4.90" = "JPEG 2000 Image Compression (Lossless Only)",
  "1.2.840.10008.1.2.4.91" = "JPEG 2000 Image Compression",
  "1.2.840.10008.1.2.4.201" = "High-Throughput JPEG 2000 (Lossless Only)",
  "1.2.840.10008.1.2.4.203" = "High-Throughput JPEG 2000",
  "1.2.840.10008.1.2.5" = "RLE Lossless"
)

.new_dataset <- function(elements = list()) {
  if (length(elements) == 0) names(elements) <- character(0)
  structure(elements, class = "dicom_dataset")
}

#' Read a DICOM file
#'
#' Parses a DICOM Part 10 file (or a raw data set without the 128-byte
#' preamble when `force = TRUE`) into a `dicom_dataset`.
#'
#' Values are converted to R types following pydicom: text VRs become
#' character vectors (split on backslash), `IS` becomes integer, `DS`, `FL`,
#' `FD`, `UL` become double, `US`, `SS`, `SL` become integer, `AT` becomes a
#' hex tag string, `SQ` becomes a list of data sets, and `OB`/`OW`/`UN` and
#' other binary VRs stay as raw vectors. Encapsulated (compressed) pixel data
#' is kept as its fragments; decode it with [dcm_pixel_array()].
#'
#' @param file Path to the file, or a raw vector holding the file's bytes.
#' @param stop_before_pixels If `TRUE`, stop reading at the Pixel Data
#'   element, which is much faster when you only need the header.
#' @param force Read files that lack the "DICM" prefix.
#' @return A `dicom_dataset`. File level information is available with
#'   [dcm_file_meta()] and [dcm_transfer_syntax()].
#' @examples
#' f <- system.file("extdata", "echo_synthetic.dcm", package = "dicomr")
#' ds <- read_dicom(f)
#' ds$Modality
#' ds[["(0028,0008)"]]          # Number of Frames, by tag
#' hdr <- read_dicom(f, stop_before_pixels = TRUE)
#' dcm_has(hdr, "PixelData")
#' @export
read_dicom <- function(file, stop_before_pixels = FALSE, force = FALSE) {
  if (is.raw(file)) {
    bytes <- file
    fname <- NA_character_
  } else {
    file <- path.expand(file)
    if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)
    bytes <- readBin(file, "raw", n = file.info(file)$size)
    fname <- normalizePath(file)
  }
  res <- .parse_dicom_raw(bytes, stop_before_pixels, force)
  for (w in unique(res$warnings)) warning(w, call. = FALSE)
  ds <- res$dataset
  attr(ds, "file_meta") <- res$meta
  attr(ds, "preamble") <- res$preamble
  attr(ds, "encoding") <- list(
    transfer_syntax = if (nzchar(res$transfer_syntax)) res$transfer_syntax else NA_character_,
    little_endian = res$little_endian,
    implicit_vr = res$implicit_vr,
    encapsulated = res$encapsulated,
    deflated = res$deflated,
    has_meta = res$has_meta
  )
  attr(ds, "filename") <- fname
  ds
}

#' File level information of a data set read from disk
#'
#' @param ds A `dicom_dataset` from [read_dicom()].
#' @return `dcm_file_meta()` returns the File Meta Information (group 0002) as
#'   a `dicom_dataset`. `dcm_transfer_syntax()` returns the transfer syntax
#'   UID with its name as the `name` attribute.
#' @examples
#' ct <- read_dicom(system.file("extdata", "ct_synthetic.dcm", package = "dicomr"))
#' dcm_transfer_syntax(ct)
#' dcm_file_meta(ct)
#' @export
dcm_file_meta <- function(ds) {
  m <- attr(ds, "file_meta")
  if (is.null(m)) .new_dataset() else m
}

#' @rdname dcm_file_meta
#' @export
dcm_transfer_syntax <- function(ds) {
  enc <- attr(ds, "encoding")
  ts <- if (is.null(enc)) NA_character_ else enc$transfer_syntax
  if (is.na(ts)) {
    if (!is.null(enc)) {
      ts <- if (enc$implicit_vr) "1.2.840.10008.1.2" else if (enc$little_endian) "1.2.840.10008.1.2.1" else "1.2.840.10008.1.2.2"
    } else {
      ts <- "1.2.840.10008.1.2.1"
    }
  }
  structure(ts, name = unname(.ts_names[ts]))
}

# --- element access ----------------------------------------------------------

#' Get or set data elements
#'
#' `dcm_get()` returns an element's value, `dcm_element()` the whole element
#' (group, element, VR and value), `dcm_has()` tests for presence, and
#' `dcm_set()` adds or replaces an element. The usual R operators also work:
#' `ds$PatientName`, `ds[["PatientName"]]`, `ds[["0010,0010"]]`,
#' `ds[[c(0x0010, 0x0010)]]`, and the same forms on the left of `<-`.
#' Assigning `NULL` removes the element.
#'
#' @param ds A `dicom_dataset`.
#' @param tag A keyword or tag (see [dcm_tag()]).
#' @param value The new value.
#' @param vr The value representation. Taken from the existing element or
#'   the dictionary when missing; required for private tags.
#' @param default Returned when the element is absent.
#' @return `dcm_get()` the value; `dcm_element()` a `dicom_element` or
#'   `NULL`; `dcm_has()` a logical; `dcm_set()` the modified data set.
#' @examples
#' ct <- read_dicom(system.file("extdata", "ct_synthetic.dcm", package = "dicomr"))
#' dcm_get(ct, "Modality")
#' dcm_get(ct, "(0028,1050)")   # Window Center
#' dcm_get(ct, "BodyPartExamined", default = NA)
#' dcm_has(ct, "PixelData")
#' dcm_element(ct, "PixelSpacing")
#' ct <- dcm_set(ct, "BodyPartExamined", "HEAD")
#' ct$PatientName <- "Anonymous"
#' ct$PatientName
#' @export
dcm_get <- function(ds, tag, default = NULL) {
  el <- dcm_element(ds, tag)
  if (is.null(el)) default else el$value
}

#' @rdname dcm_get
#' @export
dcm_element <- function(ds, tag) {
  hex <- dcm_tag(tag)
  if (!hex %in% names(ds)) return(NULL)
  .subset2(ds, hex)
}

#' @rdname dcm_get
#' @export
dcm_has <- function(ds, tag) dcm_tag(tag) %in% names(ds)

#' @rdname dcm_get
#' @export
dcm_set <- function(ds, tag, value, vr = NULL) {
  hex <- dcm_tag(tag)
  attrs <- attributes(ds)
  els <- unclass(ds)
  if (is.null(value)) {
    els[[hex]] <- NULL
  } else if (inherits(value, "dicom_element")) {
    els[[hex]] <- value
  } else {
    if (is.null(vr)) {
      vr <- if (hex %in% names(els)) els[[hex]]$vr else .dict_vr(hex)
      if (identical(hex, "7FE00010") && !is.null(attrs$encoding) && isTRUE(attrs$encoding$encapsulated))
        vr <- "OB"
    }
    if (vr == "UN" && !is.raw(value) && !(hex %in% names(els)))
      stop("Tag ", dcm_tag_string(hex), " is not in the dictionary; supply vr", call. = FALSE)
    if (vr == "SQ") value <- .as_sequence(value)
    tagn <- .tag_from_hex(hex)
    els[[hex]] <- structure(list(group = as.integer(tagn %/% 65536),
                                 element = as.integer(tagn %% 65536),
                                 vr = vr, value = value), class = "dicom_element")
  }
  els <- els[order(names(els))]
  attrs$names <- names(els)
  attributes(els) <- attrs
  els
}

.as_sequence <- function(value) {
  if (inherits(value, "dicom_dataset")) value <- list(value)
  if (!is.list(value) || !all(vapply(value, inherits, TRUE, "dicom_dataset")))
    stop("A sequence value must be a list of dicom_dataset items", call. = FALSE)
  structure(unclass(value), class = "dicom_sequence")
}

#' Create an empty data set or a sequence item
#'
#' @param ... Named values, using DICOM keywords as names.
#' @return A `dicom_dataset`.
#' @examples
#' item <- dicom_dataset(PhysicalDeltaX = 0.03, PhysicalUnitsXDirection = 3L)
#' ds <- dicom_dataset(Modality = "US", SequenceOfUltrasoundRegions = list(item))
#' ds
#' @export
dicom_dataset <- function(...) {
  args <- list(...)
  ds <- .new_dataset()
  for (nm in names(args)) ds <- dcm_set(ds, nm, args[[nm]])
  ds
}

#' @export
`$.dicom_dataset` <- function(x, name) dcm_get(x, name)

#' @export
`[[.dicom_dataset` <- function(x, i, ...) dcm_get(x, i)

#' @export
`$<-.dicom_dataset` <- function(x, name, value) dcm_set(x, name, value)

#' @export
`[[<-.dicom_dataset` <- function(x, i, value) dcm_set(x, i, value)

#' @export
`[.dicom_dataset` <- function(x, i) {
  hex <- vapply(as.list(i), dcm_tag, "")
  hex <- hex[hex %in% names(x)]
  attrs <- attributes(x)
  els <- unclass(x)[hex]
  attrs$names <- hex
  attributes(els) <- attrs
  els
}

#' Keywords of the elements in a data set
#' @param ds A `dicom_dataset`.
#' @return A character vector (empty strings for private or unknown tags).
#' @examples
#' ct <- read_dicom(system.file("extdata", "ct_synthetic.dcm", package = "dicomr"))
#' head(dcm_keywords(ct))
#' @export
dcm_keywords <- function(ds) {
  if (length(ds) == 0) return(character(0))
  .dict_info(names(ds))$keyword
}

#' @export
as.list.dicom_dataset <- function(x, ...) {
  kw <- dcm_keywords(x)
  out <- lapply(unclass(x), function(el) {
    v <- el$value
    if (inherits(v, "dicom_sequence")) lapply(v, as.list.dicom_dataset) else v
  })
  names(out) <- ifelse(nzchar(kw), kw, dcm_tag_string(names(x)))
  out
}

# --- printing ----------------------------------------------------------------

.format_value <- function(v, vr, width = 60) {
  if (is.null(v) || length(v) == 0) return("")
  if (inherits(v, "dicom_encapsulated"))
    return(sprintf("Encapsulated data: %d fragment(s), %s bytes", length(v$fragments),
                   format(sum(lengths(v$fragments)), big.mark = ",")))
  if (is.raw(v)) {
    if (length(v) <= 16) return(paste0("0x", paste(format(v), collapse = "")))
    return(sprintf("Array of %s bytes", format(length(v), big.mark = ",")))
  }
  if (is.character(v)) {
    s <- paste0("'", v, "'")
  } else {
    s <- format(v, digits = 10, trim = TRUE)
  }
  if (length(s) > 1) {
    if (length(s) > 10) s <- c(s[1:10], "...")
    s <- paste0("[", paste(s, collapse = ", "), "]")
  }
  s <- gsub("[\r\n\t]+", " ", s)
  if (nchar(s, type = "width") > width) s <- paste0(substr(s, 1, width - 3), "...")
  s
}

.print_lines <- function(ds, indent = 0) {
  if (length(ds) == 0) return(character(0))
  info <- .dict_info(names(ds))
  pad <- strrep("   ", indent)
  out <- character(0)
  for (k in seq_along(ds)) {
    el <- .subset2(ds, k)
    hex <- names(ds)[k]
    nm <- info$name[k]
    if (nm == "Unknown") nm <- "Private tag data"
    tagstr <- dcm_tag_string(hex)
    if (el$vr == "SQ") {
      out <- c(out, sprintf("%s%s  %s  %d item(s) ----", pad, tagstr, nm, length(el$value)))
      for (item in el$value) {
        out <- c(out, .print_lines(item, indent + 1), paste0(pad, "   ---------"))
      }
    } else {
      out <- c(out, sprintf("%s%s %-*s %s: %s", pad, tagstr, max(1, 40 - 3 * indent), nm,
                            el$vr, .format_value(el$value, el$vr)))
    }
  }
  out
}

#' @export
print.dicom_dataset <- function(x, max_lines = 200, ...) {
  meta <- attr(x, "file_meta")
  lines <- character(0)
  if (!is.null(meta) && length(meta) > 0) {
    lines <- c("Dataset.file_meta -------------------------------", .print_lines(meta),
               "-------------------------------------------------")
  }
  lines <- c(lines, .print_lines(x))
  if (length(lines) > max_lines) {
    cat(lines[seq_len(max_lines)], sep = "\n")
    cat(sprintf("... %d more lines (print(ds, max_lines = Inf) to see all)\n", length(lines) - max_lines))
  } else {
    cat(lines, sep = "\n")
  }
  invisible(x)
}

#' @export
print.dicom_element <- function(x, ...) {
  hex <- dcm_tag(x)
  info <- .dict_info(hex)
  if (x$vr == "SQ") {
    cat(sprintf("%s %s  SQ: %d item(s)\n", dcm_tag_string(hex), info$name, length(x$value)))
  } else {
    cat(sprintf("%s %s  %s: %s\n", dcm_tag_string(hex), info$name, x$vr, .format_value(x$value, x$vr)))
  }
  invisible(x)
}

#' @export
print.dicom_sequence <- function(x, ...) {
  cat(sprintf("<dicom_sequence with %d item(s)>\n", length(x)))
  for (i in seq_along(x)) {
    cat(sprintf("--- item %d ---\n", i))
    print(x[[i]])
  }
  invisible(x)
}

#' @export
print.dicom_encapsulated <- function(x, ...) {
  cat(.format_value(x), "\n")
  invisible(x)
}

#' Flatten a data set into a data frame
#'
#' One row per element, with sequence items flattened and marked by `depth`
#' and `path`. Values are formatted as text in the `value` column.
#'
#' @param x A `dicom_dataset`.
#' @param row.names,optional Unused, for compatibility with the generic.
#' @param ... Unused.
#' @return A data frame with columns `tag`, `keyword`, `name`, `vr`, `vm`,
#'   `value`, `depth` and `path`.
#' @examples
#' echo <- read_dicom(system.file("extdata", "echo_synthetic.dcm", package = "dicomr"))
#' df <- as.data.frame(echo)
#' head(df[, c("tag", "keyword", "vr", "value")])
#' @export
as.data.frame.dicom_dataset <- function(x, row.names = NULL, optional = FALSE, ...) {
  rows <- list()
  walk <- function(ds, depth, path) {
    if (length(ds) == 0) return(invisible())
    info <- .dict_info(names(ds))
    for (k in seq_along(ds)) {
      el <- .subset2(ds, k)
      v <- el$value
      vm <- if (el$vr == "SQ") length(v) else if (is.raw(v) || inherits(v, "dicom_encapsulated")) 1L else length(v)
      rows[[length(rows) + 1]] <<- data.frame(
        tag = dcm_tag_string(names(ds)[k]), keyword = info$keyword[k], name = info$name[k],
        vr = el$vr, vm = vm,
        value = if (el$vr == "SQ") sprintf("%d item(s)", length(v)) else .format_value(v, el$vr, width = 200),
        depth = depth, path = path, stringsAsFactors = FALSE)
      if (el$vr == "SQ") {
        for (i in seq_along(v)) walk(v[[i]], depth + 1L, paste0(path, if (nzchar(path)) "/" else "",
                                                               info$keyword[k], "[", i, "]"))
      }
    }
  }
  walk(x, 0L, "")
  if (length(rows) == 0) return(data.frame(tag = character(0), keyword = character(0)))
  do.call(rbind, rows)
}
