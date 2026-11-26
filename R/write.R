# Writing ---------------------------------------------------------------------

.impl_uid <- "2.25.147690282439102293357200469405707667985"
.impl_version <- "DICOMR_010"

.long_vrs <- c("OB", "OD", "OF", "OL", "OV", "OW", "SQ", "UC", "UN", "UR", "UT", "SV", "UV")
.text_vrs <- c("AE", "AS", "CS", "DA", "DT", "LO", "LT", "PN", "SH", "ST", "TM", "UC", "UI", "UR", "UT")
.word_size <- c(OW = 2L, OF = 4L, OL = 4L, OD = 8L, OV = 8L)

.u16 <- function(x, endian) writeBin(as.integer(x), raw(), size = 2, endian = endian)
.u32 <- function(x, endian) {
  x <- as.numeric(x)
  x[x >= 2^31] <- x[x >= 2^31] - 2^32
  writeBin(as.integer(x), raw(), size = 4, endian = endian)
}
.u64 <- function(x, endian, signed = FALSE) {
  x <- as.numeric(x)
  x[x < 0] <- x[x < 0] + 2^64
  out <- raw(0)
  for (v in x) {
    hi <- floor(v / 2^32)
    lo <- v - hi * 2^32
    b <- if (endian == "little") c(.u32(lo, endian), .u32(hi, endian)) else c(.u32(hi, endian), .u32(lo, endian))
    out <- c(out, b)
  }
  out
}

.format_ds <- function(x) {
  vapply(x, function(v) {
    if (is.na(v)) return("")
    for (d in 15:1) {
      s <- format(v, digits = d, scientific = FALSE, trim = TRUE)
      if (nchar(s) > 16) s <- format(v, digits = d, scientific = TRUE, trim = TRUE)
      if (nchar(s) <= 16) return(s)
    }
    stop("Can't fit DS value ", v, " in 16 characters", call. = FALSE)
  }, "")
}

.text_bytes <- function(s, charset) {
  s <- enc2utf8(s)
  if (!is.null(charset) && any(charset == "ISO_IR 192")) return(charToRaw(s))
  # Default repertoire and ISO_IR 100: write Latin-1, as pydicom does
  if (any(charToRaw(s) > as.raw(127))) s <- iconv(s, "UTF-8", "latin1", sub = "?")
  charToRaw(s)
}

.encode_value <- function(v, vr, endian, charset, swap_from = NULL, word = NULL) {
  if (is.null(v) || length(v) == 0) return(raw(0))
  if (vr %in% .text_vrs || ((vr %in% c("IS", "DS")) && is.character(v))) {
    # Encode each value separately: paste() would re-encode non-ASCII text
    # to the native locale.
    v <- as.character(v)
    v[is.na(v)] <- ""
    pieces <- lapply(v, .text_bytes, charset = charset)
    b <- pieces[[1]]
    for (p in pieces[-1]) b <- c(b, as.raw(0x5C), p)
    if (length(b) %% 2) b <- c(b, if (vr == "UI") as.raw(0) else charToRaw(" "))
    return(b)
  }
  switch(vr,
    IS = { s <- paste(ifelse(is.na(v), "", format(as.numeric(v), scientific = FALSE, trim = TRUE)), collapse = "\\")
           b <- charToRaw(s); if (length(b) %% 2) b <- c(b, charToRaw(" ")); b },
    DS = { b <- charToRaw(paste(.format_ds(as.numeric(v)), collapse = "\\"))
           if (length(b) %% 2) b <- c(b, charToRaw(" ")); b },
    US = , SS = .u16(v, endian),
    UL = .u32(v, endian),
    SL = writeBin(as.integer(v), raw(), size = 4, endian = endian),
    SV = .u64(v, endian, TRUE),
    UV = .u64(v, endian),
    FL = writeBin(as.numeric(v), raw(), size = 4, endian = endian),
    FD = writeBin(as.numeric(v), raw(), size = 8, endian = endian),
    AT = {
      tags <- .tag_from_hex(vapply(as.list(v), dcm_tag, ""))
      unlist(lapply(tags, function(t) c(.u16(t %/% 65536, endian), .u16(t %% 65536, endian))))
    },
    {
      if (!is.raw(v)) stop("Values with VR ", vr, " must be raw vectors", call. = FALSE)
      b <- v
      ws <- if (!is.null(word)) word else .word_size[vr]
      if (!is.null(swap_from) && swap_from != endian && !is.na(ws) && length(b) %% ws == 0) {
        m <- matrix(b, nrow = ws)
        b <- as.vector(m[ws:1, , drop = FALSE])
      }
      if (length(b) %% 2) b <- c(b, as.raw(0))
      b
    })
}

.element_header <- function(group, element, vr, len, endian, implicit) {
  tag <- c(.u16(group, endian), .u16(element, endian))
  if (implicit) return(c(tag, .u32(len, endian)))
  if (vr %in% .long_vrs) return(c(tag, charToRaw(vr), as.raw(c(0, 0)), .u32(len, endian)))
  if (len > 65535) stop(sprintf("Value of %s (VR %s) is too long for its VR", dcm_tag_string(.hex8(group * 65536 + element)), vr), call. = FALSE)
  c(tag, charToRaw(vr), .u16(len, endian))
}

.item_tag <- function(element, len, endian) c(.u16(0xFFFE, endian), .u16(element, endian), .u32(len, endian))

.encode_dataset <- function(ds, endian, implicit, charset = NULL, swap_from = NULL) {
  if (length(ds) == 0) return(raw(0))
  ba <- if ("00280100" %in% names(ds)) .subset2(ds, "00280100")$value else NULL
  # Swap pixel data by its sample size (OW is nominally 16-bit words)
  pix_word <- if (length(ba) && ba[1] > 16) as.integer(ba[1] %/% 8) else NULL
  hex <- sort(names(ds))
  cs_el <- if ("00080005" %in% hex) .subset2(ds, "00080005")$value else charset
  parts <- vector("list", length(hex))
  for (k in seq_along(hex)) {
    el <- .subset2(ds, hex[k])
    vr <- el$vr
    v <- el$value
    if (vr == "SQ") {
      items <- lapply(v, function(item) {
        c(.item_tag(0xE000, 0xFFFFFFFF, endian),
          .encode_dataset(item, endian, implicit, cs_el, swap_from),
          .item_tag(0xE00D, 0, endian))
      })
      body <- c(unlist(items), .item_tag(0xE0DD, 0, endian))
      if (is.null(body)) body <- .item_tag(0xE0DD, 0, endian)
      parts[[k]] <- c(.element_header(el$group, el$element, vr, 0xFFFFFFFF, endian, implicit), body)
    } else if (inherits(v, "dicom_encapsulated")) {
      if (implicit) stop("Encapsulated pixel data needs an explicit VR transfer syntax", call. = FALSE)
      bot <- if (length(v$offsets)) .u32(v$offsets, endian) else raw(0)
      frs <- lapply(v$fragments, function(f) {
        if (length(f) %% 2) f <- c(f, as.raw(0))
        c(.item_tag(0xE000, length(f), endian), f)
      })
      body <- c(.item_tag(0xE000, length(bot), endian), bot, unlist(frs), .item_tag(0xE0DD, 0, endian))
      parts[[k]] <- c(.element_header(el$group, el$element, if (vr %in% c("OB", "OW")) "OB" else vr,
                                      0xFFFFFFFF, endian, implicit), body)
    } else {
      b <- .encode_value(v, vr, endian, cs_el, swap_from,
                         if (hex[k] == "7FE00010" && vr == "OW") pix_word else NULL)
      parts[[k]] <- c(.element_header(el$group, el$element, vr, length(b), endian, implicit), b)
    }
  }
  unlist(parts)
}

#' Write a data set to a DICOM file
#'
#' Writes a DICOM Part 10 file: 128-byte preamble, "DICM", the File Meta
#' Information group (rebuilt with a correct group length), then the data
#' set. Sequences are written with undefined lengths. Compressed pixel data is
#' written back as it was read; changing between compressed and uncompressed
#' transfer syntaxes is not supported.
#'
#' @param ds A `dicom_dataset`.
#' @param file Output path.
#' @param transfer_syntax Transfer syntax UID. Defaults to the one the data
#'   set was read with, or Explicit VR Little Endian. Uncompressed data sets
#'   can be written as Implicit VR Little Endian (`"1.2.840.10008.1.2"`),
#'   Explicit VR Little Endian (`"1.2.840.10008.1.2.1"`), Explicit VR Big
#'   Endian (`"1.2.840.10008.1.2.2"`) or Deflated (`"1.2.840.10008.1.2.1.99"`).
#' @return The path, invisibly.
#' @examples
#' ct <- read_dicom(system.file("extdata", "ct_synthetic.dcm", package = "dicomr"))
#' out <- tempfile(fileext = ".dcm")
#' write_dicom(ct, out, transfer_syntax = "1.2.840.10008.1.2.1")
#' dcm_transfer_syntax(read_dicom(out))
#' unlink(out)
#' @export
write_dicom <- function(ds, file, transfer_syntax = NULL) {
  enc <- attr(ds, "encoding")
  src_ts <- as.character(dcm_transfer_syntax(ds))
  if (is.null(transfer_syntax)) transfer_syntax <- src_ts
  uncompressed <- c("1.2.840.10008.1.2", "1.2.840.10008.1.2.1", "1.2.840.10008.1.2.2", "1.2.840.10008.1.2.1.99")
  pix <- dcm_get(ds, "PixelData")
  is_encap <- inherits(pix, "dicom_encapsulated")
  if (is_encap && transfer_syntax %in% uncompressed)
    stop("Pixel data is compressed; it can only be written with its own transfer syntax (", src_ts, ")", call. = FALSE)
  if (!is_encap && !is.null(pix) && !(transfer_syntax %in% uncompressed))
    stop("Re-encoding uncompressed pixel data as ", transfer_syntax, " is not supported", call. = FALSE)

  implicit <- transfer_syntax == "1.2.840.10008.1.2"
  endian <- if (transfer_syntax == "1.2.840.10008.1.2.2") "big" else "little"
  src_endian <- if (!is.null(enc) && !isTRUE(enc$little_endian)) "big" else "little"

  # File meta
  meta <- dcm_file_meta(ds)
  meta <- dcm_set(meta, "00020001", as.raw(c(0, 1)), vr = "OB")
  sop_class <- dcm_get(ds, "SOPClassUID")
  sop_inst <- dcm_get(ds, "SOPInstanceUID")
  if (!is.null(sop_class)) meta <- dcm_set(meta, "00020002", sop_class, vr = "UI")
  if (!is.null(sop_inst)) meta <- dcm_set(meta, "00020003", sop_inst, vr = "UI")
  meta <- dcm_set(meta, "00020010", transfer_syntax, vr = "UI")
  if (!dcm_has(meta, "00020012")) meta <- dcm_set(meta, "00020012", .impl_uid, vr = "UI")
  if (!dcm_has(meta, "00020013")) meta <- dcm_set(meta, "00020013", .impl_version, vr = "SH")
  meta <- dcm_set(meta, "00020000", NULL)
  meta_bytes <- .encode_dataset(meta, "little", FALSE)
  meta_bytes <- c(.element_header(0x0002, 0x0000, "UL", 4, "little", FALSE), .u32(length(meta_bytes), "little"), meta_bytes)

  body_ds <- ds
  body_ds <- body_ds[names(body_ds)[substr(names(body_ds), 1, 4) != "0002"]]
  body <- .encode_dataset(body_ds, endian, implicit, NULL, src_endian)
  if (transfer_syntax == "1.2.840.10008.1.2.1.99") body <- .deflate_raw(body)

  pre <- attr(ds, "preamble")
  if (!is.raw(pre) || length(pre) != 128) pre <- raw(128)
  con <- file(file, "wb")
  on.exit(close(con))
  writeBin(c(pre, charToRaw("DICM"), meta_bytes, body), con)
  invisible(file)
}
