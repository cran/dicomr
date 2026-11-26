# Tags, keywords and the data dictionary --------------------------------------

.hex8 <- function(tag) {
  tag <- as.numeric(tag)
  sprintf("%04X%04X", as.integer(tag %/% 65536), as.integer(tag %% 65536))
}

.tag_from_hex <- function(hex) as.numeric(paste0("0x", hex))

#' Convert a tag specification to a DICOM tag
#'
#' Accepts a keyword (`"PatientName"`), a tag string (`"(0010,0010)"`,
#' `"0010,0010"`, `"00100010"`, `"0x00100010"`), a single number
#' (`0x00100010`) or a group/element pair (`c(0x0010, 0x0010)`).
#'
#' @param x A tag specification.
#' @return `dcm_tag()` returns the tag as an 8-character hex string such as
#'   `"00100010"`, which is how elements are named inside a data set.
#'   `dcm_tag_number()` returns the tag as a number.
#' @examples
#' dcm_tag("PatientName")
#' dcm_tag(c(0x0028, 0x0010))
#' dcm_tag("(7FE0,0010)")
#' dcm_tag_number("Rows")
#' @export
dcm_tag <- function(x) {
  if (inherits(x, "dicom_element")) return(.hex8(x$group * 65536 + x$element))
  if (is.numeric(x)) {
    if (length(x) == 2) return(.hex8(x[1] * 65536 + x[2]))
    if (length(x) == 1) return(.hex8(x))
    stop("A numeric tag must be one number or a (group, element) pair", call. = FALSE)
  }
  if (!is.character(x) || length(x) != 1 || is.na(x))
    stop("A tag must be a keyword, a tag string or a number", call. = FALSE)
  s <- gsub("[[:space:]()]", "", x)
  s <- sub("^0[xX]", "", s)
  if (grepl("^[0-9A-Fa-f]{4},[0-9A-Fa-f]{4}$", s)) return(toupper(sub(",", "", s)))
  tag <- .dict_keyword_to_tag(x)
  if (!is.na(tag)) return(.hex8(tag))
  if (grepl("^[0-9A-Fa-f]{8}$", s)) return(toupper(s))
  stop(sprintf("'%s' is not a DICOM keyword or tag", x), call. = FALSE)
}

#' @rdname dcm_tag
#' @export
dcm_tag_number <- function(x) .tag_from_hex(dcm_tag(x))

#' Format a tag as "(gggg,eeee)"
#' @param x A tag specification (see [dcm_tag()]).
#' @return A character string.
#' @examples
#' dcm_tag_string("PatientName")
#' dcm_tag_string(0x7FE00010)
#' @export
dcm_tag_string <- function(x) {
  h <- if (is.character(x) && all(grepl("^[0-9A-F]{8}$", x))) x else vapply(x, dcm_tag, "")
  sprintf("(%s,%s)", substr(h, 1, 4), substr(h, 5, 8))
}

#' Look up tags or keywords in the DICOM data dictionary
#'
#' @param x Character keywords or tag specifications. If missing, the whole
#'   dictionary is returned.
#' @return A data frame with columns `tag`, `vr`, `vm`, `name`, `keyword` and
#'   `retired`.
#' @examples
#' dicom_dictionary(c("PatientName", "PhysicalDeltaX", "(7FE0,0010)"))
#' nrow(dicom_dictionary())
#' @export
dicom_dictionary <- function(x) {
  if (missing(x)) {
    d <- .dict_table()
    d$tag <- dcm_tag_string(.hex8(d$tag))
    return(d)
  }
  hex <- vapply(as.list(x), dcm_tag, "")
  d <- .dict_lookup_tags(.tag_from_hex(hex))
  cbind(data.frame(tag = dcm_tag_string(hex), stringsAsFactors = FALSE), d)
}

.dict_info <- function(hex) .dict_lookup_tags(.tag_from_hex(hex))

.dict_vr <- function(hex) {
  vr <- .dict_info(hex)$vr
  if (grepl(" or ", vr)) vr <- substr(vr, 1, 2)
  vr
}
