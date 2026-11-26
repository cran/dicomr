# Compare dicomr with pydicom on every DICOM file in a folder.
# Usage: Rscript compare_all.R <folder> [output.csv]
library(dicomr)
args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[1] else {
  system2("python3", c("-c", shQuote("import pydicom.data, os; print(os.path.dirname(pydicom.data.get_testdata_file('CT_small.dcm')))")), stdout = TRUE)
}
files <- sort(list.files(dir, pattern = "\\.dcm$", full.names = TRUE))
rows <- lapply(files, function(f) {
  r <- tryCatch(dcm_compare_pydicom(f), error = function(e) e)
  if (inherits(r, "error")) {
    return(data.frame(file = basename(f), ts = NA, elements = NA, tag_issues = NA,
                      pixels = "error", max_abs_diff = NA, fraction_equal = NA,
                      note = substr(conditionMessage(r), 1, 120)))
  }
  ts <- tryCatch(attr(dcm_transfer_syntax(suppressWarnings(read_dicom(f, TRUE, force = TRUE))), "name"), error = function(e) NA)
  note <- c(if (nrow(r$issues)) paste(unique(paste(r$issues$tag, r$issues$issue))[1:min(3, nrow(r$issues))], collapse = "; "),
            r$pixels$dicomr_error, r$pixels$pydicom_error)
  data.frame(file = basename(f), ts = if (is.null(ts)) NA else ts, elements = r$elements,
             tag_issues = nrow(r$issues), pixels = r$pixels$status,
             max_abs_diff = r$pixels$max_abs_diff, fraction_equal = r$pixels$fraction_equal,
             note = substr(paste(note, collapse = " | "), 1, 160))
})
res <- do.call(rbind, rows)
print(res[, c("file", "elements", "tag_issues", "pixels", "max_abs_diff", "fraction_equal")], row.names = FALSE)
if (length(args) > 1) write.csv(res, args[2], row.names = FALSE)
