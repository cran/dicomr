test_that("elements can be read, set and removed", {
  ds <- dicom_dataset(PatientName = "Doe^Jane", Rows = 10L)
  expect_equal(ds$PatientName, "Doe^Jane")
  expect_equal(ds[["0028,0010"]], 10L)
  expect_equal(dcm_element(ds, "Rows")$vr, "US")
  ds$PatientID <- "123"
  ds[["StudyDate"]] <- "20261101"
  expect_true(dcm_has(ds, "PatientID"))
  expect_equal(names(ds), sort(names(ds)))
  ds$PatientID <- NULL
  expect_false(dcm_has(ds, "PatientID"))
  expect_null(ds$PatientID)
  expect_error(dcm_set(ds, "(0009,1001)", "x"), "supply vr")
  ds <- dcm_set(ds, "(0009,1001)", "x", vr = "LO")
  expect_equal(ds[["(0009,1001)"]], "x")
})

test_that("printing and data frame conversion work, including sequences", {
  ds <- make_us_dataset()
  out <- capture.output(print(ds))
  expect_true(any(grepl("Sequence of Ultrasound Regions  2 item\\(s\\)", out)))
  expect_true(any(grepl("Physical Delta X", out)))
  df <- as.data.frame(ds)
  expect_true(all(c("tag", "keyword", "vr", "value", "depth", "path") %in% names(df)))
  expect_equal(sum(df$depth == 1 & df$keyword == "PhysicalDeltaX"), 2)
})
