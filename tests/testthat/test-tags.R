test_that("tag specifications all resolve to the same tag", {
  expect_equal(dcm_tag("PatientName"), "00100010")
  expect_equal(dcm_tag("(0010,0010)"), "00100010")
  expect_equal(dcm_tag("0010,0010"), "00100010")
  expect_equal(dcm_tag("0x00100010"), "00100010")
  expect_equal(dcm_tag(0x00100010), "00100010")
  expect_equal(dcm_tag(c(0x0010, 0x0010)), "00100010")
  expect_equal(dcm_tag("7fe00010"), "7FE00010")
  expect_equal(dcm_tag_string("PixelData"), "(7FE0,0010)")
  expect_error(dcm_tag("NotAKeyword"), "not a DICOM keyword")
})

test_that("the dictionary knows VRs, repeaters and private tags", {
  d <- dicom_dictionary(c("PhysicalDeltaX", "SequenceOfUltrasoundRegions", "(6000,3000)", "(0009,0010)", "(0009,1001)"))
  expect_equal(d$vr, c("FD", "SQ", "OB or OW", "LO", "UN"))
  expect_equal(d$keyword[1:3], c("PhysicalDeltaX", "SequenceOfUltrasoundRegions", "OverlayData"))
  expect_gt(nrow(dicom_dictionary()), 5000)
})
