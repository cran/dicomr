make_us_dataset <- function(frames = 3, rows = 20, cols = 30) {
  arr <- array(0L, c(frames, rows, cols, 3))
  for (f in seq_len(frames)) arr[f, , , 1] <- f
  arr[, 5, 7, 2] <- 200L
  arr[, , , 3] <- outer(rep(1, frames), outer(seq_len(rows), seq_len(cols), "+"))
  region <- dicom_dataset(
    RegionSpatialFormat = 1L, RegionDataType = 1L, RegionFlags = 0,
    RegionLocationMinX0 = 2, RegionLocationMinY0 = 3,
    RegionLocationMaxX1 = 27, RegionLocationMaxY1 = 18,
    PhysicalUnitsXDirection = 3L, PhysicalUnitsYDirection = 3L,
    PhysicalDeltaX = 0.05, PhysicalDeltaY = 0.05)
  doppler <- dicom_dataset(
    RegionSpatialFormat = 3L, RegionDataType = 3L,
    RegionLocationMinX0 = 0, RegionLocationMinY0 = 19,
    RegionLocationMaxX1 = 29, RegionLocationMaxY1 = 19,
    ReferencePixelX0 = 0L, ReferencePixelY0 = 0L,
    PhysicalUnitsXDirection = 4L, PhysicalUnitsYDirection = 7L,
    ReferencePixelPhysicalValueX = 0, ReferencePixelPhysicalValueY = 0,
    PhysicalDeltaX = 0.01, PhysicalDeltaY = -2.5)
  ds <- dicom_dataset(
    SOPClassUID = "1.2.840.10008.5.1.4.1.1.3.1",
    SOPInstanceUID = "2.25.1234567890",
    Modality = "US",
    PatientName = "Test^Echo",
    FrameTime = 33.3,
    ImageType = c("ORIGINAL", "PRIMARY"),
    SequenceOfUltrasoundRegions = list(region, doppler))
  dcm_set_pixel_data(ds, arr)
}
