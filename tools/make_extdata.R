# Builds the synthetic example files in inst/extdata. They contain no
# patient data. Run from the package root with dicomr installed and
# Python 3 with Pillow and numpy available (for the JPEG frames):
#   Rscript tools/make_extdata.R
library(dicomr)

# 1. echo_synthetic.dcm: 4-frame colour ultrasound cine, JPEG Baseline
#    (YBR_FULL_422), with a 2D region and a spectral Doppler region.
frames <- 4L; rows <- 48L; cols <- 64L
truth <- array(0L, c(frames, rows, cols, 3L))
for (f in seq_len(frames)) {
  truth[f, , , 1] <- outer(seq_len(rows), seq_len(cols), function(r, c) 40L + 3L * r + 10L * f)
  truth[f, , , 2] <- outer(seq_len(rows), seq_len(cols), function(r, c) 30L + 2L * c)
  truth[f, , , 3] <- 120L + 20L * f
}
storage.mode(truth) <- "integer"
tmp <- tempfile("frames"); dir.create(tmp)
for (f in seq_len(frames))
  writeBin(as.raw(aperm(truth[f, , , ], c(3, 2, 1))), file.path(tmp, sprintf("f%d.rgb", f)))
py <- "
import sys, glob, os, numpy as np
from PIL import Image
d, rows, cols = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
for p in sorted(glob.glob(os.path.join(d, '*.rgb'))):
    a = np.fromfile(p, dtype=np.uint8).reshape(rows, cols, 3)
    Image.fromarray(a, 'RGB').save(p[:-4] + '.jpg', quality=95, subsampling=2)
"
stopifnot(system2("python3", c("-c", shQuote(py), tmp, rows, cols)) == 0)
frags <- lapply(sprintf(file.path(tmp, "f%d.jpg"), seq_len(frames)),
                function(p) readBin(p, "raw", file.size(p)))

region <- dicom_dataset(
  RegionSpatialFormat = 1L, RegionDataType = 1L, RegionFlags = 0,
  RegionLocationMinX0 = 0, RegionLocationMinY0 = 0,
  RegionLocationMaxX1 = 63, RegionLocationMaxY1 = 31,
  ReferencePixelX0 = 32L, ReferencePixelY0 = 0L,
  PhysicalUnitsXDirection = 3L, PhysicalUnitsYDirection = 3L,
  PhysicalDeltaX = 0.05, PhysicalDeltaY = 0.05)
doppler <- dicom_dataset(
  RegionSpatialFormat = 3L, RegionDataType = 3L, RegionFlags = 0,
  RegionLocationMinX0 = 0, RegionLocationMinY0 = 32,
  RegionLocationMaxX1 = 63, RegionLocationMaxY1 = 47,
  ReferencePixelX0 = 0L, ReferencePixelY0 = 40L,
  PhysicalUnitsXDirection = 4L, PhysicalUnitsYDirection = 7L,
  ReferencePixelPhysicalValueX = 0, ReferencePixelPhysicalValueY = 0,
  PhysicalDeltaX = 0.01, PhysicalDeltaY = -2.5)
ds <- dicom_dataset(
  SOPClassUID = "1.2.840.10008.5.1.4.1.1.3.1",
  SOPInstanceUID = "2.25.311960424838391650129385862651349530918",
  StudyInstanceUID = "2.25.98213846611120342190374560221883160351",
  SeriesInstanceUID = "2.25.205476372381904772180633130861478201862",
  Modality = "US", Manufacturer = "dicomr synthetic example",
  PatientName = "Synthetic^Echo", PatientID = "DICOMR-0001",
  ImageType = c("ORIGINAL", "PRIMARY"),
  FrameTime = 33.3, FrameIncrementPointer = "FrameTime",
  NumberOfFrames = frames, Rows = rows, Columns = cols,
  SamplesPerPixel = 3L, PhotometricInterpretation = "YBR_FULL_422",
  PlanarConfiguration = 0L, BitsAllocated = 8L, BitsStored = 8L,
  HighBit = 7L, PixelRepresentation = 0L, LossyImageCompression = "01",
  SequenceOfUltrasoundRegions = list(region, doppler))
enc <- structure(list(offsets = numeric(0), fragments = frags), class = "dicom_encapsulated")
ds <- dcm_set(ds, "PixelData", enc, vr = "OB")
attr(ds, "encoding") <- list(transfer_syntax = "1.2.840.10008.1.2.4.50", little_endian = TRUE,
                             implicit_vr = FALSE, encapsulated = TRUE, deflated = FALSE,
                             has_meta = TRUE)
write_dicom(ds, "inst/extdata/echo_synthetic.dcm")
unlink(tmp, recursive = TRUE)

# 2. ct_synthetic.dcm: 32 x 32 signed 16-bit greyscale slice with a
#    rescale to Hounsfield units, implicit VR little endian.
m <- outer(1:32, 1:32, function(r, c) as.integer(1024 - 1000 * ((r - 16.5)^2 + (c - 16.5)^2 > 150)))
storage.mode(m) <- "integer"
ct <- dicom_dataset(
  SOPClassUID = "1.2.840.10008.5.1.4.1.1.2",
  SOPInstanceUID = "2.25.137402914860427738460178035870322497171",
  Modality = "CT", Manufacturer = "dicomr synthetic example",
  PatientName = "Synthetic^Phantom", PatientID = "DICOMR-0002",
  PixelSpacing = c(0.5, 0.5), SliceThickness = 1,
  RescaleSlope = 1, RescaleIntercept = -1024,
  WindowCenter = 40, WindowWidth = 400)
ct <- dcm_set_pixel_data(ct, m)
write_dicom(ct, "inst/extdata/ct_synthetic.dcm", transfer_syntax = "1.2.840.10008.1.2")
