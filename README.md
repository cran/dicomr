# dicomr

Read, write and decode DICOM files natively in R, modelled on Python's
[pydicom](https://github.com/pydicom/pydicom). The file parser and image
decoders are C++ (via Rcpp and libjpeg-turbo); everything else is plain R.
Built echo-first: multi-frame JPEG cine loops in YBR colour, ultrasound
calibration regions and frame timing all work out of the box.

## Install

dicomr needs the libjpeg and zlib development files when installed from
source. Windows (Rtools) and the CRAN macOS builds already include them.

```sh
# Debian / Ubuntu
sudo apt-get install libjpeg-dev zlib1g-dev
# Fedora / RHEL
sudo dnf install libjpeg-turbo-devel zlib-devel
# macOS (Homebrew), only when building from source
brew install jpeg-turbo
```

Then from R, pointing at the package source:

```r
install.packages("Rcpp")
install.packages("path/to/dicomr_0.1.0.tar.gz", repos = NULL, type = "source")
```

If libjpeg lives somewhere unusual, set `JPEG_CFLAGS` and `JPEG_LIBS` (or
`INCLUDE_DIR` and `LIB_DIR`) before installing.

Two small synthetic files ship with the package for trying things out:

```r
echo <- read_dicom(system.file("extdata", "echo_synthetic.dcm", package = "dicomr"))
ct   <- read_dicom(system.file("extdata", "ct_synthetic.dcm", package = "dicomr"))
```

## Quick start (echo)

```r
library(dicomr)

ds <- read_dicom("echo.dcm")
ds                                   # pydicom-style listing of every element
ds$Manufacturer                      # by keyword
ds[["(0028,0008)"]]                  # by tag
dcm_transfer_syntax(ds)              # e.g. JPEG Baseline (Process 1)

clip <- dcm_pixel_array(ds)          # frames x rows x cols x 3, RGB
dim(clip)                            # 30 240 320 3
first <- dcm_pixel_array(ds, frames = 1)   # decode only what you need
dcm_show(ds, frame = 1)

dcm_us_regions(ds)                   # region table with cm/pixel calibration
dcm_us_distance(ds, c(120, 150), c(180, 150))   # caliper distance in cm
dcm_frame_times(ds)                  # ms per frame, frame rate attribute

ds$PatientName <- "Anon^Anon"        # edit, then write back out
write_dicom(ds, "echo_anon.dcm")
```

Header only (fast, skips pixel data):

```r
hdr <- read_dicom("echo.dcm", stop_before_pixels = TRUE)
as.data.frame(hdr)                   # one row per element, sequences flattened
```

## What it covers

| Area | Support |
|---|---|
| File format | Preamble or none, File Meta, implicit/explicit VR, little/big endian, deflated, nested sequences (defined and undefined length), UN sequences, encapsulated fragments, basic and extended offset tables |
| Dictionary | Full DICOM PS3.6 dictionary (5,091 tags plus repeating groups), generated from pydicom's |
| Values | Text split on `\`, IS to integer, DS to double, binary as raw; ISO_IR 100 (Latin-1) and ISO_IR 192 (UTF-8) text |
| Pixel data | Uncompressed 1, 8, 16, 32, 64-bit, signed/unsigned, float; planar configuration; RLE Lossless; JPEG Baseline and 8-bit Extended through libjpeg-turbo; YBR_FULL and YBR_FULL_422 to RGB |
| Writing | Implicit/explicit little endian, explicit big endian, deflated; compressed pixel data written back unchanged; `dcm_set_pixel_data()` stores an R array |
| Ultrasound | Sequence of Ultrasound Regions as a tidy table, pixel to physical conversion, distances, frame times |
| Checking | `dcm_compare_pydicom(file)` runs pydicom on your own file and reports any difference |

Arrays follow pydicom's axis order, so `arr[f, r, c, s]` in R is
`arr[f-1, r-1, c-1, s-1]` in Python. By default JPEG is decoded to match
pydicom with its Pillow plugin bit for bit, so a model trained on pydicom
output sees the same pixels from R.

## Validation against pydicom

`inst/validation/compare_all.R` runs every sample file shipped with pydicom
3.0.2 (78 files) through both packages:

- Data elements: 6,733 compared, including nested sequence items, 0 differences.
- Pixel arrays: identical on all 44 files that both decode, including the
  30-frame YBR_FULL_422 JPEG echo clip, all 14 JPEG Baseline files, RLE
  8/16/32-bit, big endian and deflated files. The 12 files only pydicom
  decodes are JPEG 2000, JPEG-LS and JPEG Lossless.
- Speed: an 80-frame 600 x 800 colour JPEG clip decodes in about 2.5 s,
  in line with pydicom with Pillow (about 2 s).
- Round trips: every sample written back out in each uncompressed transfer
  syntax reads back in pydicom with equal elements and identical pixels
  (only expected differences: byte-swapped binary values after an endian
  change, and odd-length values padded to even length).

## Not yet supported

- JPEG 2000, JPEG-LS, JPEG Lossless (processes 14 and SV1) and 12-bit JPEG
  Extended decoding. These files still read; `dcm_pixel_array()` reports
  that the codec isn't available.
- Re-encoding (compressing) pixel data on write.
- Character sets beyond Latin-1 and UTF-8 (ISO 2022 escapes, Japanese,
  Chinese, Korean, Cyrillic and so on); such text is read as Latin-1.
- Palette colour, VOI LUT and windowing helpers, overlays as arrays,
  DICOMDIR, networking.
- Pixel arrays are R integers (4 bytes per sample), so a long colour clip
  takes about four times pydicom's memory; use `frames =` to decode part.

## Licence

MIT. The data dictionary is generated from pydicom's (MIT); see
`LICENSE.note`.
