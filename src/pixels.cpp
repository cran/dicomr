// Pixel data decoding: native (uncompressed), RLE Lossless and JPEG baseline
// via libjpeg-turbo, plus YBR_FULL / YBR_FULL_422 handling.
//
// Every decoder writes into an R array of dim (frames, rows, cols, samples),
// matching the axis order of pydicom's pixel_array.

#include <Rcpp.h>
#include <cstdio>
#include <cstring>
#include <csetjmp>
#include <cmath>
#include <vector>
#include <string>
#include <jpeglib.h>

using namespace Rcpp;

namespace {

// Output buffer: frames x rows x cols x spp, R column-major.
struct Out {
  int nf, rows, cols, spp;
  inline R_xlen_t idx(int f, int r, int c, int s) const {
    return (R_xlen_t) f + (R_xlen_t) nf * (r + (R_xlen_t) rows * (c + (R_xlen_t) cols * s));
  }
};

// pydicom's float32 YBR_FULL -> RGB conversion (PS3.3 C.7.6.3.1.2)
inline void ybr_to_rgb(float y, float cb, float cr, int& r, int& g, int& b) {
  const float m10 = (float) (-0.114 * 1.772 / 0.587);
  const float m11 = (float) 1.772;
  const float m20 = (float) 1.402;
  const float m21 = (float) (-0.299 * 1.402 / 0.587);
  cb -= 128.0f; cr -= 128.0f;
  float fr = y + cr * m20;
  float fg = y + cb * m10 + cr * m21;
  float fb = y + cb * m11;
  auto conv = [](float v) {
    v = std::floor(v + 0.5f);
    if (v < 0) v = 0;
    if (v > 255) v = 255;
    return (int) v;
  };
  r = conv(fr); g = conv(fg); b = conv(fb);
}

inline double read_sample(const uint8_t* p, int bytes, bool is_signed, bool little) {
  uint64_t u = 0;
  if (little) for (int k = bytes - 1; k >= 0; --k) u = (u << 8) | p[k];
  else for (int k = 0; k < bytes; ++k) u = (u << 8) | p[k];
  if (!is_signed) return (double) u;
  switch (bytes) {
    case 1: return (double) (int8_t) u;
    case 2: return (double) (int16_t) u;
    case 4: return (double) (int32_t) u;
    default: return (double) (int64_t) u;
  }
}

SEXP alloc_out(const Out& o, bool as_double) {
  R_xlen_t n = (R_xlen_t) o.nf * o.rows * o.cols * o.spp;
  SEXP arr = PROTECT(Rf_allocVector(as_double ? REALSXP : INTSXP, n));
  SEXP dim = PROTECT(Rf_allocVector(INTSXP, 4));
  INTEGER(dim)[0] = o.nf; INTEGER(dim)[1] = o.rows; INTEGER(dim)[2] = o.cols; INTEGER(dim)[3] = o.spp;
  Rf_setAttrib(arr, R_DimSymbol, dim);
  UNPROTECT(2);
  return arr;
}

// Store one frame held as interleaved little-endian samples in `buf`.
void store_frame(SEXP arr, const Out& o, int f, const uint8_t* buf, int bytes,
                 bool is_signed, bool little, bool planar, bool as_double,
                 int bits_stored, bool correct_bits) {
  int shift = bytes * 8 - bits_stored;
  bool do_shift = correct_bits && shift > 0 && bytes <= 4;
  for (int r = 0; r < o.rows; ++r) {
    for (int c = 0; c < o.cols; ++c) {
      for (int s = 0; s < o.spp; ++s) {
        size_t k = planar ? ((size_t) s * o.rows * o.cols + (size_t) r * o.cols + c)
                          : (((size_t) r * o.cols + c) * o.spp + s);
        double v = read_sample(buf + k * bytes, bytes, false, little);
        uint64_t u = (uint64_t) v;
        if (do_shift) {
          // emulate pydicom: left shift then (arithmetic) right shift
          int width = bytes * 8;
          uint64_t mask = (width == 64) ? ~0ULL : ((1ULL << width) - 1);
          u = (u << shift) & mask;
          if (is_signed) {
            int64_t sv;
            if (width == 8) sv = (int8_t) u; else if (width == 16) sv = (int16_t) u; else sv = (int32_t) u;
            v = (double) (sv >> shift);
          } else {
            v = (double) (u >> shift);
          }
        } else {
          v = read_sample(buf + k * bytes, bytes, is_signed, little);
        }
        R_xlen_t i = o.idx(f, r, c, s);
        if (as_double) REAL(arr)[i] = v; else INTEGER(arr)[i] = (int) v;
      }
    }
  }
}

// ---------------- JPEG (libjpeg-turbo) ----------------

struct JErr {
  struct jpeg_error_mgr pub;
  jmp_buf jb;
  char msg[JMSG_LENGTH_MAX];
};

static void jerr_exit(j_common_ptr cinfo) {
  JErr* e = (JErr*) cinfo->err;
  (*cinfo->err->format_message)(cinfo, e->msg);
  longjmp(e->jb, 1);
}

static void jerr_silent(j_common_ptr, int) {}

struct JpegResult {
  std::vector<uint8_t> data;
  int width = 0, height = 0, comps = 0;
  bool rgb_ids = false;
  bool jfif = false;
  bool adobe = false;
  bool raw_ycbcr = false;  // output left in YCbCr
  std::string error;
};

// Decode one JPEG baseline/extended (8-bit) codestream. Colour handling
// mirrors pydicom + Pillow: with no Adobe APP14 marker the samples are
// returned without colour transformation.
//
// libjpeg reports fatal errors through error_exit, which here longjmps back
// to the setjmp below. The function that calls setjmp keeps no C++ objects
// of its own: results go to `res`, which lives in the caller's frame, and
// the error message is copied into it only after the jump.
static void decode_jpeg_into(const uint8_t* src, size_t n, bool convert_in_libjpeg,
                             JpegResult& res) {
  struct jpeg_decompress_struct cinfo;
  JErr jerr;
  cinfo.err = jpeg_std_error(&jerr.pub);
  jerr.pub.error_exit = jerr_exit;
  jerr.pub.emit_message = jerr_silent;
  jerr.msg[0] = '\0';
  if (setjmp(jerr.jb)) {
    jpeg_destroy_decompress(&cinfo);
    res.error = jerr.msg[0] ? jerr.msg : "JPEG decoding failed";
    return;
  }
  jpeg_create_decompress(&cinfo);
  // The cast suits both the older non-const and newer const prototypes.
  jpeg_mem_src(&cinfo, (unsigned char*) src, (unsigned long) n);
  jpeg_read_header(&cinfo, TRUE);
  if (cinfo.data_precision != 8) {
    int prec = cinfo.data_precision;
    jpeg_destroy_decompress(&cinfo);
    res.error = "Only 8-bit JPEG is supported (this codestream has " +
                std::to_string(prec) + "-bit precision)";
    return;
  }
  res.jfif = cinfo.saw_JFIF_marker;
  res.adobe = cinfo.saw_Adobe_marker;
  if (cinfo.num_components == 3) {
    int a = cinfo.comp_info[0].component_id, b = cinfo.comp_info[1].component_id,
        c = cinfo.comp_info[2].component_id;
    res.rgb_ids = (a == 82 && b == 71 && c == 66) || (a == 114 && b == 103 && c == 98);
    if (convert_in_libjpeg || cinfo.saw_Adobe_marker) {
      cinfo.out_color_space = JCS_RGB;
      res.raw_ycbcr = false;
    } else {
      cinfo.out_color_space = cinfo.jpeg_color_space;  // no transform
      res.raw_ycbcr = (cinfo.jpeg_color_space == JCS_YCbCr);
    }
  }
  jpeg_start_decompress(&cinfo);
  res.width = (int) cinfo.output_width;
  res.height = (int) cinfo.output_height;
  res.comps = cinfo.output_components;
  size_t stride = (size_t) res.width * (size_t) res.comps;
  res.data.resize(stride * (size_t) res.height);
  uint8_t* base = res.data.data();
  while (cinfo.output_scanline < cinfo.output_height) {
    JSAMPROW row = base + stride * cinfo.output_scanline;
    jpeg_read_scanlines(&cinfo, &row, 1);
  }
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
}

JpegResult decode_jpeg(const uint8_t* src, size_t n, bool convert_in_libjpeg) {
  JpegResult res;
  decode_jpeg_into(src, n, convert_in_libjpeg, res);
  if (!res.error.empty()) res.data.clear();
  return res;
}

// ---------------- RLE Lossless (PS3.5 Annex G) ----------------

void rle_segment(const uint8_t* src, size_t n, std::vector<uint8_t>& out, size_t want) {
  out.clear();
  out.reserve(want);
  size_t i = 0;
  while (i < n && out.size() < want) {
    int8_t h = (int8_t) src[i++];
    if (h >= 0) {
      size_t cnt = (size_t) h + 1;
      if (i + cnt > n) cnt = n - i;
      out.insert(out.end(), src + i, src + i + cnt);
      i += cnt;
    } else if (h != -128) {
      size_t cnt = (size_t) (1 - h);
      if (i >= n) break;
      out.insert(out.end(), cnt, src[i++]);
    }
  }
  if (out.size() > want) out.resize(want);
}

}  // namespace

// [[Rcpp::export(.decode_native)]]
SEXP decode_native(RawVector data, int rows, int cols, int spp, int bits_alloc,
                   int bits_stored, int pixel_rep, int planar, int nframes,
                   IntegerVector frames, bool little, bool swap8, bool is_float,
                   bool correct_bits) {
  Out o{(int) frames.size(), rows, cols, spp};
  const uint8_t* b = data.begin();
  size_t n = data.size();
  size_t px = (size_t) rows * cols * spp;

  if (bits_alloc == 1) {
    SEXP arr = PROTECT(alloc_out(o, false));
    for (int fi = 0; fi < o.nf; ++fi) {
      size_t bit0 = (size_t) frames[fi] * px;
      for (int r = 0; r < rows; ++r)
        for (int c = 0; c < cols; ++c)
          for (int s = 0; s < spp; ++s) {
            size_t bit = bit0 + ((size_t) r * cols + c) * spp + s;
            if (bit / 8 >= n) Rcpp::stop("Pixel data is shorter than expected");
            INTEGER(arr)[o.idx(fi, r, c, s)] = (b[bit / 8] >> (bit % 8)) & 1;
          }
    }
    UNPROTECT(1);
    return arr;
  }

  int bytes = bits_alloc / 8;
  size_t fbytes = px * bytes;
  bool as_double = is_float || (bytes == 4 && pixel_rep == 0) || bytes == 8;
  SEXP arr = PROTECT(alloc_out(o, as_double));
  std::vector<uint8_t> tmp;
  for (int fi = 0; fi < o.nf; ++fi) {
    size_t off = (size_t) frames[fi] * fbytes;
    if (off + fbytes > n) Rcpp::stop("Pixel data is shorter than expected for frame %d", frames[fi] + 1);
    const uint8_t* fb = b + off;
    if (swap8 && bytes == 1) {
      // 8-bit data stored as OW in big endian: bytes are swapped in pairs
      tmp.assign(b, b + n + (n & 1));
      if (n & 1) tmp[n] = 0;
      for (size_t i = 0; i + 1 < tmp.size(); i += 2) std::swap(tmp[i], tmp[i + 1]);
      fb = tmp.data() + off;
    }
    if (is_float) {
      for (int r = 0; r < rows; ++r)
        for (int c = 0; c < cols; ++c)
          for (int s = 0; s < spp; ++s) {
            size_t k = planar ? ((size_t) s * rows * cols + (size_t) r * cols + c)
                              : (((size_t) r * cols + c) * spp + s);
            uint8_t w[8];
            for (int q = 0; q < bytes; ++q) w[q] = little ? fb[k * bytes + q] : fb[k * bytes + bytes - 1 - q];
            double v;
            if (bytes == 4) { float f; std::memcpy(&f, w, 4); v = f; } else std::memcpy(&v, w, 8);
            REAL(arr)[o.idx(fi, r, c, s)] = v;
          }
    } else {
      store_frame(arr, o, fi, fb, bytes, pixel_rep == 1, little, planar == 1,
                  as_double, bits_stored, correct_bits);
    }
  }
  UNPROTECT(1);
  return arr;
}

// Expand YBR_FULL_422 native data (Y1 Y2 Cb Cr per pixel pair) to full YBR.
// [[Rcpp::export(.decode_ybr422)]]
SEXP decode_ybr422(RawVector data, int rows, int cols, IntegerVector frames) {
  Out o{(int) frames.size(), rows, cols, 3};
  size_t fbytes = (size_t) rows * cols * 2;
  SEXP arr = PROTECT(alloc_out(o, false));
  int* a = INTEGER(arr);
  for (int fi = 0; fi < o.nf; ++fi) {
    size_t off = (size_t) frames[fi] * fbytes;
    if (off + fbytes > (size_t) data.size()) Rcpp::stop("Pixel data is shorter than expected");
    const uint8_t* p = data.begin() + off;
    for (int r = 0; r < rows; ++r) {
      for (int c = 0; c < cols; c += 2) {
        size_t k = ((size_t) r * cols + c) * 2;
        int y1 = p[k], y2 = p[k + 1], cb = p[k + 2], cr = p[k + 3];
        a[o.idx(fi, r, c, 0)] = y1; a[o.idx(fi, r, c, 1)] = cb; a[o.idx(fi, r, c, 2)] = cr;
        if (c + 1 < cols) {
          a[o.idx(fi, r, c + 1, 0)] = y2; a[o.idx(fi, r, c + 1, 1)] = cb; a[o.idx(fi, r, c + 1, 2)] = cr;
        }
      }
    }
  }
  UNPROTECT(1);
  return arr;
}

// In-place YBR_FULL -> RGB on an integer array of dim (frames, rows, cols, 3).
// [[Rcpp::export(.ybr_to_rgb)]]
IntegerVector ybr_to_rgb_array(IntegerVector arr) {
  IntegerVector dim = arr.attr("dim");
  R_xlen_t plane = (R_xlen_t) dim[0] * dim[1] * dim[2];
  IntegerVector out = clone(arr);
  int* a = out.begin();
  for (R_xlen_t i = 0; i < plane; ++i) {
    int r, g, b;
    ybr_to_rgb((float) a[i], (float) a[i + plane], (float) a[i + 2 * plane], r, g, b);
    a[i] = r; a[i + plane] = g; a[i + 2 * plane] = b;
  }
  return out;
}

// Decode a list of JPEG frames (raw vectors) into one array. `pi` is the
// data set's Photometric Interpretation; it is conformed to the codestream
// the way pydicom does (RGB component IDs mean RGB, a JFIF marker means YBR),
// and YBR output is converted to RGB with pydicom's formula when `as_rgb`.
// [[Rcpp::export(.decode_jpeg_frames)]]
SEXP decode_jpeg_frames(List frames, int rows, int cols, int spp, std::string pi,
                        bool baseline, bool as_rgb, bool accurate) {
  Out o{(int) frames.size(), rows, cols, spp};
  SEXP arr = PROTECT(alloc_out(o, false));
  int* a = INTEGER(arr);
  bool convert = false;
  size_t npx = (size_t) rows * cols, fsize = npx * spp;
  std::vector<uint8_t> all((size_t) o.nf * fsize);  // frame-contiguous decode buffer
  for (int fi = 0; fi < o.nf; ++fi) {
    RawVector src = frames[fi];
    JpegResult jr = decode_jpeg(src.begin(), src.size(), accurate && as_rgb);
    if (!jr.error.empty()) {
      UNPROTECT(1);
      Rcpp::stop("JPEG decoding failed for frame %d: %s", fi + 1, jr.error.c_str());
    }
    if (jr.width != cols || jr.height != rows || jr.comps != spp) {
      UNPROTECT(1);
      Rcpp::stop("JPEG frame %d is %dx%dx%d but the data set says %dx%dx%d", fi + 1,
                 jr.height, jr.width, jr.comps, rows, cols, spp);
    }
    if (fi == 0 && spp == 3) {
      if (jr.rgb_ids && pi != "RGB") pi = "RGB";
      else if (jr.jfif && pi.find("YBR") == std::string::npos) pi = baseline ? "YBR_FULL_422" : "YBR_FULL";
      if (!jr.raw_ycbcr) pi = "RGB";  // decoder already produced RGB
      convert = as_rgb && (pi == "YBR_FULL" || pi == "YBR_FULL_422");
      if (convert) pi = "RGB";
    }
    uint8_t* p = jr.data.data();
    if (convert) {
      for (size_t k = 0; k < npx; ++k) {
        int r, g, b;
        ybr_to_rgb((float) p[3 * k], (float) p[3 * k + 1], (float) p[3 * k + 2], r, g, b);
        p[3 * k] = (uint8_t) r; p[3 * k + 1] = (uint8_t) g; p[3 * k + 2] = (uint8_t) b;
      }
    }
    std::memcpy(all.data() + (size_t) fi * fsize, p, fsize);
  }
  // Transpose to (frames, rows, cols, samples), writing the output in order
  R_xlen_t i = 0;
  const uint8_t* base = all.data();
  for (int s = 0; s < spp; ++s)
    for (int c = 0; c < cols; ++c)
      for (int r = 0; r < rows; ++r) {
        const uint8_t* q = base + ((size_t) r * cols + c) * spp + s;
        for (int fi = 0; fi < o.nf; ++fi) a[i++] = q[(size_t) fi * fsize];
      }
  Rf_setAttrib(arr, Rf_install("photometric_interpretation"), Rf_mkString(pi.c_str()));
  UNPROTECT(1);
  return arr;
}

// [[Rcpp::export(.decode_rle_frames)]]
SEXP decode_rle_frames(List frames, int rows, int cols, int spp, int bits_alloc,
                       int bits_stored, int pixel_rep, bool correct_bits) {
  Out o{(int) frames.size(), rows, cols, spp};
  int bytes = bits_alloc / 8;
  if (bytes < 1 || bytes > 4 || bits_alloc % 8) Rcpp::stop("RLE: unsupported Bits Allocated %d", bits_alloc);
  bool as_double = (bytes == 4 && pixel_rep == 0);
  SEXP arr = PROTECT(alloc_out(o, as_double));
  size_t px = (size_t) rows * cols;
  std::vector<uint8_t> frame(px * spp * bytes), seg;
  for (int fi = 0; fi < o.nf; ++fi) {
    RawVector src = frames[fi];
    const uint8_t* s = src.begin();
    size_t n = src.size();
    if (n < 64) { UNPROTECT(1); Rcpp::stop("RLE frame %d has no valid header", fi + 1); }
    auto u32 = [&](size_t p) {
      return (uint32_t) s[p] | ((uint32_t) s[p + 1] << 8) | ((uint32_t) s[p + 2] << 16) | ((uint32_t) s[p + 3] << 24);
    };
    uint32_t nseg = u32(0);
    if ((int) nseg != spp * bytes) {
      UNPROTECT(1);
      Rcpp::stop("RLE frame %d has %d segments, expected %d", fi + 1, (int) nseg, spp * bytes);
    }
    std::fill(frame.begin(), frame.end(), 0);
    for (uint32_t k = 0; k < nseg; ++k) {
      size_t start = u32(4 + 4 * k);
      size_t stop = (k + 1 < nseg) ? u32(4 + 4 * (k + 1)) : n;
      if (start > n || stop > n || stop < start) { UNPROTECT(1); Rcpp::stop("RLE segment offsets are invalid"); }
      rle_segment(s + start, stop - start, seg, px);
      int sample = k / bytes;
      int byte_in_sample = bytes - 1 - (k % bytes);  // first segment is the MSB
      for (size_t i = 0; i < seg.size(); ++i)
        frame[(i * spp + sample) * bytes + byte_in_sample] = seg[i];
    }
    store_frame(arr, o, fi, frame.data(), bytes, pixel_rep == 1, true, false, as_double,
                bits_stored, correct_bits);
  }
  UNPROTECT(1);
  return arr;
}

// Find JPEG frame boundaries: a new frame starts at a fragment that begins
// with the SOI marker (FF D8).
// [[Rcpp::export(.fragment_starts_soi)]]
LogicalVector fragment_starts_soi(List fragments) {
  LogicalVector out(fragments.size());
  for (R_xlen_t i = 0; i < fragments.size(); ++i) {
    RawVector f = fragments[i];
    out[i] = f.size() >= 2 && f[0] == 0xFF && f[1] == 0xD8;
  }
  return out;
}
