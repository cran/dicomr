// DICOM Part 10 file parser.
//
// Handles the 128-byte preamble and "DICM" prefix (or files without them),
// the File Meta Information group, implicit and explicit VR, little and big
// endian, deflated explicit VR little endian, sequences with defined and
// undefined lengths, UN elements that hold undefined-length sequences, and
// encapsulated (compressed) pixel data split into fragments.

#include <Rcpp.h>
#include <zlib.h>
#include <string>
#include <vector>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <cctype>
#include "dicomr.h"

using namespace Rcpp;

namespace {

const uint32_t TAG_ITEM = 0xFFFEE000;
const uint32_t TAG_ITEM_DELIM = 0xFFFEE00D;
const uint32_t TAG_SEQ_DELIM = 0xFFFEE0DD;
const uint32_t UNDEFINED = 0xFFFFFFFF;
const uint32_t TAG_PIXEL_DATA = 0x7FE00010;

enum Charset { CS_DEFAULT = 0, CS_LATIN1 = 1, CS_UTF8 = 2, CS_OTHER = 3 };

bool is_long_vr(const std::string& vr) {
  static const char* L[] = {"OB", "OD", "OF", "OL", "OV", "OW", "SQ", "UC",
                            "UN", "UR", "UT", "SV", "UV"};
  for (const char* v : L) if (vr == v) return true;
  return false;
}

bool is_known_vr(const std::string& vr) {
  static const char* K[] = {"AE", "AS", "AT", "CS", "DA", "DS", "DT", "FL", "FD",
                            "IS", "LO", "LT", "OB", "OD", "OF", "OL", "OV", "OW",
                            "PN", "SH", "SL", "SQ", "SS", "ST", "SV", "TM", "UC",
                            "UI", "UL", "UN", "UR", "US", "UT", "UV"};
  for (const char* v : K) if (vr == v) return true;
  return false;
}

std::string hex8(uint32_t tag) {
  char b[9];
  std::snprintf(b, sizeof(b), "%08X", tag);
  return std::string(b);
}

struct Ctx {
  const uint8_t* buf;
  size_t n;
  bool little;
  bool implicit;
  bool encapsulated;     // transfer syntax uses encapsulated pixel data
  bool stop_before_pixels;
  int pixel_rep;         // for resolving "US or SS"; -1 unknown
  std::vector<std::string> warnings;
};

uint16_t rd16(const Ctx& c, size_t p) {
  if (p + 2 > c.n) Rcpp::stop("Unexpected end of file at byte %d", (int) p);
  const uint8_t* b = c.buf + p;
  return c.little ? (uint16_t)(b[0] | (b[1] << 8)) : (uint16_t)((b[0] << 8) | b[1]);
}

uint32_t rd32(const Ctx& c, size_t p) {
  if (p + 4 > c.n) Rcpp::stop("Unexpected end of file at byte %d", (int) p);
  const uint8_t* b = c.buf + p;
  if (c.little) return (uint32_t) b[0] | ((uint32_t) b[1] << 8) | ((uint32_t) b[2] << 16) | ((uint32_t) b[3] << 24);
  return ((uint32_t) b[0] << 24) | ((uint32_t) b[1] << 16) | ((uint32_t) b[2] << 8) | (uint32_t) b[3];
}

uint32_t rdtag(const Ctx& c, size_t p) {
  return ((uint32_t) rd16(c, p) << 16) | rd16(c, p + 2);
}

// Copy `len` bytes in file order, swapping each `width`-byte word when the
// data is big endian, so the result is always little endian.
std::vector<uint8_t> le_bytes(const Ctx& c, size_t p, size_t len, int width) {
  std::vector<uint8_t> v(c.buf + p, c.buf + p + len);
  if (!c.little && width > 1) {
    for (size_t i = 0; i + width <= len; i += width)
      for (int a = 0, b = width - 1; a < b; ++a, --b) std::swap(v[i + a], v[i + b]);
  }
  return v;
}

SEXP mkstr(const std::string& s, int cs) {
  bool ascii = true;
  for (unsigned char ch : s) if (ch > 127) { ascii = false; break; }
  cetype_t enc = CE_NATIVE;
  if (!ascii) enc = (cs == CS_UTF8) ? CE_UTF8 : CE_LATIN1;
  return Rf_mkCharLenCE(s.data(), (int) s.size(), enc);
}

std::string strip(const std::string& s, const char* chars) {
  size_t b = s.find_first_not_of(chars);
  if (b == std::string::npos) return std::string();
  size_t e = s.find_last_not_of(chars);
  return s.substr(b, e - b + 1);
}

std::vector<std::string> split_bs(const std::string& s) {
  std::vector<std::string> out;
  size_t start = 0;
  while (true) {
    size_t k = s.find('\\', start);
    if (k == std::string::npos) { out.push_back(s.substr(start)); break; }
    out.push_back(s.substr(start, k - start));
    start = k + 1;
  }
  return out;
}

SEXP strings_sexp(const std::vector<std::string>& v, int cs) {
  SEXP out = PROTECT(Rf_allocVector(STRSXP, v.size()));
  for (size_t i = 0; i < v.size(); ++i) SET_STRING_ELT(out, i, mkstr(v[i], cs));
  UNPROTECT(1);
  return out;
}

// Decode the value of a non-sequence element into an R vector.
SEXP decode_value(Ctx& c, const std::string& vr, size_t p, size_t len, int cs) {
  if (len == 0) return R_NilValue;
  std::string raw((const char*) c.buf + p, len);

  if (vr == "AE") {
    std::vector<std::string> parts = split_bs(raw);
    for (auto& s : parts) s = strip(s, " ");
    return strings_sexp(parts, cs);
  }
  if (vr == "AS" || vr == "CS" || vr == "DA" || vr == "DT" || vr == "TM" ||
      vr == "PN" || vr == "UI") {
    size_t e = raw.size();
    while (e > 0 && (raw[e - 1] == ' ' || raw[e - 1] == '\0')) --e;
    std::string s = raw.substr(0, e);
    return strings_sexp(split_bs(s), cs);
  }
  if (vr == "SH" || vr == "LO" || vr == "UC") {
    std::vector<std::string> parts = split_bs(raw);
    for (auto& s : parts) {
      size_t e = s.size();
      while (e > 0 && (s[e - 1] == ' ' || s[e - 1] == '\0')) --e;
      s = s.substr(0, e);
    }
    return strings_sexp(parts, cs);
  }
  if (vr == "LT" || vr == "ST" || vr == "UT" || vr == "UR") {
    size_t e = raw.size();
    while (e > 0 && (raw[e - 1] == ' ' || raw[e - 1] == '\0')) --e;
    std::vector<std::string> one(1, raw.substr(0, e));
    return strings_sexp(one, cs);
  }
  if (vr == "IS" || vr == "DS") {
    // strip whole value, split, parse each; keep as character if any fail
    size_t b = 0, e = raw.size();
    while (b < e && (raw[b] == ' ' || raw[b] == '\0')) ++b;
    while (e > b && (raw[e - 1] == ' ' || raw[e - 1] == '\0')) --e;
    std::string s = raw.substr(b, e - b);
    if (s.empty()) return R_NilValue;
    std::vector<std::string> parts = split_bs(s);
    std::vector<double> vals(parts.size());
    bool ok = true, all_int = true;
    for (size_t i = 0; i < parts.size(); ++i) {
      std::string t = strip(parts[i], " ");
      if (t.empty()) { vals[i] = NA_REAL; continue; }
      char* endp = nullptr;
      double d = std::strtod(t.c_str(), &endp);
      if (endp == t.c_str() || *endp != '\0') { ok = false; break; }
      vals[i] = d;
      if (vr == "IS" && (d != std::floor(d) || std::fabs(d) > 2147483647.0)) all_int = false;
    }
    if (!ok) {
      c.warnings.push_back("Invalid " + vr + " value '" + s + "' kept as character");
      return strings_sexp(parts, cs);
    }
    if (vr == "IS" && all_int) {
      IntegerVector out(vals.size());
      for (size_t i = 0; i < vals.size(); ++i)
        out[i] = ISNA(vals[i]) ? NA_INTEGER : (int) vals[i];
      return out;
    }
    return NumericVector(vals.begin(), vals.end());
  }
  if (vr == "US" || vr == "SS") {
    size_t m = len / 2;
    IntegerVector out(m);
    for (size_t i = 0; i < m; ++i) {
      uint16_t u = rd16(c, p + 2 * i);
      out[i] = (vr == "US") ? (int) u : (int) (int16_t) u;
    }
    return out;
  }
  if (vr == "UL" || vr == "SL") {
    size_t m = len / 4;
    if (vr == "SL") {
      IntegerVector out(m);
      for (size_t i = 0; i < m; ++i) out[i] = (int) (int32_t) rd32(c, p + 4 * i);
      return out;
    }
    NumericVector out(m);
    for (size_t i = 0; i < m; ++i) out[i] = (double) rd32(c, p + 4 * i);
    return out;
  }
  if (vr == "SV" || vr == "UV") {
    size_t m = len / 8;
    NumericVector out(m);
    std::vector<uint8_t> v = le_bytes(c, p, m * 8, 8);
    for (size_t i = 0; i < m; ++i) {
      uint64_t u = 0;
      for (int k = 7; k >= 0; --k) u = (u << 8) | v[8 * i + k];
      out[i] = (vr == "UV") ? (double) u : (double) (int64_t) u;
    }
    return out;
  }
  if (vr == "FL") {
    size_t m = len / 4;
    NumericVector out(m);
    std::vector<uint8_t> v = le_bytes(c, p, m * 4, 4);
    for (size_t i = 0; i < m; ++i) {
      float f; std::memcpy(&f, &v[4 * i], 4);  // host assumed little endian
      out[i] = (double) f;
    }
    return out;
  }
  if (vr == "FD") {
    size_t m = len / 8;
    NumericVector out(m);
    std::vector<uint8_t> v = le_bytes(c, p, m * 8, 8);
    for (size_t i = 0; i < m; ++i) {
      double d; std::memcpy(&d, &v[8 * i], 8);
      out[i] = d;
    }
    return out;
  }
  if (vr == "AT") {
    size_t m = len / 4;
    std::vector<std::string> v(m);
    for (size_t i = 0; i < m; ++i) v[i] = hex8(rdtag(c, p + 4 * i));
    return strings_sexp(v, CS_DEFAULT);
  }
  // OB, OW, OF, OD, OL, OV, UN and anything else: raw bytes as stored
  RawVector out(len);
  std::memcpy(&out[0], c.buf + p, len);
  return out;
}

int charset_code(SEXP val) {
  if (Rf_isNull(val) || TYPEOF(val) != STRSXP || Rf_length(val) == 0) return CS_DEFAULT;
  int code = CS_DEFAULT;
  for (int i = 0; i < Rf_length(val); ++i) {
    std::string s = CHAR(STRING_ELT(val, i));
    if (s.empty() || s == "ISO_IR 6" || s == "ISO 2022 IR 6") continue;
    if (s == "ISO_IR 192") return CS_UTF8;
    if (s == "ISO_IR 100" || s == "ISO 2022 IR 100") code = CS_LATIN1;
    else code = CS_OTHER;
  }
  return code;
}

List parse_dataset(Ctx& c, size_t& pos, size_t end, int cs, bool stop_at_meta_end);

List make_element(uint32_t tag, const std::string& vr, SEXP value, bool undefined_length) {
  List el = List::create(_["group"] = (int) (tag >> 16), _["element"] = (int) (tag & 0xFFFF),
                         _["vr"] = vr, _["value"] = value);
  if (undefined_length) el.attr("undefined_length") = true;
  el.attr("class") = "dicom_element";
  return el;
}

// Parse a sequence's items, from pos to end (end may be SIZE_MAX for undefined).
List parse_sequence(Ctx& c, size_t& pos, size_t end, int cs) {
  std::vector<List> keep;
  while (pos < end && pos + 8 <= c.n) {
    uint32_t tag = rdtag(c, pos);
    uint32_t len = rd32(c, pos + 4);
    pos += 8;
    if (tag == TAG_SEQ_DELIM) break;
    if (tag != TAG_ITEM) {
      c.warnings.push_back("Expected an item tag in a sequence but found " + hex8(tag));
      // skip unexpected element
      if (len != UNDEFINED) pos += len;
      continue;
    }
    size_t item_end = (len == UNDEFINED) ? (size_t) -1 : pos + len;
    if (item_end != (size_t) -1 && item_end > c.n) {
      c.warnings.push_back("Sequence item runs past end of file");
      item_end = c.n;
    }
    List ds = parse_dataset(c, pos, item_end, cs, false);
    keep.push_back(ds);
    if (len != UNDEFINED) pos = item_end;
  }
  List seq(keep.size());
  for (size_t i = 0; i < keep.size(); ++i) seq[i] = keep[i];
  seq.attr("class") = "dicom_sequence";
  return seq;
}

// Encapsulated pixel data: first item is the Basic Offset Table.
List parse_encapsulated(Ctx& c, size_t& pos) {
  std::vector<RawVector> frags;
  NumericVector offsets;
  bool first = true;
  while (pos + 8 <= c.n) {
    uint32_t tag = rdtag(c, pos);
    uint32_t len = rd32(c, pos + 4);
    pos += 8;
    if (tag == TAG_SEQ_DELIM) break;
    if (tag != TAG_ITEM || len == UNDEFINED) {
      c.warnings.push_back("Malformed encapsulated pixel data at tag " + hex8(tag));
      break;
    }
    if (pos + len > c.n) {
      c.warnings.push_back("Pixel data fragment runs past end of file; truncated");
      len = (uint32_t) (c.n - pos);
    }
    if (first) {
      size_t m = len / 4;
      offsets = NumericVector(m);
      for (size_t i = 0; i < m; ++i) offsets[i] = (double) rd32(c, pos + 4 * i);
      first = false;
    } else {
      RawVector f(len);
      if (len) std::memcpy(&f[0], c.buf + pos, len);
      frags.push_back(f);
    }
    pos += len;
  }
  List fl(frags.size());
  for (size_t i = 0; i < frags.size(); ++i) fl[i] = frags[i];
  List out = List::create(_["offsets"] = offsets, _["fragments"] = fl);
  out.attr("class") = "dicom_encapsulated";
  return out;
}

List parse_dataset(Ctx& c, size_t& pos, size_t end, int cs, bool stop_at_meta_end) {
  std::vector<List> elems;
  std::vector<std::string> names;
  int local_cs = cs;

  while (pos < end && pos + 4 <= c.n) {
    if (pos + 8 > c.n) {  // trailing junk shorter than a header
      c.warnings.push_back("Trailing bytes at end of data set ignored");
      pos = c.n;
      break;
    }
    uint32_t tag = rdtag(c, pos);
    uint16_t group = tag >> 16;
    if (stop_at_meta_end && group != 0x0002) break;

    if (tag == TAG_ITEM_DELIM) { pos += 8; break; }  // end of undefined-length item
    if (tag == TAG_SEQ_DELIM) { pos += 8; break; }   // stray; end this data set
    if (group == 0xFFFE) {                            // stray item tag
      uint32_t len = rd32(c, pos + 4);
      c.warnings.push_back("Unexpected item tag in data set skipped");
      pos += 8 + (len == UNDEFINED ? 0 : len);
      continue;
    }

    std::string vr;
    uint32_t len;
    size_t hdr = pos + 4;
    if (!c.implicit) {
      if (hdr + 2 > c.n) Rcpp::stop("Unexpected end of file");
      vr = std::string((const char*) c.buf + hdr, 2);
      if (is_known_vr(vr)) {
        if (is_long_vr(vr)) { len = rd32(c, hdr + 4); hdr += 8; }
        else { len = rd16(c, hdr + 2); hdr += 4; }
      } else {
        bool letters = std::isupper((unsigned char) vr[0]) && std::isupper((unsigned char) vr[1]);
        if (letters) {
          // An unrecognised VR: read with a 2-byte length and keep the bytes
          c.warnings.push_back("Unknown VR '" + vr + "' for tag " + hex8(tag) + "; value kept as UN");
          len = rd16(c, hdr + 2);
          hdr += 4;
          vr = "UN";
        } else {
          // Element encoded with implicit VR inside an explicit VR data set
          c.warnings.push_back("Tag " + hex8(tag) + " has no explicit VR; read as implicit VR");
          vr = dict_vr(tag);
          len = rd32(c, hdr);
          hdr += 4;
        }
      }
    } else {
      vr = dict_vr(tag);
      len = rd32(c, hdr);
      hdr += 4;
    }

    // Resolve ambiguous dictionary VRs (implicit VR, or explicit fallback)
    if (vr.find(" or ") != std::string::npos) {
      if (tag == TAG_PIXEL_DATA) vr = c.encapsulated ? "OB" : "OW";
      else if (vr == "OB or OW") vr = "OW";
      else if (vr == "US or SS or OW" || vr == "US or OW") vr = "US";
      else if (vr == "US or SS") vr = (c.pixel_rep == 1) ? "SS" : "US";
      else vr = vr.substr(0, 2);
    }
    if (vr == "NONE") vr = "UN";

    pos = hdr;
    if (tag == TAG_PIXEL_DATA && c.stop_before_pixels) { pos = c.n; break; }

    SEXP value;
    bool undef = (len == UNDEFINED);
    if (vr == "SQ" || (undef && vr == "UN")) {
      if (undef && vr == "UN") {
        // UN with undefined length is a sequence encoded implicit VR LE
        bool sl = c.little, si = c.implicit;
        c.little = true; c.implicit = true;
        value = parse_sequence(c, pos, (size_t) -1, local_cs);
        c.little = sl; c.implicit = si;
        vr = "SQ";
      } else {
        size_t seq_end = undef ? (size_t) -1 : pos + len;
        if (!undef && seq_end > c.n) seq_end = c.n;
        value = parse_sequence(c, pos, seq_end, local_cs);
        if (!undef) pos = seq_end;
      }
    } else if (undef) {
      // Encapsulated pixel data (or other undefined-length binary data)
      value = parse_encapsulated(c, pos);
    } else if (vr == "UN" && dict_vr(tag) != "UN") {
      // Like pydicom, replace UN with the dictionary VR when it is known
      std::string dvr = dict_vr(tag);
      if (dvr.find(" or ") != std::string::npos) {
        if (dvr == "US or SS") dvr = (c.pixel_rep == 1) ? "SS" : "US";
        else if (dvr == "OB or OW") dvr = "OW";
        else dvr = dvr.substr(0, 2);
      }
      if (pos + len > c.n) len = (uint32_t) (c.n - pos);
      if (dvr == "SQ") {
        bool sl = c.little, si = c.implicit;
        c.little = true; c.implicit = true;
        size_t p2 = pos;
        value = parse_sequence(c, p2, pos + len, local_cs);
        c.little = sl; c.implicit = si;
      } else {
        value = decode_value(c, dvr, pos, len, local_cs);
      }
      vr = dvr;
      pos += len;
    } else {
      if (pos + len > c.n) {
        c.warnings.push_back("Value of tag " + hex8(tag) + " runs past end of file; truncated");
        len = (uint32_t) (c.n - pos);
      }
      value = decode_value(c, vr, pos, len, local_cs);
      pos += len;
    }

    if (tag == 0x00080005) local_cs = charset_code(value);
    if (tag == 0x00280103 && TYPEOF(value) == INTSXP && Rf_length(value) > 0)
      c.pixel_rep = INTEGER(value)[0];

    names.push_back(hex8(tag));
    elems.push_back(make_element(tag, vr, value, undef));
  }

  List ds(elems.size());
  for (size_t i = 0; i < elems.size(); ++i) ds[i] = elems[i];
  ds.attr("names") = wrap(names);
  ds.attr("class") = "dicom_dataset";
  return ds;
}

bool looks_explicit(const uint8_t* b, size_t n, size_t p) {
  if (p + 6 > n) return false;
  std::string vr((const char*) b + p + 4, 2);
  return is_known_vr(vr);
}

std::vector<uint8_t> inflate_raw(const uint8_t* data, size_t n) {
  std::vector<uint8_t> out;
  z_stream zs;
  std::memset(&zs, 0, sizeof(zs));
  if (inflateInit2(&zs, -MAX_WBITS) != Z_OK) Rcpp::stop("zlib init failed");
  zs.next_in = const_cast<Bytef*>(data);
  zs.avail_in = (uInt) n;
  std::vector<uint8_t> chunk(1 << 16);
  int ret;
  do {
    zs.next_out = chunk.data();
    zs.avail_out = (uInt) chunk.size();
    ret = inflate(&zs, Z_NO_FLUSH);
    if (ret != Z_OK && ret != Z_STREAM_END) {
      inflateEnd(&zs);
      Rcpp::stop("Failed to inflate deflated transfer syntax data");
    }
    out.insert(out.end(), chunk.data(), chunk.data() + (chunk.size() - zs.avail_out));
  } while (ret != Z_STREAM_END && zs.avail_in > 0);
  inflateEnd(&zs);
  return out;
}

}  // namespace

// [[Rcpp::export(.parse_dicom_raw)]]
List parse_dicom_raw(RawVector data, bool stop_before_pixels = false, bool force = false) {
  const uint8_t* b = data.begin();
  size_t n = data.size();
  Ctx c{b, n, true, false, false, stop_before_pixels, -1, {}};

  SEXP preamble = R_NilValue;
  size_t pos = 0;
  bool has_prefix = n >= 132 && std::memcmp(b + 128, "DICM", 4) == 0;
  if (has_prefix) {
    RawVector pre(128);
    std::memcpy(&pre[0], b, 128);
    preamble = pre;
    pos = 132;
  } else if (n >= 4 && std::memcmp(b, "DICM", 4) == 0) {
    pos = 4;  // prefix without preamble
  } else if (!force) {
    // Accept files without preamble if they start with group 0002 or 0008
    // in a recognisable way; otherwise ask for force = TRUE.
    if (n < 8) Rcpp::stop("File is too small to be DICOM");
    uint16_t g = (uint16_t) (b[0] | (b[1] << 8));
    uint16_t gb = (uint16_t) ((b[0] << 8) | b[1]);
    if (!(g == 0x0002 || g == 0x0008 || gb == 0x0008)) {
      Rcpp::stop("File is missing the 'DICM' prefix. Use force = TRUE to read it anyway.");
    }
  }

  // File Meta Information: always explicit VR little endian.
  List meta;
  std::string ts;
  bool has_meta = false;
  if (pos + 6 <= n && b[pos] == 0x02 && b[pos + 1] == 0x00) {
    c.little = true;
    c.implicit = !looks_explicit(b, n, pos);  // some broken files use implicit
    meta = parse_dataset(c, pos, n, CS_DEFAULT, true);
    has_meta = true;
    if (meta.containsElementNamed("00020010")) {
      List tse = meta["00020010"];
      SEXP v = tse["value"];
      if (TYPEOF(v) == STRSXP && Rf_length(v) > 0) ts = CHAR(STRING_ELT(v, 0));
    }
  } else {
    meta = List(0);
    meta.attr("names") = CharacterVector(0);
    meta.attr("class") = "dicom_dataset";
  }

  // Work out the encoding of the main data set.
  bool deflated = false;
  bool guessed = false;
  if (ts == "1.2.840.10008.1.2") { c.little = true; c.implicit = true; }
  else if (ts == "1.2.840.10008.1.2.2") { c.little = false; c.implicit = false; }
  else if (ts == "1.2.840.10008.1.2.1") { c.little = true; c.implicit = false; }
  else if (ts == "1.2.840.10008.1.2.1.99") { c.little = true; c.implicit = false; deflated = true; }
  else if (!ts.empty()) { c.little = true; c.implicit = false; c.encapsulated = true; }
  else {
    guessed = true;
    if (pos + 6 > n) Rcpp::stop("No data set found after the file meta information");
    uint16_t gl = (uint16_t) (b[pos] | (b[pos + 1] << 8));
    uint16_t gb = (uint16_t) ((b[pos] << 8) | b[pos + 1]);
    c.little = !(gl > 0x00FF && gb <= 0x00FF);
    c.implicit = !looks_explicit(b, n, pos);
  }

  // Sanity check: some files declare one VR encoding but use the other.
  if (!deflated && pos + 6 <= n) {
    bool expl = looks_explicit(b, n, pos);
    if (c.implicit && expl) { c.implicit = false; c.warnings.push_back("Expected implicit VR, but found explicit VR - using explicit VR for reading"); }
    else if (!c.implicit && !expl) { c.implicit = true; c.warnings.push_back("Expected explicit VR, but found implicit VR - using implicit VR for reading"); }
  }

  List ds;
  if (deflated) {
    std::vector<uint8_t> inflated = inflate_raw(b + pos, n - pos);
    Ctx d{inflated.data(), inflated.size(), true, false, false, stop_before_pixels, -1, {}};
    size_t p2 = 0;
    ds = parse_dataset(d, p2, inflated.size(), CS_DEFAULT, false);
    c.warnings.insert(c.warnings.end(), d.warnings.begin(), d.warnings.end());
  } else {
    ds = parse_dataset(c, pos, n, CS_DEFAULT, false);
  }

  return List::create(_["preamble"] = preamble,
                      _["meta"] = meta,
                      _["has_meta"] = has_meta,
                      _["dataset"] = ds,
                      _["transfer_syntax"] = ts,
                      _["little_endian"] = c.little,
                      _["implicit_vr"] = c.implicit,
                      _["encapsulated"] = c.encapsulated,
                      _["deflated"] = deflated,
                      _["guessed_encoding"] = guessed,
                      _["warnings"] = wrap(c.warnings));
}

// [[Rcpp::export(.deflate_raw)]]
RawVector deflate_raw(RawVector data) {
  z_stream zs;
  std::memset(&zs, 0, sizeof(zs));
  if (deflateInit2(&zs, Z_DEFAULT_COMPRESSION, Z_DEFLATED, -MAX_WBITS, 8, Z_DEFAULT_STRATEGY) != Z_OK)
    Rcpp::stop("zlib init failed");
  zs.next_in = data.begin();
  zs.avail_in = (uInt) data.size();
  std::vector<uint8_t> out;
  std::vector<uint8_t> chunk(1 << 16);
  int ret;
  do {
    zs.next_out = chunk.data();
    zs.avail_out = (uInt) chunk.size();
    ret = deflate(&zs, Z_FINISH);
    out.insert(out.end(), chunk.data(), chunk.data() + (chunk.size() - zs.avail_out));
  } while (ret == Z_OK);
  deflateEnd(&zs);
  RawVector r(out.size());
  if (!out.empty()) std::memcpy(&r[0], out.data(), out.size());
  return r;
}
