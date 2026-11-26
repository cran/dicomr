// Data dictionary lookups shared by the parser and the R API.
#include <Rcpp.h>
#include <unordered_map>
#include <string>
#include <cstring>
#include "dict_data.h"
#include "dicomr.h"

using namespace Rcpp;

static std::unordered_map<uint32_t, int>& dict_index() {
  static std::unordered_map<uint32_t, int> idx;
  if (idx.empty()) {
    idx.reserve(DICOM_DICT_N * 2);
    for (int i = 0; i < DICOM_DICT_N; ++i) idx[DICOM_DICT[i].tag] = i;
  }
  return idx;
}

static std::unordered_map<std::string, int>& keyword_index() {
  static std::unordered_map<std::string, int> idx;
  if (idx.empty()) {
    idx.reserve(DICOM_DICT_N * 2);
    for (int i = 0; i < DICOM_DICT_N; ++i) {
      if (DICOM_DICT[i].keyword[0] != '\0') idx[DICOM_DICT[i].keyword] = i;
    }
  }
  return idx;
}

static int repeater_match(uint32_t tag) {
  char hex[9];
  std::snprintf(hex, sizeof(hex), "%08X", tag);
  for (int i = 0; i < DICOM_REPEATERS_N; ++i) {
    const char* m = DICOM_REPEATERS[i].mask;
    bool ok = true;
    for (int k = 0; k < 8; ++k) {
      if (m[k] != 'X' && m[k] != hex[k]) { ok = false; break; }
    }
    if (ok) return i;
  }
  return -1;
}

// Returns the dictionary VR string for a tag ("UN" if unknown).
std::string dict_vr(uint32_t tag) {
  uint16_t group = tag >> 16, elem = tag & 0xFFFF;
  if (elem == 0x0000) return "UL";                       // group length
  if (group & 1) {                                        // private group
    if (elem >= 0x0010 && elem <= 0x00FF) return "LO";    // private creator
    return "UN";
  }
  auto& idx = dict_index();
  auto it = idx.find(tag);
  if (it != idx.end()) return DICOM_DICT[it->second].vr;
  int r = repeater_match(tag);
  if (r >= 0) return DICOM_REPEATERS[r].vr;
  return "UN";
}

// [[Rcpp::export(.dict_lookup_tags)]]
DataFrame dict_lookup_tags(NumericVector tags) {
  R_xlen_t n = tags.size();
  CharacterVector vr(n), vm(n), name(n), keyword(n);
  LogicalVector retired(n);
  auto& idx = dict_index();
  for (R_xlen_t i = 0; i < n; ++i) {
    if (NumericVector::is_na(tags[i])) {
      vr[i] = NA_STRING; vm[i] = NA_STRING; name[i] = NA_STRING; keyword[i] = NA_STRING;
      retired[i] = NA_LOGICAL; continue;
    }
    uint32_t tag = (uint32_t) tags[i];
    uint16_t group = tag >> 16, elem = tag & 0xFFFF;
    auto it = idx.find(tag);
    if (it != idx.end()) {
      const DictEntry& e = DICOM_DICT[it->second];
      vr[i] = e.vr; vm[i] = e.vm; name[i] = e.name; keyword[i] = e.keyword;
      retired[i] = e.retired[0] != '\0';
      continue;
    }
    if (elem == 0 ) {
      vr[i] = "UL"; vm[i] = "1"; name[i] = "Group Length"; keyword[i] = ""; retired[i] = false; continue;
    }
    if (group & 1) {
      if (elem >= 0x10 && elem <= 0xFF) {
        vr[i] = "LO"; vm[i] = "1"; name[i] = "Private Creator"; keyword[i] = "";
      } else {
        vr[i] = "UN"; vm[i] = "1"; name[i] = "Private tag data"; keyword[i] = "";
      }
      retired[i] = false; continue;
    }
    int r = repeater_match(tag);
    if (r >= 0) {
      const RepeaterEntry& e = DICOM_REPEATERS[r];
      vr[i] = e.vr; vm[i] = e.vm; name[i] = e.name; keyword[i] = e.keyword;
      retired[i] = e.retired[0] != '\0';
      continue;
    }
    vr[i] = "UN"; vm[i] = NA_STRING; name[i] = "Unknown"; keyword[i] = ""; retired[i] = false;
  }
  return DataFrame::create(_["vr"] = vr, _["vm"] = vm, _["name"] = name,
                           _["keyword"] = keyword, _["retired"] = retired,
                           _["stringsAsFactors"] = false);
}

// [[Rcpp::export(.dict_keyword_to_tag)]]
NumericVector dict_keyword_to_tag(CharacterVector keywords) {
  auto& idx = keyword_index();
  NumericVector out(keywords.size());
  for (R_xlen_t i = 0; i < keywords.size(); ++i) {
    if (CharacterVector::is_na(keywords[i])) { out[i] = NA_REAL; continue; }
    auto it = idx.find(std::string(keywords[i]));
    out[i] = (it == idx.end()) ? NA_REAL : (double) DICOM_DICT[it->second].tag;
  }
  return out;
}

// [[Rcpp::export(.dict_table)]]
DataFrame dict_table() {
  int n = DICOM_DICT_N;
  NumericVector tag(n);
  CharacterVector vr(n), vm(n), name(n), keyword(n);
  LogicalVector retired(n);
  for (int i = 0; i < n; ++i) {
    const DictEntry& e = DICOM_DICT[i];
    tag[i] = e.tag; vr[i] = e.vr; vm[i] = e.vm; name[i] = e.name;
    keyword[i] = e.keyword; retired[i] = e.retired[0] != '\0';
  }
  return DataFrame::create(_["tag"] = tag, _["vr"] = vr, _["vm"] = vm, _["name"] = name,
                           _["keyword"] = keyword, _["retired"] = retired,
                           _["stringsAsFactors"] = false);
}
