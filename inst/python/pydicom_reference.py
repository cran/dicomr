#!/usr/bin/env python3
"""Dump a DICOM file as seen by pydicom, for comparison with dicomr.

Usage: python3 pydicom_reference.py <input.dcm> <output_dir> [jpeg_plugin]

jpeg_plugin (default "pillow") picks pydicom's decoder for JPEG Baseline and
Extended data; use "default" to let pydicom choose.

Writes <output_dir>/tags.json (every element, with nested sequences) and,
when pydicom can decode the pixel data, <output_dir>/pixels.bin holding
pixel_array as float64 little endian in Fortran (column-major) order so R can
read it straight into an array of the same shape.
"""
import hashlib
import json
import os
import sys
import warnings

warnings.simplefilter("ignore")

import numpy as np
import pydicom
from pydicom.encaps import generate_fragments
from pydicom.multival import MultiValue
from pydicom.tag import BaseTag
from pydicom.valuerep import PersonName


def md5(b):
    return hashlib.md5(bytes(b)).hexdigest()


def scalar(x):
    if isinstance(x, PersonName):
        return str(x)
    if isinstance(x, BaseTag):
        return "%08X" % int(x)
    if isinstance(x, float):
        return float(x)
    if isinstance(x, int):
        return int(x)
    if x is None:
        return None
    return str(x)


def norm(elem, encapsulated):
    vr = elem.VR
    v = elem.value
    if vr == "SQ":
        return [dump(item, encapsulated) for item in v]
    if elem.tag == 0x7FE00010 and encapsulated:
        return {"fragments": [md5(f) for f in list(generate_fragments(v))[1:]]}  # skip offset table
    if v is None:
        return None
    if isinstance(v, (bytes, bytearray)):
        return {"md5": md5(v), "len": len(v)}
    if isinstance(v, (MultiValue, list, tuple)):
        out = [scalar(x) for x in v]
    else:
        out = [scalar(v)]
    if out == [""] or out == []:
        return None
    return out


def dump(ds, encapsulated):
    out = {}
    for elem in ds:
        out["%08X" % int(elem.tag)] = {"vr": elem.VR, "value": norm(elem, encapsulated)}
    return out


JPEG_8BIT = ("1.2.840.10008.1.2.4.50", "1.2.840.10008.1.2.4.51")


def main(path, outdir, jpeg_plugin="pillow"):
    os.makedirs(outdir, exist_ok=True)
    ds = pydicom.dcmread(path, force=True)
    ts = getattr(getattr(ds, "file_meta", None), "TransferSyntaxUID", None)
    encapsulated = bool(ts is not None and ts.is_encapsulated)
    result = {"meta": dump(ds.file_meta, False) if hasattr(ds, "file_meta") else {},
              "dataset": dump(ds, encapsulated), "pixel_shape": None, "pixel_error": None}
    if "PixelData" in ds:
        try:
            if ts in JPEG_8BIT and jpeg_plugin != "default":
                from pydicom.pixels import pixel_array
                arr = pixel_array(ds, decoding_plugin=jpeg_plugin)
                result["pixel_plugin"] = jpeg_plugin
            else:
                arr = ds.pixel_array
            result["pixel_shape"] = list(arr.shape)
            result["pixel_dtype"] = str(arr.dtype)
            np.asarray(arr, dtype="<f8").ravel(order="F").tofile(os.path.join(outdir, "pixels.bin"))
        except Exception as e:  # noqa: BLE001
            result["pixel_error"] = "%s: %s" % (type(e).__name__, e)
    with open(os.path.join(outdir, "tags.json"), "w") as f:
        json.dump(result, f)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], *(sys.argv[3:4]))
