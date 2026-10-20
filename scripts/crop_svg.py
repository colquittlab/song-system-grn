#!/usr/bin/env python3
"""Crop an SVG to the bounding box of its drawing (viewBox, width and height are rewritten; the content is not touched).

Figures drawn on a generous canvas (R devices, ComplexHeatmap, ...) end up with large margins that Inkscape does not remove on export; this measures the drawing with
Inkscape (`--query-all`, the union of every object's visual bounding box except the root) and sets the document to that box plus a small pad. Make the source
background transparent (e.g. svglite(bg = "transparent")): a full-canvas background rectangle is itself an object and defeats the crop.

    python scripts/crop_svg.py in.svg out.svg [--pad 0.5]      # pad in the file's own user units (pt for svglite output)
    prints the cropped size in inches
"""
import argparse
import re
import subprocess
import sys

from lxml import etree

UNIT_PX = {"px": 1.0, "pt": 96 / 72, "in": 96.0, "mm": 96 / 25.4, "cm": 96 / 2.54, "": 1.0}


def to_px(v):
    m = re.fullmatch(r"\s*([-0-9.eE]+)\s*([a-z%]*)\s*", v)
    return float(m.group(1)) * UNIT_PX[m.group(2)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("input")
    ap.add_argument("output")
    ap.add_argument("--pad", type=float, default=0.5, help="margin around the drawing, in the file's user units (default 0.5)")
    a = ap.parse_args()
    tree = etree.parse(a.input, etree.XMLParser(huge_tree=True))
    root = tree.getroot()
    vb = [float(x) for x in re.split(r"[\s,]+", root.get("viewBox").strip())]
    scale = to_px(root.get("width")) / vb[2]                      # document px per user unit
    out = subprocess.run(["inkscape", a.input, "--query-all"], capture_output=True, text=True, check=True).stdout.strip().splitlines()
    boxes = []
    for line in out[1:]:                                          # the first line is the root element
        p = line.split(",")
        if len(p) == 5:
            x, y, w, h = (float(v) for v in p[1:])
            page_w, page_h = vb[2] * scale, vb[3] * scale
            if w > 0 and h > 0 and not (w >= 0.98 * page_w and h >= 0.98 * page_h):      # skip full-canvas rectangles (device backgrounds, clip frames)
                boxes.append((x, y, x + w, y + h))
    if not boxes:
        sys.exit("no drawing found")
    x0 = min(b[0] for b in boxes) / scale - a.pad + vb[0]
    y0 = min(b[1] for b in boxes) / scale - a.pad + vb[1]
    x1 = max(b[2] for b in boxes) / scale + a.pad + vb[0]
    y1 = max(b[3] for b in boxes) / scale + a.pad + vb[1]
    w, h = x1 - x0, y1 - y0
    unit = re.sub(r"[-0-9.eE\s]", "", root.get("width")) or "px"
    root.set("viewBox", "%.4f %.4f %.4f %.4f" % (x0, y0, w, h))
    root.set("width", "%.4f%s" % (w * (scale / UNIT_PX[unit]), unit))
    root.set("height", "%.4f%s" % (h * (scale / UNIT_PX[unit]), unit))
    tree.write(a.output, xml_declaration=True, encoding="UTF-8")
    print("%.3f x %.3f in" % (w * scale / 96, h * scale / 96))


if __name__ == "__main__":
    main()
