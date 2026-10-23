#!/usr/bin/env python3
"""Genome-browser PDF (or an SVG made from one) -> a lean SVG whose track bars are merged, so it stays responsive in Inkscape.

A genome-browser track is drawn as thousands of tiny filled rectangles (one per bin). Opened in Inkscape each is its own object, and zoom insets are
often a whole track repeated and clipped to a small window. This script:

  1. converts a PDF to SVG (Inkscape's importer by default, or pdftocairo); an .svg input skips this step;
  2. drops paths that lie entirely outside the clip window of an ancestor (the repeated, clipped copies in zoom insets), and removes a clip when everything left
     already lies inside it;
  3. merges sibling paths with identical attributes (same fill, stroke, ...) into ONE compound path (all path data is converted to absolute coordinates first,
     since a relative "m" is not safe to concatenate). Within a stretch of adjacent plain paths this is done across different fills (the usual case: bars shaded in
     a few colors, interleaved, often layered light-under-dark) by ordering the fills so every overlap keeps its original stacking order (one merged path per
     fill, drawn in that order); if the stacking is contradictory, only consecutive runs are merged;
  4. removes clipPath definitions nothing uses any more, and groups emptied by the trimming.

Nothing else is changed: text, images, groups and transforms stay as they were. The fill rule is checked: paths using evenodd are never merged.

    python scripts/browser_pdf_to_svg.py tracks.pdf                    # -> tracks.svg
    python scripts/browser_pdf_to_svg.py tracks.pdf -o out.svg --converter pdftocairo
    python scripts/browser_pdf_to_svg.py overview.svg -o overview_merged.svg   # clean up an existing SVG (e.g. one saved from Inkscape)
    python scripts/browser_pdf_to_svg.py tracks.pdf --no-trim          # merge only; keep every clipped copy as is
    python scripts/browser_pdf_to_svg.py figure.svg -o figure_v2.svg --within layer2   # only inside one group (e.g. one panel's layer) of a composite

Check the result against the source after any change (render both, e.g. `inkscape f.svg --export-type=png`, and compare).
"""
import argparse
import math
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from lxml import etree

SVG = "http://www.w3.org/2000/svg"
NS = "{%s}" % SVG
ARITY = {"M": 2, "L": 2, "H": 1, "V": 1, "C": 6, "S": 4, "Q": 4, "T": 2, "A": 7, "Z": 0}
TOK = re.compile(r"([MmLlHhVvCcSsQqTtAaZz])|([-+]?(?:\d*\.\d+|\d+\.?)(?:[eE][-+]?\d+)?)")
EPS = 1e-6
OVERLAP_TOL = 0.02   # overlaps shallower than this (user units) count as touching: PDF/SVG coordinates are rounded to ~3 decimals, so neighboring bars overlap by ~0.001
GRAPHIC_OTHER = {"text", "image", "use", "foreignObject", "switch", "rect", "circle", "ellipse", "line", "polyline", "polygon"}


# ---- path data -----------------------------------------------------------------------------------------------------------------------------------
def to_absolute(d):
    """Parse path data into absolute segments [(CMD, [numbers])], or None if it uses something not handled here."""
    items = [m.group(1) or float(m.group(2)) for m in TOK.finditer(d or "")]
    out, i, n, cmd = [], 0, len(items), None
    cx = cy = sx = sy = 0.0
    while i < n:
        it = items[i]
        if isinstance(it, str):
            cmd = it
            i += 1
            if cmd in "Zz":
                out.append(("Z", []))
                cx, cy = sx, sy
                continue
        elif cmd is None:
            return None
        up, rel = cmd.upper(), cmd.islower()
        if up == "Z":
            return None                                   # numbers after a closepath
        k = ARITY[up]
        if i + k > n:
            return None
        a = items[i:i + k]
        if any(isinstance(x, str) for x in a):
            return None
        i += k
        if up == "H":
            x = a[0] + (cx if rel else 0.0)
            out.append(("H", [x])); cx = x
        elif up == "V":
            y = a[0] + (cy if rel else 0.0)
            out.append(("V", [y])); cy = y
        elif up == "A":
            x, y = a[5] + (cx if rel else 0.0), a[6] + (cy if rel else 0.0)
            out.append(("A", a[:5] + [x, y])); cx, cy = x, y
        else:                                              # M L C S Q T: all pairs shift by the segment's start point
            o = (cx, cy) if rel else (0.0, 0.0)
            pts = []
            for j in range(0, k, 2):
                pts += [a[j] + o[0], a[j + 1] + o[1]]
            out.append((up, pts)); cx, cy = pts[-2], pts[-1]
            if up == "M":
                sx, sy = cx, cy
                cmd = "l" if rel else "L"                  # further coordinate pairs after a moveto are linetos
    return out


def fmt(x):
    s = "%.9g" % x
    return "0" if s in ("-0", "0") else s


def segs_to_d(segs):
    return "".join(c + ",".join(fmt(v) for v in a) for c, a in segs)


def bbox_of(segs):
    """Conservative bounding box (control points included); None when an arc makes it unreliable."""
    xs, ys, cx, cy = [], [], 0.0, 0.0
    for c, a in segs:
        if c == "A":
            return None
        if c == "H":
            cx = a[0]; xs.append(cx); ys.append(cy)
        elif c == "V":
            cy = a[0]; xs.append(cx); ys.append(cy)
        elif c != "Z":
            xs += a[0::2]; ys += a[1::2]; cx, cy = a[-2], a[-1]
    return (min(xs), min(ys), max(xs), max(ys)) if xs else None


# ---- transforms ----------------------------------------------------------------------------------------------------------------------------------
IDENT = (1.0, 0.0, 0.0, 1.0, 0.0, 0.0)


def mmul(m, n):
    a1, b1, c1, d1, e1, f1 = m
    a2, b2, c2, d2, e2, f2 = n
    return (a1 * a2 + c1 * b2, b1 * a2 + d1 * b2, a1 * c2 + c1 * d2, b1 * c2 + d1 * d2, a1 * e2 + c1 * f2 + e1, b1 * e2 + d1 * f2 + f1)


def parse_transform(s):
    m = IDENT
    for name, args in re.findall(r"(\w+)\s*\(([^)]*)\)", s or ""):
        v = [float(x) for x in re.split(r"[\s,]+", args.strip()) if x]
        if name == "matrix" and len(v) == 6:
            t = tuple(v)
        elif name == "translate":
            t = (1, 0, 0, 1, v[0], v[1] if len(v) > 1 else 0.0)
        elif name == "scale":
            t = (v[0], 0, 0, v[1] if len(v) > 1 else v[0], 0, 0)
        elif name == "rotate":
            r = math.radians(v[0])
            t = (math.cos(r), math.sin(r), -math.sin(r), math.cos(r), 0, 0)
            if len(v) == 3:
                t = mmul(mmul((1, 0, 0, 1, v[1], v[2]), t), (1, 0, 0, 1, -v[1], -v[2]))
        elif name == "skewX":
            t = (1, 0, math.tan(math.radians(v[0])), 1, 0, 0)
        elif name == "skewY":
            t = (1, math.tan(math.radians(v[0])), 0, 1, 0, 0)
        else:
            raise ValueError("unsupported transform: " + name)
        m = mmul(m, t)
    return m


def bbox_transformed(bb, m):
    pts = [(bb[0], bb[1]), (bb[2], bb[1]), (bb[0], bb[3]), (bb[2], bb[3])]
    xs = [m[0] * x + m[2] * y + m[4] for x, y in pts]
    ys = [m[1] * x + m[3] * y + m[5] for x, y in pts]
    return (min(xs), min(ys), max(xs), max(ys))


# ---- clip windows --------------------------------------------------------------------------------------------------------------------------------
def clip_ref(el):
    v = el.get("clip-path")
    if not v and el.get("style"):
        m = re.search(r"clip-path:\s*([^;]+)", el.get("style"))
        v = m.group(1) if m else None
    m = re.search(r"url\(\s*#([^)\s]+)\s*\)", v or "")
    return m.group(1) if m else None


def clip_bbox(cp):
    """Union bounding box of a clipPath's shapes in the owner's coordinate system, or None if it holds anything but rects/paths (or objectBoundingBox units)."""
    if cp is None or cp.get("clipPathUnits") == "objectBoundingBox":
        return None
    outer = parse_transform(cp.get("transform"))
    boxes = []
    for ch in cp:
        if not isinstance(ch.tag, str):
            continue
        tag = etree.QName(ch).localname
        if tag == "rect":
            x, y, w, h = (float(ch.get(k, 0)) for k in ("x", "y", "width", "height"))
            bb = (x, y, x + w, y + h)
        elif tag == "path":
            segs = to_absolute(ch.get("d"))
            bb = bbox_of(segs) if segs else None
        else:
            return None
        if bb is None:
            return None
        boxes.append(bbox_transformed(bb, mmul(outer, parse_transform(ch.get("transform")))))
    if not boxes:
        return None
    return (min(b[0] for b in boxes), min(b[1] for b in boxes), max(b[2] for b in boxes), max(b[3] for b in boxes))


def has_visible_stroke(el):
    st = el.get("style", "")
    m = re.search(r"stroke:\s*([^;]+)", st)
    v = m.group(1).strip() if m else el.get("stroke", "none")
    return v not in ("none", "")


# ---- trimming ------------------------------------------------------------------------------------------------------------------------------------
class Owner:
    def __init__(self, el, box):
        self.el, self.box = el, box
        self.n_boundary = 0          # kept paths that are not fully inside the window
        self.other = False           # text, images, ... below the owner: the clip must stay
        self.stroked = False


def trim_clipped_copies(root, doc=None):
    """Trim below `root`; clipPath ids are looked up in `doc` (the whole document; default: root itself)."""
    ids = {e.get("id"): e for e in (root if doc is None else doc).iter() if isinstance(e.tag, str) and e.get("id")}
    owners, removed, emptied = [], [], set()

    def visit(el, active):
        """active: list of [Owner, matrix from the owner's coordinate system to el's PARENT coordinate system]."""
        if not isinstance(el.tag, str):
            return
        tag = etree.QName(el).localname
        if tag in ("defs", "clipPath", "mask", "pattern", "marker", "symbol", "style", "metadata", "namedview"):
            return
        t = parse_transform(el.get("transform"))
        here = [[o, mmul(m, t)] for o, m in active]          # now relative to el's own (inner) coordinates
        ref = clip_ref(el)
        if ref:
            box = clip_bbox(ids.get(ref))
            if box is not None:
                o = Owner(el, box)
                owners.append(o)
                here = here + [[o, IDENT]]                    # the window lives in el's coordinate system (after its own transform)
            else:
                pass
        if tag == "path":
            segs = to_absolute(el.get("d"))
            bb = bbox_of(segs) if segs else None
            drop, inside_all = False, True
            for o, m in here:
                if o.el is el:
                    continue                                   # a path carrying its own clip is not trimmed by it
                if bb is None:
                    inside_all = False
                    o.n_boundary += 1
                    continue
                b = bbox_transformed(bb, m)
                if b[2] < o.box[0] - EPS or b[0] > o.box[2] + EPS or b[3] < o.box[1] - EPS or b[1] > o.box[3] + EPS:
                    drop = True
                    break
            if drop:
                removed.append(el)
                return
            for o, m in here:
                if bb is not None:                             # also for a path carrying its own clip: it may reach past its window
                    b = bbox_transformed(bb, m)
                    if not (b[0] >= o.box[0] - EPS and b[1] >= o.box[1] - EPS and b[2] <= o.box[2] + EPS and b[3] <= o.box[3] + EPS):
                        o.n_boundary += 1
                if has_visible_stroke(el):
                    o.stroked = True
            return
        if tag in GRAPHIC_OTHER:
            for o, m in here:
                o.other = True                                 # includes el's own clip: a rect/text/image clipped by itself must keep that clip
            return
        for ch in el:
            visit(ch, here)

    visit(root, [])
    for el in removed:
        par = el.getparent()
        par.remove(el)
        emptied.add(par)
    # groups emptied by the trimming go too (repeatedly, upwards); only groups that held something before
    while emptied:
        nxt = set()
        for g in emptied:
            if etree.QName(g).localname == "g" and len([c for c in g if isinstance(c.tag, str)]) == 0 and g.getparent() is not None:
                p = g.getparent()
                p.remove(g)
                nxt.add(p)
        emptied = nxt
    dropped_clips = 0
    for o in owners:
        if o.el.getparent() is None and o.el is not root:
            continue
        if o.n_boundary == 0 and not o.other and not o.stroked:
            if o.el.get("clip-path"):
                del o.el.attrib["clip-path"]
                dropped_clips += 1
            elif o.el.get("style"):
                o.el.set("style", re.sub(r"clip-path:\s*[^;]+;?", "", o.el.get("style")))
                dropped_clips += 1
    return len(removed), dropped_clips


# ---- merging -------------------------------------------------------------------------------------------------------------------------------------
SKIP_ATTRS = {"d", "id"}


def attr_key(el):
    return tuple(sorted((k, v) for k, v in el.attrib.items() if k not in SKIP_ATTRS and "sodipodi" not in k and "inkscape:" not in k and "{http://www.inkscape.org" not in k))


def mergeable(el):
    if not isinstance(el.tag, str) or etree.QName(el).localname != "path" or len(el):
        return None
    if "evenodd" in (el.get("style", "") + el.get("fill-rule", "")):
        return None
    if any("marker" in k or "marker" in v for k, v in el.attrib.items() if k != "d"):
        return None
    segs = to_absolute(el.get("d"))
    return segs or None


def merge_segments(root, min_run):
    """Merge same-attribute paths across a stretch of plain (no transform) sibling paths, keeping every overlap's stacking order (see key_order)."""
    removed = 0
    for parent in list(root.iter()):
        if not isinstance(parent.tag, str):
            continue
        kids = [k for k in parent if isinstance(k.tag, str)]
        i = 0
        while i < len(kids):
            if kids[i].get("transform") or mergeable(kids[i]) is None:
                i += 1
                continue
            j, items = i, []
            while j < len(kids) and not kids[j].get("transform"):
                segs = mergeable(kids[j])
                if segs is None:
                    break
                items.append((kids[j], segs, bbox_of(segs)))
                j += 1
            keys = [attr_key(e) for e, _, _ in items]
            distinct = {k: n for n, k in enumerate(dict.fromkeys(keys))}
            if len(items) >= min_run and len(distinct) > 1 and all(bb is not None for _, _, bb in items):
                order = key_order(items, [distinct[k] for k in keys], len(distinct))
                if order is not None:
                    by_key = {}
                    for (e, segs, _), k in zip(items, keys):
                        by_key.setdefault(distinct[k], []).append((e, segs))
                    anchor = parent.index(items[0][0])
                    for e, _, _ in items:
                        parent.remove(e)
                    for pos, kid in enumerate(order):
                        members = by_key[kid]
                        members[0][0].set("d", "".join(segs_to_d(sg) for _, sg in members))
                        parent.insert(anchor + pos, members[0][0])
                        removed += len(members) - 1
            i = j
    return removed


def key_order(items, key_ids, n_keys):
    """Order of the distinct fills (key ids) that reproduces the original stacking: for every overlapping pair of paths with different fills, the fill of the
    earlier (lower) path must come before the later one's. Returns the order, or None if those constraints contradict each other."""
    order = sorted(range(len(items)), key=lambda n: items[n][2][0])
    active, edges, tol = [], set(), OVERLAP_TOL
    for n in order:
        x0, y0, x1, y1 = items[n][2]
        active = [m for m in active if items[m][2][2] > x0 + tol]
        for m in active:
            if key_ids[m] != key_ids[n]:
                a = items[m][2]
                if y0 < a[3] - tol and y1 > a[1] + tol:
                    lo, hi = (m, n) if m < n else (n, m)
                    edges.add((key_ids[lo], key_ids[hi]))
        active.append(n)
    indeg = {k: 0 for k in range(n_keys)}
    nxt = {k: [] for k in range(n_keys)}
    for a, b in edges:
        nxt[a].append(b)
        indeg[b] += 1
    # Kahn's algorithm, preferring the fills that appear first in the document
    ready = sorted(k for k, v in indeg.items() if v == 0)
    out = []
    while ready:
        k = ready.pop(0)
        out.append(k)
        for b in nxt[k]:
            indeg[b] -= 1
            if indeg[b] == 0:
                ready.append(b)
        ready.sort()
    return out if len(out) == n_keys else None


def merge_runs(root, min_run):
    n_before = n_after = 0
    for parent in list(root.iter()):
        if not isinstance(parent.tag, str):
            continue
        kids = [k for k in parent if isinstance(k.tag, str)]
        i = 0
        while i < len(kids):
            segs = mergeable(kids[i])
            if segs is None:
                i += 1
                continue
            key, run, runsegs = attr_key(kids[i]), [kids[i]], [segs]
            j = i + 1
            while j < len(kids):
                s2 = mergeable(kids[j])
                if s2 is None or attr_key(kids[j]) != key:
                    break
                run.append(kids[j]); runsegs.append(s2); j += 1
            n_before += len(run)
            if len(run) >= min_run:
                run[0].set("d", "".join(segs_to_d(s) for s in runsegs))
                for el in run[1:]:
                    parent.remove(el)
                n_after += 1
            else:
                n_after += len(run)
            i = j
    return n_before, n_after


def drop_unused_clippaths(root):
    used = {clip_ref(e) for e in root.iter() if isinstance(e.tag, str)}
    n = 0
    for cp in [e for e in root.iter() if isinstance(e.tag, str) and etree.QName(e).localname == "clipPath"]:
        if cp.get("id") not in used:
            cp.getparent().remove(cp)
            n += 1
    return n


# ---- driver --------------------------------------------------------------------------------------------------------------------------------------
def count(root):
    tags = [etree.QName(e).localname for e in root.iter() if isinstance(e.tag, str)]
    return tags.count("path"), tags.count("clipPath")


def pdf_to_svg(pdf, out, converter, page):
    if converter == "pdftocairo":
        cmd = ["pdftocairo", "-svg", "-f", str(page), "-l", str(page), str(pdf), str(out)]
    else:
        cmd = ["inkscape", str(pdf), "--export-type=svg", "--export-plain-svg", "--export-filename=" + str(out)]
        if page != 1:
            cmd.append("--pdf-page=%d" % page)
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter, epilog="See the module docstring for details.")
    ap.add_argument("input", type=Path, help="genome-browser .pdf, or an .svg to clean up")
    ap.add_argument("-o", "--output", type=Path, help="output .svg (default: input with .svg, or _merged.svg when the input is an svg)")
    ap.add_argument("--converter", choices=["inkscape", "pdftocairo"], default="inkscape", help="PDF -> SVG converter (default: inkscape; keeps text as text)")
    ap.add_argument("--page", type=int, default=1, help="PDF page to convert (default 1)")
    ap.add_argument("--no-trim", action="store_true", help="do not drop paths outside their clip window")
    ap.add_argument("--within", metavar="ID", help="only trim and merge inside the element with this id (e.g. one panel's layer of a composite); the rest is left as is")
    ap.add_argument("--min-run", type=int, default=2, help="merge runs of at least this many identical consecutive paths (default 2)")
    a = ap.parse_args()
    src = a.input
    out = a.output or (src.with_name(src.stem + "_merged.svg") if src.suffix.lower() == ".svg" else src.with_suffix(".svg"))
    if out.resolve() == src.resolve():
        sys.exit("refusing to overwrite the input; pass -o")
    tmp = None
    if src.suffix.lower() == ".pdf":
        tmp = Path(tempfile.mkdtemp()) / "converted.svg"
        pdf_to_svg(src, tmp, a.converter, a.page)
        svg_path = tmp
    else:
        svg_path = src
    tree = etree.parse(str(svg_path), etree.XMLParser(huge_tree=True))
    root = tree.getroot()
    work = root
    if a.within:
        hits = [e for e in root.iter() if isinstance(e.tag, str) and e.get("id") == a.within]
        if not hits:
            sys.exit("no element with id " + a.within)
        work = hits[0]
    p0, c0 = count(root)
    w0, _ = count(work)
    n_removed = n_clip_dropped = 0
    if not a.no_trim:
        # clip windows and paths are compared inside the subtree, so the transforms of its ancestors (shared by both) do not matter
        n_removed, n_clip_dropped = trim_clipped_copies(work, root)
    n_seg = merge_segments(work, a.min_run)
    mb, ma = merge_runs(work, a.min_run)
    n_cp = drop_unused_clippaths(root)
    p1, c1 = count(root)
    if a.within:
        print(f"  within #{a.within}: paths {w0:,} -> {count(work)[0]:,}")
    tree.write(str(out), xml_declaration=True, encoding="UTF-8")
    print(f"{src.name} -> {out.name}")
    print(f"  paths {p0:,} -> {p1:,}  (outside a clip window: -{n_removed:,}; merged across interleaved fills: -{n_seg:,}; consecutive runs merged: {mb:,} -> {ma:,})")
    print(f"  clipPaths {c0:,} -> {c1:,}  (clips dropped as redundant: {n_clip_dropped:,}; unused definitions removed: {n_cp:,})")
    print(f"  file {src.stat().st_size / 1e6:.1f} MB -> {out.stat().st_size / 1e6:.1f} MB")
    if tmp is not None:
        shutil.rmtree(tmp.parent, ignore_errors=True)


if __name__ == "__main__":
    main()
