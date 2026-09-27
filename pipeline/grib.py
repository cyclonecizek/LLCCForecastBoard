"""Read what we need from remote GRIB2 files without downloading them whole.

- Records are located from the .idx inventory and fetched with HTTP Range.
  NOMADS accepts several ranges in one request (one hit per file); AWS does
  not, so there adjacent records are merged and fetched range by range.
- Every field is reduced to its maximum over the grid points within each radius
  of each pad's flight path (masks are computed once per grid and reused).
"""
from __future__ import annotations

import math
import os
import re
import tempfile
import threading
import time

import numpy as np

from .common import SESSION, log

try:
    import eccodes
    eccodes.codes_grib_multi_support_on()
except ImportError:
    eccodes = None

# Area fields: the maximum within each radius of each pad's flight path is kept.
AREA = {
    "ltng": r":LTNG:entire atmosphere",
    "refd263": r":REFD:263 K level:",
    "retop": r":RETOP:",
    "hcdc": r":HCDC:high cloud layer:",
    "mcdc": r":MCDC:middle cloud layer:",
    "refc": r":REFC:entire atmosphere",
}
_SKIP = re.compile(r"(\bmax\b|\bmin\b|\bave\b|\bacc\b|prob|%)", re.I)
NM = 1852.0


class Missing(Exception):
    """File or inventory not published (yet)."""


class _RateLimit:
    """Evenly spaced requests; NOMADS blocks IPs above ~120 hits/minute."""

    def __init__(self, per_minute: float):
        self.gap = 60.0 / per_minute
        self.next = 0.0
        self.lock = threading.Lock()

    def wait(self):
        with self.lock:
            now = time.monotonic()
            t = max(now, self.next)
            self.next = t + self.gap
        if t > now:
            time.sleep(t - now)


NOMADS_LIMIT = _RateLimit(90)


def _nomads(url: str) -> bool:
    return "nomads.ncep.noaa.gov" in url


def _get(url: str, **kw):
    if _nomads(url):
        NOMADS_LIMIT.wait()
    return SESSION.get(url, **kw)


def exists(url: str) -> bool:
    try:
        if _nomads(url):
            NOMADS_LIMIT.wait()
        return SESSION.head(url, timeout=15, allow_redirects=True).status_code == 200
    except Exception:
        return False


def read_idx(idx_url: str):
    r = _get(idx_url, timeout=30)
    if r.status_code in (403, 404):
        raise Missing(idx_url)
    r.raise_for_status()
    recs = []
    for line in r.text.splitlines():
        parts = line.split(":")
        if len(parts) < 3:
            continue
        try:
            start = int(parts[1])
        except ValueError:
            continue
        recs.append((parts[0], start, ":" + ":".join(parts[2:])))
    out = []
    for i, (no, start, desc) in enumerate(recs):
        end = next((recs[j][1] - 1 for j in range(i + 1, len(recs)) if recs[j][1] > start), None)
        out.append((no, start, end, desc))
    return out


def match_fields(inv, wanted: dict[str, str]) -> dict[str, tuple]:
    found = {}
    for name, pat in wanted.items():
        rx = re.compile(pat)
        for rec in inv:
            m = rx.search(rec[3])
            if m and not _SKIP.search(rec[3][m.end():]):
                found[name] = rec
                break
    return found


# ---------------------------------------------------------------- decoding
_MASKS: dict = {}
_MASK_LOCK = threading.Lock()


def _grid_key(gid):
    return tuple(eccodes.codes_get(gid, k) for k in (
        "gridType", "Ni", "Nj", "latitudeOfFirstGridPointInDegrees", "longitudeOfFirstGridPointInDegrees"))


def _masks(gid, pads, radii_nm):
    """{pad_id: {radius: grid indices within radius of the pad's flight-path segment}}."""
    key = (_grid_key(gid), tuple((p["id"], p["lat"], p["lon"], p["azimuth_deg"], p["length_nmi"]) for p in pads),
           tuple(radii_nm))
    with _MASK_LOCK:
        if key in _MASKS:
            return _MASKS[key]
    lats = np.asarray(eccodes.codes_get_array(gid, "latitudes"))
    lons = np.asarray(eccodes.codes_get_array(gid, "longitudes"))
    out = {}
    for p in pads:
        # local plane in nmi around the pad; segment from pad along the azimuth
        x = (((lons - p["lon"] + 180) % 360) - 180) * 60.0 * math.cos(math.radians(p["lat"]))
        y = (lats - p["lat"]) * 60.0
        az = math.radians(p["azimuth_deg"])
        sx, sy = p["length_nmi"] * math.sin(az), p["length_nmi"] * math.cos(az)
        t = np.clip((x * sx + y * sy) / (sx * sx + sy * sy), 0.0, 1.0)
        d = np.hypot(x - t * sx, y - t * sy)
        near = d <= max(radii_nm)
        out[p["id"]] = {r: np.nonzero(near & (d <= r))[0] for r in radii_nm}
    with _MASK_LOCK:
        _MASKS[key] = out
    return out


def _nearest(gid, lat, lon):
    for lo in (lon, lon % 360.0):
        try:
            r = eccodes.codes_grib_find_nearest(gid, lat, lo)[0]
            return float(r["value"] if isinstance(r, dict) else getattr(r, "value"))
        except Exception:
            continue
    raise RuntimeError("nearest-point lookup failed")


def _meta(gid) -> dict:
    meta = {"grid": eccodes.codes_get(gid, "gridType")}
    try:
        meta["rel"] = int(eccodes.codes_get(gid, "uvRelativeToGrid"))
    except Exception:
        meta["rel"] = 0
    if meta["grid"] == "lambert":
        meta["lov"] = float(eccodes.codes_get(gid, "LoVInDegrees"))
        meta["latin1"] = float(eccodes.codes_get(gid, "Latin1InDegrees"))
    return meta


def _decode(blob: bytes, jobs: list, pads, radii) -> list:
    """jobs[i] = 'point' | 'area' for the i-th message in blob."""
    out = []
    fd, path = tempfile.mkstemp(suffix=".grib2")
    try:
        with os.fdopen(fd, "wb") as f:
            f.write(blob)
        with open(path, "rb") as f:
            i = 0
            while True:
                gid = eccodes.codes_grib_new_from_file(f)
                if gid is None:
                    break
                try:
                    job = jobs[i] if i < len(jobs) else None
                    miss = eccodes.codes_get(gid, "missingValue")
                    if job == "area":
                        vals = np.asarray(eccodes.codes_get_values(gid))
                        res = {}
                        for pid, byr in _masks(gid, pads, radii).items():
                            res[pid] = {}
                            for r, idx in byr.items():
                                sub = vals[idx]
                                sub = sub[(np.abs(sub - miss) > 1e-6) & (np.abs(sub) < 1e10)]
                                res[pid][r] = float(sub.max()) if sub.size else None
                        out.append((res, {}))
                    else:
                        out.append((None, {}))
                    i += 1
                finally:
                    eccodes.codes_release(gid)
    finally:
        os.unlink(path)
    return out


def _multipart(resp) -> list[tuple[int, bytes]]:
    """Parse a multipart/byteranges body into [(start, bytes)]."""
    ctype = resp.headers.get("Content-Type", "")
    if "multipart/byteranges" not in ctype:
        cr = resp.headers.get("Content-Range", "")
        m = re.search(r"bytes (\d+)-", cr)
        return [(int(m.group(1)) if m else 0, resp.content)]
    boundary = re.search(r"boundary=\"?([^\";]+)\"?", ctype).group(1).encode()
    parts = []
    for chunk in resp.content.split(b"--" + boundary):
        if b"\r\n\r\n" not in chunk:
            continue
        head, body = chunk.split(b"\r\n\r\n", 1)
        m = re.search(rb"Content-Range:\s*bytes (\d+)-(\d+)", head, re.I)
        if not m:
            continue
        a, b = int(m.group(1)), int(m.group(2))
        parts.append((a, body[: b - a + 1]))
    return parts


def fetch(grib_url: str, pads: list, radii_nm: list[float]) -> dict:
    """{field: {pad_id: {radius: max}}} for the AREA fields found. Raises Missing."""
    if eccodes is None:
        raise RuntimeError("eccodes not installed")
    inv = read_idx(grib_url + ".idx")
    found = sorted(match_fields(inv, AREA).items(), key=lambda kv: kv[1][1])
    vals, _ = _fetch_found(grib_url, found, {k: "area" for k in AREA}, pads, radii_nm)
    return vals


def fetch_select(grib_url: str, select, pads: list, radii_nm: list[float]) -> dict:
    """select(inventory) -> [(name, record)]; area maxima for each. Raises Missing."""
    if eccodes is None:
        raise RuntimeError("eccodes not installed")
    inv = read_idx(grib_url + ".idx")
    found = sorted(select(inv), key=lambda kv: kv[1][1])
    vals, _ = _fetch_found(grib_url, found, {n: "area" for n, _ in found}, pads, radii_nm)
    return vals


def _fetch_found(grib_url, found, kinds, pads, radii_nm):
    if not found:
        return {}, {}

    # contiguous runs of whole messages; a sub-message record stands alone
    groups: list[list] = []
    for name, rec in found:
        whole = "." not in rec[0]
        g = groups[-1] if groups else None
        if g and whole and g[-1][2] and g[-1][1][2] is not None and rec[1] == g[-1][1][2] + 1:
            g.append((name, rec, whole))
        else:
            groups.append([(name, rec, whole)])

    def rng(g):
        a, b = g[0][1][1], g[-1][1][2]
        return f"{a}-{b}" if b is not None else f"{a}-"

    blobs: dict[int, bytes] = {}
    if _nomads(grib_url) and len(groups) > 1 and all(g[-1][1][2] is not None for g in groups):
        r = _get(grib_url, headers={"Range": "bytes=" + ",".join(rng(g) for g in groups)}, timeout=90)
        if r.status_code in (403, 404, 416):
            raise Missing(grib_url)
        if r.status_code == 200:
            raise RuntimeError("server ignored byte ranges")
        r.raise_for_status()
        for a, body in _multipart(r):
            blobs[a] = body
    for g in groups:
        if g[0][1][1] in blobs:
            continue
        r = _get(grib_url, headers={"Range": "bytes=" + rng(g)}, timeout=90)
        if r.status_code in (403, 404, 416):
            raise Missing(grib_url)
        r.raise_for_status()
        blobs[g[0][1][1]] = r.content

    vals, meta = {}, {}
    for g in groups:
        blob = blobs.get(g[0][1][1])
        if blob is None:
            continue
        if len(g) == 1 and not g[0][2]:
            k = int(g[0][1][0].split(".")[1]) - 1
            jobs = [None] * k + [kinds[g[0][0]]]
            dec = _decode(blob, jobs, pads, radii_nm)
            pairs = [(g[0][0], dec[k] if k < len(dec) else (None, {}))]
        else:
            dec = _decode(blob, [kinds[n] for n, _, _ in g], pads, radii_nm)
            pairs = [(n, dec[i] if i < len(dec) else (None, {})) for i, (n, _, _) in enumerate(g)]
        for n, (v, m) in pairs:
            vals[n], meta[n] = v, m
    return vals, meta
