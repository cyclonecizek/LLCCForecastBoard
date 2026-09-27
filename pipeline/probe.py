"""Check each source: where it resolved, latest cycle, and which fields exist.

    python -m pipeline.probe
"""
from __future__ import annotations

import os
import time

import yaml

from .common import SESSION, floor_hour, iso
from .grib import AREA, Missing, match_fields, read_idx
from .run import ROOT
from .src_ncep import _fmt
from .src_web import ENS_URL, VARS


def probe_grib(bases, tmpl, cycles, show):
    for base in bases:
        for c in cycles:
            url = f"{base}/{_fmt(tmpl, c, 1 if 'blend' not in tmpl else 6)}"
            try:
                inv = read_idx(url + ".idx")
            except Missing:
                continue
            except Exception as e:
                print(f"    {base}: {e}")
                break
            print(f"    OK {iso(c)}  {url}")
            show(inv)
            return
        print(f"    -- nothing in {base}")


def show_area(inv):
    found = match_fields(inv, AREA)
    for k in AREA:
        print(f"       {k:8s} {found[k][3] if k in found else '-- not in inventory (rules using it are skipped for this model)'}")


def show_tstm(inv):
    lines = [r[3] for r in inv if ":TSTM:" in r[3]]
    print("\n".join(f"       {x}" for x in lines) or "       -- no TSTM records")


def main():
    with open(os.path.join(ROOT, "config.yaml")) as f:
        cfg = yaml.safe_load(f)
    now = int(time.time())
    for s in cfg["sources"]:
        if not s.get("enabled", True):
            continue
        print(f"\n[{s['id']}] {s.get('label', '')} ({s['kind']}, {s.get('tier', 'direct')})")
        if s["kind"] == "tle":
            probe_grib([s["base"]], s["file"], [floor_hour(now) - k * 3600 for k in range(10)], show_area)
        elif s["kind"] == "multi_model":
            for comp in s["components"]:
                print(f"  {comp['id']}:")
                cyc = [c for c in (floor_hour(now) - k * 3600 for k in range(36)) if time.gmtime(c).tm_hour in set(comp.get("cycles", [0, 6, 12, 18]))]
                probe_grib([comp["base"]], comp["file"], cyc, show_area)
        elif s["kind"] == "nbm_prob":
            cyc = [c for c in (floor_hour(now) - k * 3600 for k in range(36)) if time.gmtime(c).tm_hour in set(s.get("cycles", [0, 6, 12, 18]))]
            probe_grib([s["base"]], s["file"], cyc, show_tstm)
        elif s["kind"] == "openmeteo_ens":
            p = cfg["pads"][0]
            r = SESSION.get(ENS_URL, timeout=60, params={"latitude": p["lat"], "longitude": p["lon"], "models": s["model"],
                            "hourly": ",".join(VARS), "forecast_days": 2, "timeformat": "unixtime"})
            if r.status_code != 200:
                print(f"    HTTP {r.status_code}: {r.text[:200]}")
                continue
            h = r.json()["hourly"]
            for v in VARS:
                n = sum(x is not None for x in h.get(v, []))
                print(f"    {v:18s} {'-- missing' if v not in h else f'{n} non-null hours'}")


if __name__ == "__main__":
    main()
