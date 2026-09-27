"""Build docs/data/board.json for the LLCC probability board.

    python -m pipeline.run [--only hrrr,gefs]
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time

import yaml

from . import src_ncep, src_web
from .common import Context, PointCache, SourceResult, floor_hour, iso, log

KINDS = {
    "tle": src_ncep.tle,
    "multi_model": src_ncep.multi_model,
    "nbm_prob": src_ncep.nbm_prob,
    "openmeteo_ens": src_web.openmeteo_ens,
}
RULES = ("lt", "cu", "an", "db", "tk", "dw")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def add_debris(series: dict) -> dict:
    """Debris (4.1.6, rough): mid/high cloud now, no deep convection now, and deep
    convection within 10 nmi at some point in the previous 3 hours."""
    times = sorted(series)
    suffixes = {k[3:] for rec in series.values() for k in rec if k.startswith("d10")}
    for suf in suffixes:
        d, cl, db = f"d10{suf}", f"cl{suf}", f"db{suf}"
        for t in times:
            rec = series[t]
            if cl not in rec or d not in rec:
                continue
            past = [series[t - h * 3600].get(d) for h in (1, 2, 3) if (t - h * 3600) in series]
            past = [x for x in past if x is not None]
            if not past:
                continue
            rec[db] = int(bool(rec[cl]) and not rec[d] and any(past))
    for rec in series.values():
        for k in [k for k in rec if k.startswith(("d10", "cl"))]:
            del rec[k]
    return series


def build(cfg: dict, only: set | None = None) -> dict:
    now = int(time.time())
    days = int(cfg["days"])
    t_start = floor_hour(now) - 12 * 3600
    t_end = floor_hour(now) + (days + 1) * 86400
    timeline = list(range(t_start, t_end + 1, 3600))
    idx = {t: i for i, t in enumerate(timeline)}
    pads = cfg["pads"]
    cache = PointCache(os.path.join(ROOT, "cache", "points.json"))
    cache.prune(now - 4 * 86400)
    lat = sum(p["lat"] for p in pads) / len(pads)
    lon = sum(p["lon"] for p in pads) / len(pads)
    ctx = Context(now=now, lat=lat, lon=lon, t_start=t_start, t_end=t_end, cache=cache,
                  pads=pads, rules=cfg["rules"], root=ROOT,
                  om_refresh_h=float(cfg.get("openmeteo_refresh_hours", 3)))

    sources = []
    for scfg in cfg["sources"]:
        if (only and scfg["id"] not in only) or not scfg.get("enabled", True):
            continue
        t0 = time.time()
        try:
            res = KINDS[scfg["kind"]](scfg, ctx)
        except Exception as e:
            log.exception("%s failed", scfg["id"])
            res = SourceResult({}, status="error", note=f"{type(e).__name__}: {e}"[:200])
        members = []
        for mid, series in sorted(res.members.items()):
            series = add_debris({int(t): dict(r) for t, r in series.items()})
            v = {}
            for t, rec in series.items():
                i = idx.get(int(t))
                if i is None:
                    continue
                for k, x in rec.items():
                    v.setdefault(k, [None] * len(timeline))[i] = round(float(x), 2) if isinstance(x, float) else int(x)
            if v:
                members.append({"id": mid, "v": v})
        el = round(time.time() - t0, 1)
        log.info("%-10s %-8s %3d members %6.1fs  %s", scfg["id"], res.status, len(members), el, res.note)
        sources.append({"id": scfg["id"], "label": scfg.get("label", scfg["id"]),
                        "family": scfg.get("family", "global"), "tier": scfg.get("tier", "direct"),
                        "weight": float(scfg.get("weight", 1)), "status": res.status, "cycle": res.cycle,
                        "note": res.note, "seconds": el, "members": members})
    cache.save()
    return {"generated": iso(now), "generated_unix": now, "display_tz": cfg["display_tz"],
            "days": days, "window_hours": int(cfg.get("window_hours", 3)),
            "min_window_coverage": cfg.get("min_window_coverage", 0.5),
            "pads": pads, "rules": cfg["rules"], "times": timeline, "sources": sources}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=os.path.join(ROOT, "config.yaml"))
    ap.add_argument("--only", default="")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    with open(a.config) as f:
        cfg = yaml.safe_load(f)
    data = build(cfg, set(filter(None, a.only.split(","))) or None)
    if not any(s["members"] for s in data["sources"]):
        log.error("no members from any source; keeping previous output")
        return 1
    out = os.path.join(ROOT, "docs", "data")
    os.makedirs(out, exist_ok=True)
    tmp = os.path.join(out, "board.json.tmp")
    with open(tmp, "w") as f:
        json.dump(data, f, separators=(",", ":"))
    os.replace(tmp, os.path.join(out, "board.json"))
    log.info("wrote board.json (%.0f kB)", os.path.getsize(os.path.join(out, "board.json")) / 1024)
    return 0


if __name__ == "__main__":
    sys.exit(main())
