"""Global-ensemble stand-ins (Open-Meteo), used only where no direct evidence exists.

Each ensemble is sampled at the midpoint of the pads plus rings at 10 and 25 nmi.
"Deep convection" at a point = precip >= precip_in with CAPE >= cape_jkg. Records are
pad-agnostic (keys without "@pad"): at 25 km resolution the two pads are the same place.
    cu   deep convection within 10 nmi          (Cumulus stand-in)
    an   deep convection within 25 nmi + high cloud >= anvil_high_cloud_pct
    tk   mid cloud >= thick_mid_cloud_pct + precip within 10 nmi
    dw   mid cloud >= thick_mid_cloud_pct + moderate precip within 10 nmi
    d10, cl   helpers for the debris rule
No lightning: these models have no explicit lightning parameter.
"""
from __future__ import annotations

import json
import math
import os
import re
import time

from .common import SESSION, Context, SourceResult, log

ENS_URL = "https://ensemble-api.open-meteo.com/v1/ensemble"
VARS = ["precipitation", "cape", "cloud_cover_mid", "cloud_cover_high"]
_KEY = re.compile(r"^(" + "|".join(VARS) + r")(?:_member(\d+))?$")


def ring_points(lat, lon):
    pts = [(lat, lon, 0.0)]
    for dist, start in ((10, 0), (25, 30)):
        for k in range(6):
            b = math.radians(start + 60 * k)
            pts.append((round(lat + dist / 60.0 * math.cos(b), 4),
                        round(lon + dist / 60.0 * math.sin(b) / math.cos(math.radians(lat)), 4), float(dist)))
    return pts


def _center(ctx):
    return (sum(p["lat"] for p in ctx.pads) / len(ctx.pads), sum(p["lon"] for p in ctx.pads) / len(ctx.pads))


def _fetch(model, ctx, pts):
    variables = list(VARS)
    days = max(1, math.ceil((ctx.t_end - ctx.now) / 86400) + 1)
    for _ in range(len(VARS)):
        r = SESSION.get(ENS_URL, timeout=120, params={
            "latitude": ",".join(str(p[0]) for p in pts), "longitude": ",".join(str(p[1]) for p in pts),
            "models": model, "hourly": ",".join(variables), "precipitation_unit": "inch",
            "timeformat": "unixtime", "timezone": "GMT", "past_days": 1, "forecast_days": min(days, 16)})
        if r.status_code != 400:
            break
        bad = [v for v in variables if v in r.text]
        if not bad or len(variables) == 1:
            break
        log.info("%s: dropping unsupported %s", model, bad[0])
        variables.remove(bad[0])
    r.raise_for_status()
    js = r.json()
    return js if isinstance(js, list) else [js]


def _members(locs, pts, ctx):
    c = ctx.rules
    near10 = [i for i, p in enumerate(pts) if p[2] <= 10]
    times = locs[0]["hourly"]["time"]
    per: dict[str, dict] = {}
    for li, loc in enumerate(locs):
        for key, arr in loc["hourly"].items():
            m = _KEY.match(key)
            if m:
                mid = f"m{int(m.group(2)):02d}" if m.group(2) else "m00"
                per.setdefault(mid, {}).setdefault(m.group(1), [None] * len(locs))[li] = arr

    out = {}
    for mid, f in per.items():
        def at(var, li, i):
            a = f.get(var)
            return None if not a or a[li] is None else a[li][i]
        series = {}
        for i, t in enumerate(times):
            if not ctx.in_window(t):
                continue
            pr = [at("precipitation", li, i) for li in range(len(pts))]
            cp = [at("cape", li, i) for li in range(len(pts))]
            mid_c, high_c = at("cloud_cover_mid", 0, i), at("cloud_cover_high", 0, i)
            conv = [None if pr[li] is None or cp[li] is None else
                    (pr[li] >= c["precip_in"] and cp[li] >= c["cape_jkg"]) for li in range(len(pts))]
            rec = {}
            if any(conv[li] is not None for li in near10):
                rec["cu"] = rec["d10"] = int(any(conv[li] for li in near10 if conv[li] is not None))
            if high_c is not None and any(x is not None for x in conv):
                rec["an"] = int(any(x for x in conv if x is not None) and high_c >= c["anvil_high_cloud_pct"])
            wet = [pr[li] for li in near10 if pr[li] is not None]
            if mid_c is not None and wet:
                rec["tk"] = int(mid_c >= c["thick_mid_cloud_pct"] and max(wet) >= c["precip_in"])
                rec["dw"] = int(mid_c >= c["thick_mid_cloud_pct"] and max(wet) >= c["moderate_precip_in"])
            clouds = [x for x in (mid_c, high_c) if x is not None]
            if clouds:
                rec["cl"] = int(max(clouds) >= c["debris_cloud_pct"])
            if rec:
                series[int(t)] = rec
        if series:
            out[mid] = series
    return out


def openmeteo_ens(scfg, ctx):
    path = os.path.join(ctx.root, "cache", f"om_{scfg['id']}.json")
    try:
        with open(path) as f:
            old = json.load(f)
        if ctx.now - old["fetched"] < ctx.om_refresh_h * 3600:
            mem = {m: {int(t): r for t, r in s.items() if ctx.in_window(int(t))} for m, s in old["members"].items()}
            return SourceResult(mem, cycle=time.strftime("%H:%MZ", time.gmtime(old["fetched"])) + " (cached)",
                                note=f"{len(mem)} members", status="ok" if mem else "missing")
    except (OSError, ValueError, KeyError):
        pass
    lat, lon = _center(ctx)
    pts = ring_points(lat, lon)
    mem = _members(_fetch(scfg["model"], ctx, pts), pts, ctx)
    if mem:
        with open(path, "w") as f:
            json.dump({"fetched": ctx.now, "members": mem}, f, separators=(",", ":"))
    return SourceResult(mem, cycle=time.strftime("%H:%MZ", time.gmtime(ctx.now)),
                        note=f"{len(mem)} members", status="ok" if mem else "missing")
