"""Hi-res NCEP guidance turned into per-pad LLCC evidence for each member-hour.

Record keys are "<rule>@<pad>" (0/1, or 0-1 probability for NBM):
    lt   Lightning (4.1.1): lightning field within 10 nmi of the flight path
    cu   Cumulus (4.1.3): -10 C-level reflectivity >= cu_dbz within 5 nmi,
         or echo top above the -20 C height within 10 nmi
    an   Anvil (4.1.4-4.1.5, rough): deep convection within 20 nmi + high cloud within 10 nmi
    tk   Thick cloud layers (4.1.8, rough): mid cloud near overcast + echo >= 7.5 dBZ within 5 nmi
    dw   Disturbed weather (4.1.7, rough): >= 30 dBZ within 5 nmi under an overcast mid deck
Helpers for the debris rule (worked out in run.py once each member's series is assembled):
    d10  deep convection within 10 nmi;  cl   mid/high cloud within 5 nmi
"""
from __future__ import annotations

import re
import time
from concurrent.futures import ThreadPoolExecutor

from .common import Context, SourceResult, floor_hour, iso, log
from .grib import Missing, exists, fetch, fetch_select

RADII = [5.0, 10.0, 20.0]


def _fmt(tmpl: str, cycle: int, fh: int) -> str:
    g = time.gmtime(cycle)
    return tmpl.format(ymd=time.strftime("%Y%m%d", g), hh=f"{g.tm_hour:02d}", fh=fh)


def _area(url: str, cycle: int, ctx: Context) -> dict | None:
    cached = ctx.cache.get(url)
    if cached is not None:
        return cached
    try:
        vals = fetch(url, ctx.pads, RADII)
    except Missing:
        return None
    out = {f: {pid: {str(int(r)): v for r, v in byr.items()} for pid, byr in d.items()}
           for f, d in vals.items() if isinstance(d, dict)}
    if out:
        ctx.cache.put(url, cycle, out)
    return out


def _get(p, field, pid, r):
    try:
        return p[field][pid][str(int(r))]
    except (KeyError, TypeError):
        return None


def record(p: dict, ctx: Context) -> dict | None:
    if not p:
        return None
    c = ctx.rules
    rec = {}
    for pad in ctx.pads:
        pid = pad["id"]
        g = lambda f, r: _get(p, f, pid, r)
        lt10 = g("ltng", 10)
        if lt10 is not None:
            rec[f"lt@{pid}"] = int(lt10 > c["ltng_threshold"])
        rd5, et10 = g("refd263", 5), g("retop", 10)
        if rd5 is not None or et10 is not None:
            rec[f"cu@{pid}"] = int((rd5 is not None and rd5 >= c["cu_dbz"]) or
                                   (et10 is not None and et10 >= c["h20_m"]))

        def deep(r):
            et, lt = g("retop", r), g("ltng", r)
            if et is None and lt is None:
                return None
            return int((et is not None and et >= c["h20_m"]) or (lt is not None and lt > c["ltng_threshold"]))

        d20, d10 = deep(20), deep(10)
        hc10, mc5, hc5, rc5 = g("hcdc", 10), g("mcdc", 5), g("hcdc", 5), g("refc", 5)
        if d20 is not None and hc10 is not None:
            rec[f"an@{pid}"] = int(bool(d20) and hc10 >= c["anvil_high_cloud_pct"])
        if mc5 is not None and rc5 is not None:
            rec[f"tk@{pid}"] = int(mc5 >= c["thick_mid_cloud_pct"] and rc5 >= 7.5)
            rec[f"dw@{pid}"] = int(rc5 >= c["disturbed_dbz"] and mc5 >= c["thick_mid_cloud_pct"])
        if d10 is not None:
            rec[f"d10@{pid}"] = d10
        clouds = [x for x in (mc5, hc5) if x is not None]
        if clouds:
            rec[f"cl@{pid}"] = int(max(clouds) >= c["debris_cloud_pct"])
    return rec or None


def _pick(bases, tmpl, cycles):
    for base in bases:
        for c in cycles:
            if exists(f"{base}/{_fmt(tmpl, c, 1)}.idx"):
                return base, c
    return None


def _run_tasks(scfg, ctx, tasks, workers):
    def run(task):
        mid, url, c, valid = task
        try:
            return task, record(_area(url, c, ctx), ctx)
        except Exception as e:
            log.warning("%s %s: %s", scfg["id"], url.rsplit("/", 1)[-1], e)
            return task, None

    members: dict[str, dict] = {}
    with ThreadPoolExecutor(workers) as ex:
        for (mid, _, _, valid), rec in ex.map(run, tasks):
            if rec:
                members.setdefault(mid, {})[valid] = rec
    return members


def tle(scfg: dict, ctx: Context) -> SourceResult:
    recent = [floor_hour(ctx.now) - k * 3600 for k in range(10)]
    picked = _pick([scfg["base"]], scfg["file"], recent)
    if not picked:
        return SourceResult({}, status="missing", note="no recent cycle found")
    base, latest = picked
    cycles = [latest - k * 3600 for k in range(int(scfg.get("lag_cycles", 6)))]
    cycles += [c for c in (latest - k * 3600 for k in range(30))
               if time.gmtime(c).tm_hour % 6 == 0 and c not in cycles][: int(scfg.get("synoptic_extra", 2))]
    tasks = []
    for c in cycles:
        mx = scfg.get("long_fh", 48) if time.gmtime(c).tm_hour % 6 == 0 else scfg.get("short_fh", 18)
        for fh in range(0, mx + 1):
            if ctx.in_window(c + fh * 3600):
                tasks.append((time.strftime("%d/%HZ", time.gmtime(c)), f"{base}/{_fmt(scfg['file'], c, fh)}", c, c + fh * 3600))
    members = _run_tasks(scfg, ctx, tasks, 16)
    return SourceResult(members, cycle=iso(latest), note=f"{len(members)}/{len(cycles)} cycles",
                        status="ok" if len(members) == len(cycles) else ("partial" if members else "missing"))


def multi_model(scfg: dict, ctx: Context) -> SourceResult:
    tasks, notes, missing = [], [], []
    for comp in scfg["components"]:
        hours = set(comp.get("cycles", [0, 6, 12, 18]))
        cands = [c for c in (floor_hour(ctx.now) - k * 3600 for k in range(60)) if time.gmtime(c).tm_hour in hours]
        picked = _pick(comp.get("bases") or [comp["base"]], comp["file"], cands)
        if not picked:
            missing.append(comp["id"])
            continue
        base, latest = picked
        use = [c for c in cands if c <= latest][: int(comp.get("lag", 2))]
        notes.append(f"{comp['id']} " + ", ".join(time.strftime("%HZ", time.gmtime(c)) for c in use))
        for c in use:
            for fh in range(0, int(comp.get("max_fh", 48)) + 1):
                if ctx.in_window(c + fh * 3600):
                    tasks.append((f"{comp['id']} {time.strftime('%d/%HZ', time.gmtime(c))}",
                                  f"{base}/{_fmt(comp['file'], c, fh)}", c, c + fh * 3600))
    members = _run_tasks(scfg, ctx, tasks, int(scfg.get("workers", 4)))
    note = "; ".join(notes + ([f"missing: {', '.join(missing)}"] if missing else []))
    return SourceResult(members, note=note,
                        status="missing" if not members else ("partial" if missing else "ok"))


# ---------------------------------------------------------------- NBM thunder probability
_TSTM = re.compile(r":TSTM:surface:(\d+)-(\d+) hour")


def nbm_prob(scfg: dict, ctx: Context) -> SourceResult:
    """NBM 1/3/6 h thunderstorm probability, highest value within 10 nmi of each flight
    path, as one member whose lt@pad values are probabilities (0-1)."""
    cands = [c for c in (floor_hour(ctx.now) - k * 3600 for k in range(36))
             if time.gmtime(c).tm_hour in set(scfg.get("cycles", [0, 6, 12, 18]))]
    picked = _pick([scfg["base"]], scfg["file"], cands)
    if not picked:
        return SourceResult({}, status="missing", note="no recent cycle found")
    base, cycle = picked
    maxdur = int(scfg.get("max_period_h", 6))

    def select(inv):
        out = []
        for rec in inv:
            m = _TSTM.search(rec[3])
            if m and "prob" in rec[3].lower():
                a, b = int(m.group(1)), int(m.group(2))
                if 0 < b - a <= maxdur:
                    out.append((f"t{a}-{b}", rec))
        return out

    def run(fh):
        url = f"{base}/{_fmt(scfg['file'], cycle, fh)}"
        key = url + "#llcc"
        got = ctx.cache.get(key)
        if got is None:
            try:
                vals = fetch_select(url, select, ctx.pads, [10.0])
            except Missing:
                return None
            except Exception as e:
                log.warning("%s f%03d: %s", scfg["id"], fh, e)
                return None
            got = {k: {pid: byr.get(10.0) for pid, byr in d.items()} for k, d in vals.items() if isinstance(d, dict)}
            ctx.cache.put(key, cycle, got)
        return got

    fhs = [fh for fh in range(1, int(scfg.get("max_fh", 192)) + 1) if ctx.in_window(cycle + fh * 3600)]
    series: dict[int, dict] = {}
    with ThreadPoolExecutor(16) as ex:
        for got in ex.map(run, fhs):
            for period, bypad in (got or {}).items():
                a, b = map(int, period[1:].split("-"))
                for h in range(a + 1, b + 1):
                    t = cycle + h * 3600
                    if not ctx.in_window(t):
                        continue
                    rec = series.setdefault(t, {})
                    for pid, v in bypad.items():
                        if v is not None:
                            k = f"lt@{pid}"
                            rec[k] = max(rec.get(k, 0.0), min(1.0, max(0.0, v / 100.0)))
    series = {t: r for t, r in series.items() if r}
    return SourceResult({"m00": series} if series else {}, cycle=iso(cycle),
                        note="thunder probability, highest within 10 nmi of the flight path",
                        status="ok" if series else "missing")
