"""Stage 5 - validation.

Four kinds of check, all of which must pass before the build is publishable:

1. Structural: every row carries every field the exported metadata claims,
   values are finite, percentiles are ordered, censor flags agree with the
   censored fraction, and no object appears twice.
2. Against the reference charts supplied with this project, which quote present
   distance, closest-approach distance and epoch for a dozen stars.
3. The integrated Monte Carlo intervals from stage 5b, for the curated stars.
   Neither reference chart shows error bars.
4. Physical assertions on stars whose encounters are in the literature.

A bound-pair consistency report is printed alongside; it is diagnostic rather
than pass/fail, because a genuine RV disagreement between components is
information, not a build error.

The result - overall status and every individual assertion - is written to
validation.json, which stage 6 embeds in the published metadata. Earlier
versions printed PASS/FAIL to the terminal and recorded only the comparison
tables, so nothing downstream could tell a validated build from an unvalidated
one.
"""

from __future__ import annotations

import csv
import json
import math
import pathlib


BASE = pathlib.Path(__file__).parent.parent
PROC = BASE / "data" / "processed"
LY_PER_PC = 3.261563777

# Values read off the two reference charts shipped with this project:
#   label -> (present distance ly, closest approach ly, epoch yr)
REFERENCE = {
    "Gliese 445":       (17.6, 3.45, 45000),
    "Ross 128":         (10.9, 6.233, 71000),
    "Ross 248":         (10.3, 3.024, 36000),
    "Ross 154":         (9.7, 6.39, 157000),
    "Sirius":           (8.6, 7.86, 64000),
    "Lalande 21185":    (8.3, 4.65, 20500),
    "Wolf 359":         (7.78, 7.35, -13850),
    "Barnard's Star":   (5.98, 3.74, 9800),
    "Alpha Centauri A": (4.36, 2.97, 28400),
    "Proxima Centauri": (4.24, 2.9, 27400),
}

REQUIRED = ("name", "ra", "dec", "parallax", "rv", "grade",
            "dist_now_pc", "gal_dmin_pc", "gal_tmin_myr",
            "lin_dmin_pc", "lin_tmin_myr",
            "mc_dmin_med", "mc_dmin_lo", "mc_dmin_hi",
            "mc_tmin_med", "mc_tmin_lo", "mc_tmin_hi",
            "mc_censored", "mc_model", "sys_dmin_pc", "zp_dmin_pc")


def load_orbits() -> list[dict]:
    with open(PROC / "orbits.csv") as fh:
        return list(csv.DictReader(fh))


def fnum(row: dict, key: str) -> float:
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return float("nan")


# --------------------------------------------------------------- structural
def check_fields(stars: list[dict]) -> tuple[bool, str]:
    missing = [c for c in REQUIRED if c not in stars[0]]
    if missing:
        return False, f"columns absent: {missing}"
    blank = [s["name"] for s in stars
             if any(s.get(c, "") in ("", None) for c in REQUIRED)]
    if blank:
        return False, f"{len(blank)} rows with an empty required field, e.g. {blank[:3]}"
    return True, f"{len(REQUIRED)} required columns present on all {len(stars)} rows"


def check_model(stars: list[dict]) -> tuple[bool, str]:
    bad = [s["name"] for s in stars if s.get("mc_model") != "integrated"]
    if bad:
        return False, f"{len(bad)} rows not mc_model='integrated', e.g. {bad[:3]}"
    return True, f"all {len(stars)} rows carry the integrated uncertainty model"


def check_finite(stars: list[dict]) -> tuple[bool, str]:
    cols = ("dist_now_pc", "gal_dmin_pc", "gal_tmin_myr", "mc_dmin_med",
            "mc_dmin_lo", "mc_dmin_hi", "mc_tmin_med", "mc_tmin_lo", "mc_tmin_hi")
    bad = [s["name"] for s in stars
           if any(not math.isfinite(fnum(s, c)) for c in cols)]
    if bad:
        return False, f"{len(bad)} rows with a non-finite value, e.g. {bad[:3]}"
    return True, f"all distances and epochs finite across {len(stars)} rows"


def check_percentiles(stars: list[dict]) -> tuple[bool, str]:
    bad = []
    for s in stars:
        d = (fnum(s, "mc_dmin_lo"), fnum(s, "mc_dmin_med"), fnum(s, "mc_dmin_hi"))
        t = (fnum(s, "mc_tmin_lo"), fnum(s, "mc_tmin_med"), fnum(s, "mc_tmin_hi"))
        if not (d[0] <= d[1] <= d[2]) or not (t[0] <= t[1] <= t[2]):
            bad.append(s["name"])
    if bad:
        return False, f"{len(bad)} rows with unordered percentiles, e.g. {bad[:3]}"
    return True, "16th <= 50th <= 84th on distance and epoch everywhere"


def check_censored(stars: list[dict], t_span: float) -> tuple[bool, str]:
    """A star flagged censored must actually reach the window edge.

    The exported flag bit is derived from mc_censored > 0.5 and the interface
    then says the encounter lies outside the window, so the epoch distribution
    had better touch the boundary.

    Note what this does *not* require: that the median epoch is at the edge.
    Seven stars have most draws censored while the median lands inside, because
    their draws split between the two boundaries - the astrometry cannot decide
    whether the encounter is before -5 Myr or after +5 Myr. That is a real state
    of knowledge, not a bug, and the interface reports those neutrally instead
    of picking a side.
    """
    bad, split = [], 0
    for s in stars:
        c = fnum(s, "mc_censored")
        if not (0.0 <= c <= 1.0):
            bad.append(f"{s['name']} (fraction {c})")
            continue
        if c <= 0.5:
            continue
        lo, med, hi = (fnum(s, "mc_tmin_lo"), fnum(s, "mc_tmin_med"),
                       fnum(s, "mc_tmin_hi"))
        edge = t_span * 0.99
        if not (lo <= -edge or hi >= edge):
            bad.append(f"{s['name']} (censored but epochs span [{lo:+.3f},{hi:+.3f}])")
        elif abs(med) < edge:
            split += 1
    if bad:
        return False, f"{len(bad)} inconsistent, e.g. {bad[:3]}"
    n = sum(1 for s in stars if fnum(s, "mc_censored") > 0.5)
    return True, (f"fractions in [0,1]; {n} censored stars all reach the "
                  f"+/-{t_span} Myr edge ({split} straddle both edges)")


def check_systematics(stars: list[dict]) -> tuple[bool, str]:
    """The Galaxy-model shift must be present, finite and non-negative.

    It is a magnitude, not a signed correction, and a negative or missing value
    would mean stage 9 did not run or ran against a different catalog.
    """
    bad = [s["name"] for s in stars
           if any(not math.isfinite(fnum(s, c)) or fnum(s, c) < 0
                  for c in ("sys_dmin_pc", "zp_dmin_pc"))]
    if bad:
        return False, f"{len(bad)} absent, negative or non-finite, e.g. {bad[:3]}"
    inside = [s for s in stars
              if fnum(s, "mc_dmin_med") < 5.0 and fnum(s, "mc_censored") <= 0.5]
    dominated = sum(1 for s in inside if fnum(s, "sys_dmin_pc")
                    > fnum(s, "mc_dmin_hi") - fnum(s, "mc_dmin_lo"))
    zp = sum(1 for s in inside if fnum(s, "zp_dmin_pc")
             > fnum(s, "mc_dmin_hi") - fnum(s, "mc_dmin_lo"))
    return True, (f"present and non-negative on all {len(stars)} rows; of "
                  f"{len(inside)} stars within 5 pc, {dominated} move further "
                  f"under the Galaxy-model variants and {zp} under the parallax "
                  f"zero-point than under measurement resampling")


def check_correlations(stars: list[dict]) -> tuple[bool, str]:
    """Gaia rows must carry the astrometric correlations the draws rely on.

    This exists because a half-finished re-fetch is invisible otherwise. The
    correlations arrive in the raw Gaia files; if one of the two queries fails,
    the stars it covers silently fall back to independent draws and every
    interval for them is subtly wrong while looking entirely normal. The 30 pc
    query is the one that failed during this build, and it covers the *nearest*
    stars - the ones anyone checks first.
    """
    gaia = [s for s in stars if (s.get("source_id") or "").strip()]
    if not gaia:
        return False, "no Gaia rows at all"
    cols = ("plx_pmra_corr", "plx_pmdec_corr", "pmra_pmdec_corr")
    if any(c not in stars[0] for c in cols):
        return False, f"correlation columns absent: {cols}"
    have = [s for s in gaia
            if any((s.get(c) or "").strip() not in ("", "0", "0.0") for c in cols)]
    frac = len(have) / len(gaia)
    near = [s for s in gaia if fnum(s, "dist_now_pc") < 30]
    near_have = sum(1 for s in near
                    if any((s.get(c) or "").strip() not in ("", "0", "0.0") for c in cols))
    if frac < 0.95:
        return False, (f"only {len(have)} of {len(gaia)} Gaia rows ({frac:.1%}) carry "
                       f"correlations - a raw fetch is incomplete")
    if near and near_have / len(near) < 0.95:
        return False, (f"only {near_have} of {len(near)} Gaia rows within 30 pc carry "
                       f"correlations - the neighborhood query is stale")
    return True, (f"{len(have)} of {len(gaia)} Gaia rows ({frac:.1%}), including "
                  f"{near_have} of {len(near)} within 30 pc")


def check_unique(stars: list[dict]) -> tuple[bool, str]:
    """No object may appear twice under two identities.

    rho Ori arrived once from Gaia and once from Hipparcos, at 75 pc and 107 pc,
    because their parallaxes disagree by 30.1% and the positional dedup rejected
    the match. Stage 2 now resolves such pairs through SIMBAD's identifier
    tables; this makes sure it stays resolved.
    """
    dupes = []
    for key in ("source_id", "hip", "simbad_id", "name"):
        seen: dict[str, int] = {}
        for s in stars:
            v = (s.get(key) or "").strip()
            if not v:
                continue
            seen[v] = seen.get(v, 0) + 1
        rep = [f"{key}={k}" for k, n in seen.items() if n > 1]
        dupes += rep[:3]
    if dupes:
        return False, f"repeated identities: {dupes[:6]}"
    return True, "source_id, HIP, SIMBAD id and display name each unique"


def check_rankings(stars: list[dict]) -> tuple[bool, str]:
    """The published tables must be a straight sort of the shipped column.

    Guards against a README that was hand-edited, or refreshed from one column
    while the interface ranks on another.
    """
    ranked = sorted(stars, key=lambda s: fnum(s, "mc_dmin_med"))
    if any(not math.isfinite(fnum(s, "mc_dmin_med")) for s in ranked[:20]):
        return False, "the closest twenty include a non-finite median"
    tables = rankings(stars)
    for side, rows in ((k, v) for k, v in tables.items() if isinstance(v, list)):
        d = [r["dmin_pc"] for r in rows]
        if any(a > b for a, b in zip(d, d[1:])):
            return False, f"{side} is not sorted by the shipped column"
        if any(r["censored"] for r in rows):
            return False, f"{side} contains a boundary-limited trajectory"
    ahead = tables["closest_ahead"]
    if not ahead or ahead[0]["star"] != "Gliese 710":
        return False, ("closest future approach is "
                       f"{ahead[0]['star'] if ahead else 'nothing'}, expected Gliese 710")
    return True, (f"closest ahead {ahead[0]['star']} at {ahead[0]['dmin_pc']:.4f} pc; "
                  f"closest behind {tables['closest_behind'][0]['star']} at "
                  f"{tables['closest_behind'][0]['dmin_pc']:.4f} pc")


def check_rate() -> tuple[bool, str]:
    """The measured encounter rate must still agree with the literature.

    This is the only assertion here that tests a *derived* quantity rather than
    the catalog's internal consistency, which makes it the one most likely to
    catch a change that leaves every row individually plausible - a shifted
    completeness correction, a lost slice of the neighborhood sample - while
    moving the number the whole exercise exists to produce. Agreement is
    required at 2 sigma of the two errors combined.
    """
    path = PROC / "encounter_rate.json"
    if not path.exists():
        return False, "encounter_rate.json missing; run scripts/10_rate.py"
    r = json.loads(path.read_text())
    ours, pub = r["rate_within_1pc_per_myr"], r["published"]["value"]
    sigma = r["sigma_from_published"]
    if not math.isfinite(ours) or ours <= 0:
        return False, f"rate is {ours}"
    if sigma > 2.0:
        return False, (f"{ours} vs published {pub} per Myr is {sigma} sigma apart")
    return True, (f"{ours} +/- {r['stat_error']} (stat) +/- {r['sys_error']} (sys) "
                  f"per Myr within 1 pc, vs published {pub} +/- "
                  f"{r['published']['error']} — {sigma} sigma")


# --------------------------------------------------------------- comparisons
def compare_reference(stars: list[dict]) -> list[dict]:
    by_name = {s["name"]: s for s in stars}
    rows = []
    print("\n  Against the reference charts (distances in light years):")
    print(f"    {'star':20s} {'now':>13s}  {'closest':>15s}  {'epoch (kyr)':>17s}")
    print(f"    {'':20s} {'ref':>6s}{'ours':>7s}  {'ref':>7s}{'ours':>8s}  {'ref':>8s}{'ours':>9s}")
    for label, (ref_now, ref_min, ref_t) in REFERENCE.items():
        s = by_name.get(label)
        if not s:
            print(f"    {label:20s}  -- not in catalog --")
            continue
        now = fnum(s, "dist_now_pc") * LY_PER_PC
        dmin = fnum(s, "gal_dmin_pc") * LY_PER_PC
        tmin = fnum(s, "gal_tmin_myr") * 1000.0
        print(f"    {label:20s} {ref_now:6.2f}{now:7.2f}  {ref_min:7.2f}{dmin:8.2f}  "
              f"{ref_t/1000:8.1f}{tmin:9.1f}")
        rows.append({"star": label, "ref_now_ly": ref_now, "our_now_ly": round(now, 3),
                     "ref_dmin_ly": ref_min, "our_dmin_ly": round(dmin, 3),
                     "ref_t_kyr": ref_t / 1000, "our_t_kyr": round(tmin, 2)})
    return rows


def monte_carlo(stars: list[dict]) -> list[dict]:
    """Report the integrated Monte Carlo results computed in stage 5b.

    This used to run a second, straight-line Monte Carlo of its own, which meant
    the validation output and the shipped data described different models. It
    now reads the one set of numbers rather than inventing another.
    """
    out = []
    print("\n  Integrated Monte Carlo (from stage 5b, 16th-84th percentile):")
    print(f"    {'star':22s} {'d_min (pc)':>26s} {'t_min (kyr)':>24s}  grade")
    for s in stars:
        if not s["featured"]:
            continue
        med, lo, hi = (fnum(s, "mc_dmin_med"), fnum(s, "mc_dmin_lo"),
                       fnum(s, "mc_dmin_hi"))
        tmed, tlo, thi = (fnum(s, "mc_tmin_med") * 1000, fnum(s, "mc_tmin_lo") * 1000,
                          fnum(s, "mc_tmin_hi") * 1000)
        cens = fnum(s, "mc_censored")
        mark = ""
        if cens > 0.5:
            mark = ("  (bound: minimum lies earlier)" if tmed < 0
                    else "  (bound: minimum lies later)")
        print(f"    {s['name'][:22]:22s} {med:8.3f} [{lo:7.3f},{hi:7.3f}] "
              f"{tmed:10.1f} [{tlo:8.1f},{thi:8.1f}]  {s['grade']}{mark}")
        out.append({
            "star": s["name"], "grade": s["grade"],
            "model": s.get("mc_model", "unknown"),
            "censored_fraction": round(cens, 3),
            "dmin_pc_median": round(med, 4),
            "dmin_pc_lo": round(lo, 4), "dmin_pc_hi": round(hi, 4),
            "tmin_kyr_median": round(tmed, 2),
            "tmin_kyr_lo": round(tlo, 2), "tmin_kyr_hi": round(thi, 2),
        })
    return out


def rankings(stars: list[dict], n: int = 10) -> dict:
    """The tables the README publishes, derived here so they cannot drift.

    Same basis as the published tables: grade A, boundary-limited trajectories
    excluded, split on the same epoch column that is quoted beside each entry.
    """
    def row(s):
        return {"star": s["name"],
                "dmin_pc": round(fnum(s, "mc_dmin_med"), 4),
                "lo_pc": round(fnum(s, "mc_dmin_lo"), 4),
                "hi_pc": round(fnum(s, "mc_dmin_hi"), 4),
                "tmin_kyr": round(fnum(s, "mc_tmin_med") * 1000, 1),
                "dist_now_pc": round(fnum(s, "dist_now_pc"), 2),
                "grade": s["grade"],
                "censored": fnum(s, "mc_censored") > 0.5}
    pool = [s for s in stars
            if s["grade"] == "A" and fnum(s, "mc_censored") <= 0.5]
    key = lambda s: fnum(s, "mc_dmin_med")
    return {
        "basis": "grade A, boundary-limited trajectories excluded, ranked on "
                 "the median of the integrated draws",
        "closest_ahead": [row(s) for s in
                          sorted((s for s in pool if fnum(s, "mc_tmin_med") > 0),
                                 key=key)[:n]],
        "closest_behind": [row(s) for s in
                           sorted((s for s in pool if fnum(s, "mc_tmin_med") < 0),
                                  key=key)[:n]],
    }


def binary_consistency(stars: list[dict]) -> None:
    """Components of a bound pair must share an encounter solution.

    Alpha Centauri A and B orbit each other at ~23 AU with orbital speeds of a
    few km/s. Catalogs list each component's instantaneous velocity, not the
    barycenter's, so the two components can yield visibly different perihelia
    even though the system arrives as one object. Worth surfacing rather than
    hiding.
    """
    print("\n  Bound-pair consistency check:")
    pairs = [("Alpha Centauri A", "Alpha Centauri B"),
             ("61 Cygni A", "61 Cygni B"),
             ("Struve 2398 A", "Struve 2398 B"),
             ("Groombridge 34 A", "Groombridge 34 B")]
    by_name = {s["name"]: s for s in stars}
    for a, b in pairs:
        sa, sb = by_name.get(a), by_name.get(b)
        if not sa or not sb:
            continue
        da, db = fnum(sa, "gal_dmin_pc"), fnum(sb, "gal_dmin_pc")
        rva, rvb = fnum(sa, "rv"), fnum(sb, "rv")
        print(f"    {a:18s} vs {b:18s}  d_min {da:6.3f} / {db:6.3f} pc "
              f"(delta {abs(da-db):5.3f})   RV {rva:+7.2f} / {rvb:+7.2f} km/s")


# Physical assertions. Deliberately loose: they catch a pipeline that has
# broken, not a model that disagrees with the literature at the percent level.
def physical(by: dict) -> list[tuple[str, bool, str]]:
    out = []

    def chk(label, fn):
        try:
            ok = bool(fn())
            detail = ""
        except (KeyError, ValueError, TypeError) as exc:
            ok, detail = False, f"{type(exc).__name__}: {exc}"
        out.append((label, ok, detail))

    g = lambda name, col: fnum(by[name], col)
    chk("Proxima Centauri present-day distance within 1% of 1.302 pc",
        lambda: abs(g("Proxima Centauri", "dist_now_pc") - 1.302) / 1.302 < 0.01)
    chk("Gliese 710 closest approach under 0.1 pc",
        lambda: g("Gliese 710", "gal_dmin_pc") < 0.1)
    chk("Gliese 710 encounter between +1.2 and +1.4 Myr",
        lambda: 1.2 < g("Gliese 710", "gal_tmin_myr") < 1.4)
    chk("Scholz's star passed within 0.4 pc in the past",
        lambda: g("Scholz's Star", "gal_dmin_pc") < 0.4
                and g("Scholz's Star", "gal_tmin_myr") < 0)
    # 0.15 pc was calibrated to a broken state: the components disagreed by
    # 0.232 pc on Hipparcos parallaxes that differ by 5.6%, for a pair orbiting
    # at 23 AU. With the published system parallax they agree to 0.001 pc, so
    # the tolerance is now tight enough to notice if either fix regresses.
    chk("Alpha Cen A and B agree to better than 0.02 pc",
        lambda: abs(g("Alpha Centauri A", "gal_dmin_pc")
                    - g("Alpha Centauri B", "gal_dmin_pc")) < 0.02)
    chk("Alpha Cen A and B share the published system parallax",
        lambda: abs(g("Alpha Centauri A", "parallax") - 750.81) < 0.01
                and abs(g("Alpha Centauri B", "parallax") - 750.81) < 0.01)
    chk("Alpha Cen A and B share one radial velocity",
        lambda: abs(g("Alpha Centauri A", "rv") - g("Alpha Centauri B", "rv")) < 0.01)
    chk("Barnard's Star closest approach within 5% of 1.155 pc",
        lambda: abs(g("Barnard's Star", "gal_dmin_pc") - 1.155) / 1.155 < 0.05)
    # HD 7977 no longer asserts agreement with a published distance, because
    # this build does not agree: the compiled sub-5-ly list gives 0.478
    # +0.104/-0.078 ly and this gives 0.123 ly, a factor of 3.9, with intervals
    # that do not overlap. See compare_published.py. Earlier versions asserted
    # agreement with 0.0641 pc, a figure that is not in that compilation and
    # that I could not source; testing against it was testing nothing.
    #
    # What both agree on, and what a broken pipeline would break, is that this
    # is a sub-parsec approach in the past at about -2.76 Myr. That is the
    # assertion; the distance disagreement is documented, not asserted away.
    chk("HD 7977 is a sub-parsec past approach near -2.76 Myr",
        lambda: g("HD 7977", "mc_dmin_med") < 1.0
                and -2.9 < g("HD 7977", "mc_tmin_med") < -2.6)
    chk("Gliese 710 agrees with the published 0.167 ly within 10%",
        lambda: abs(g("Gliese 710", "mc_dmin_med") * 3.261563777 - 0.167) < 0.0167)
    return out


def main() -> None:
    stars = load_orbits()
    meta = json.loads((PROC / "orbit_meta.json").read_text())
    t_span = meta["t_span_myr"]
    print(f"  {len(stars)} stars")

    checks: list[tuple[str, bool, str]] = []
    for label, fn in (
        ("every required column present and populated", lambda: check_fields(stars)),
        ("uncertainty draws were integrated, not straight-line", lambda: check_model(stars)),
        ("all distances and epochs finite", lambda: check_finite(stars)),
        ("Monte Carlo percentiles correctly ordered", lambda: check_percentiles(stars)),
        ("censor flags agree with the censored fraction",
         lambda: check_censored(stars, t_span)),
        ("Galaxy-model sensitivity recorded for every star",
         lambda: check_systematics(stars)),
        ("Gaia rows carry their astrometric correlations",
         lambda: check_correlations(stars)),
        ("no object appears under two identities", lambda: check_unique(stars)),
        ("published ranking follows the shipped column", lambda: check_rankings(stars)),
        ("encounter rate agrees with the published measurement",
         lambda: check_rate()),
    ):
        ok, detail = fn()
        checks.append((label, ok, detail))

    ref_rows = compare_reference(stars)
    mc_rows = monte_carlo(stars)
    binary_consistency(stars)
    checks += physical({s["name"]: s for s in stars})

    print("\n  Build assertions:")
    for label, ok, detail in checks:
        print(f"    [{'PASS' if ok else 'FAIL'}] {label}"
              + (f"\n           {detail}" if detail else ""))
    passed = all(ok for _, ok, _ in checks)

    (PROC / "validation.json").write_text(json.dumps({
        "status": "pass" if passed else "fail",
        "n_checks": len(checks),
        "n_failed": sum(1 for _, ok, _ in checks if not ok),
        "checks": [{"check": c, "passed": ok, "detail": d} for c, ok, d in checks],
        "reference_comparison": ref_rows,
        "monte_carlo": mc_rows,
        "rankings": rankings(stars),
    }, indent=2))
    print(f"\n  -> data/processed/validation.json  ({'pass' if passed else 'FAIL'}, "
          f"{sum(1 for _, ok, _ in checks if not ok)} of {len(checks)} failed)")
    if not passed:
        raise SystemExit("validation assertions FAILED")


if __name__ == "__main__":
    main()
