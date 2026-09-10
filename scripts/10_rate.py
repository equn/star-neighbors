"""Stage 10 - measure the stellar encounter rate.

Every other stage produces a catalog: this star, that perihelion. This one
produces a *number* - how often any star passes within a given distance of the
Sun - which is the quantity the encounter literature actually reports, because
it is what sets the rate of Oort cloud perturbation. It can be compared against
a published measurement, and it is, at the bottom of this file.

The method rests on one physical assumption, and the data are made to test it
rather than asked to assume it.

    Over +/-5 Myr the solar neighborhood is statistically steady. The Sun
    travels about 100 pc; the local stellar density and velocity distribution
    do not measurably change. So the true rate of encounters is CONSTANT IN
    TIME, and a histogram of perihelion epochs must be flat.

It is not flat. It falls off steeply with |t_ph|, and the reason is visible in
the same table: the median present-day distance of an encountering star climbs
from 27 pc for encounters happening now to 204 pc for encounters 4-5 Myr away.
A star whose perihelion is far from the present is far from the Sun today, and
distant stars are the ones missing a radial velocity and therefore missing from
this catalog entirely.

That turns the histogram into a measurement instead of a nuisance. Near t = 0
the encountering stars are nearby, the sample is nearly complete, and the rate
there is very close to the true one. Far from t = 0 the deficit is
incompleteness and nothing else. So:

  1. Find the window around t = 0 where the rate stops rising as the window
     shrinks. That plateau is the complete regime. (It is |t| < ~0.375 Myr;
     `--scan` prints the evidence.)
  2. Inside it, weight each encounter by 1 / C(d_now), where C is this
     catalog's completeness against the Gaia census at that distance -
     the standard inverse-probability estimator.
  3. Fit N = k d_ph^2, which is the scaling a homogeneous flux of stars must
     obey, and read off the rate within 1 pc.

The d^2 fit is not decoration. It is a check on the whole calculation: if the
counts did not scale as the cross-sectional area, the sample would not be
behaving like a uniform flux of stars and the rate would mean nothing.

    python3 scripts/10_rate.py [--scan] [--json]
"""

from __future__ import annotations

import csv
import json
import math
import pathlib
import sys

BASE = pathlib.Path(__file__).parent.parent
PROC = BASE / "data" / "processed"
RAW = BASE / "data" / "raw"

# The plateau window. Chosen from --scan, which shows the rate flat to within
# its Poisson error for |t| < 0.375 Myr and falling steadily beyond.
WINDOW_MYR = 0.25

# Perihelion cuts the d^2 fit runs over. The lower end is where the statistics
# run out; the upper end is where residual incompleteness starts to bite, and
# the spread between fits is reported as the dominant systematic.
FIT_CUTS = (1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0)

# Distance bins for the completeness calibration, pc.
COMP_BINS = ((0, 5), (5, 10), (10, 15), (15, 20), (20, 25), (25, 30))

# Bailer-Jones et al. (2018), A&A 616, A37: "the present rate of encounters
# within 1 pc to be 19.7 +/- 2.2 per Myr", completeness-corrected, from Gaia
# DR2 over the same +/-5 Myr window.
PUBLISHED = (19.7, 2.2)

KM_S_TO_PC_MYR = 1.0227121650537077


def fnum(row: dict, key: str):
    try:
        v = float(row[key])
        return v if v == v else None
    except (TypeError, ValueError, KeyError):
        return None


def load_encounters() -> list[dict]:
    """Uncensored encounters: a real perihelion inside the window."""
    out = []
    for r in csv.DictReader(open(PROC / "orbits.csv")):
        d, t = fnum(r, "mc_dmin_med"), fnum(r, "mc_tmin_med")
        dn = fnum(r, "dist_now_pc")
        if d is None or t is None or dn is None:
            continue
        if (fnum(r, "mc_censored") or 0.0) > 0.5:
            continue
        out.append({"d": d, "t": abs(t), "dn": dn,
                    "lo": fnum(r, "mc_dmin_lo"), "hi": fnum(r, "mc_dmin_hi"),
                    "gal": fnum(r, "gal_dmin_pc"), "lin": fnum(r, "lin_dmin_pc")})
    return out


def completeness() -> tuple[dict, list]:
    """What fraction of the Gaia neighborhood census reaches this catalog.

    Deliberately measured against the census rather than against Gaia's
    radial-velocity column: this pipeline supplements Gaia RVs with SIMBAD
    ones, so Gaia's own RV fraction understates the coverage by ten to twenty
    points and would inflate the corrected rate by the same amount.
    """
    have = {r["source_id"] for r in csv.DictReader(open(PROC / "orbits.csv"))
            if r.get("source_id")}
    rows = list(csv.DictReader(open(RAW / "gaia_30pc.csv")))
    table, C = [], {}
    for lo, hi in COMP_BINS:
        sub = []
        for r in rows:
            p = fnum(r, "parallax")
            if p and p > 0 and lo <= 1000.0 / p < hi:
                sub.append(r)
        if not sub:
            continue
        n_rv = sum(1 for r in sub if fnum(r, "radial_velocity") is not None)
        n_our = sum(1 for r in sub if r["source_id"] in have)
        C[(lo, hi)] = n_our / len(sub)
        table.append({"bin_pc": [lo, hi], "census": len(sub),
                      "gaia_rv": n_rv, "ours": n_our,
                      "c_gaia_rv": round(n_rv / len(sub), 4),
                      "c_ours": round(n_our / len(sub), 4)})
    return C, table


def comp_at(C: dict, d_pc: float) -> float:
    for (lo, hi), c in C.items():
        if lo <= d_pc < hi:
            return c
    return C[COMP_BINS[-1]]          # beyond calibration, hold the last bin


def rate_within(enc, C, window, cuts, dkey="d"):
    """Inverse-completeness weighted fit of N = k d^2; returns rate at 1 pc."""
    span = 2.0 * window                       # elapsed time the window covers
    rows, num, den, n_raw = [], 0.0, 0.0, 0
    for c in cuts:
        n = w = 0
        for e in enc:
            d = e[dkey]
            if d is None or e["t"] >= window or d >= c:
                continue
            n += 1
            w += 1.0 / max(comp_at(C, e["dn"]), 0.05)
        rows.append({"d_ph_pc": c, "n": n, "raw_per_myr": round(n / span, 2),
                     "corrected_per_myr": round(w / span, 2),
                     "per_d2": round(w / span / (c * c), 2)})
        num += w
        den += c * c
        n_raw += n
    if den == 0 or n_raw == 0:
        return 0.0, 0.0, rows, 0
    rate = num / den / span
    return rate, rate / math.sqrt(n_raw), rows, n_raw


def scan(enc, C) -> None:
    """The evidence for the plateau, and for everything the rate depends on."""
    print("\n  rate(<1 pc) against the plateau window "
          "[flat while complete, falling once not]")
    print(f"    {'|t| <':>10}{'rate':>8}{'+/-':>6}{'N':>7}")
    for w in (0.125, 0.1875, 0.25, 0.375, 0.5, 0.75, 1.0):
        r, e, _, n = rate_within(enc, C, w, FIT_CUTS)
        print(f"    {w:8.3f}  {r:8.1f}{e:6.1f}{n:7d}")

    print("\n  ... against which perihelion column is used")
    print(f"    {'column':>16}{'rate':>8}{'+/-':>6}{'N':>7}")
    for k, lab in (("lo", "16th pct"), ("d", "median"), ("hi", "84th pct"),
                   ("gal", "nominal"), ("lin", "straight line")):
        r, e, _, n = rate_within(enc, C, WINDOW_MYR, FIT_CUTS, dkey=k)
        print(f"    {lab:>16}{r:8.1f}{e:6.1f}{n:7d}")

    print("\n  ... against the perihelion range the d^2 fit is taken over")
    print(f"    {'range':>16}{'rate':>8}{'+/-':>6}{'N':>7}")
    spread = []
    for dmax in (2.0, 3.0, 4.0, 5.0):
        cuts = tuple(c for c in FIT_CUTS if c <= dmax)
        r, e, _, n = rate_within(enc, C, WINDOW_MYR, cuts)
        spread.append(r)
        print(f"    {'1-' + str(int(dmax)) + ' pc':>16}{r:8.1f}{e:6.1f}{n:7d}")
    print(f"\n    the fit range is the dominant systematic: "
          f"{min(spread):.1f} to {max(spread):.1f} per Myr")


def kinetic_check(enc) -> dict:
    """Γ = π d² n <v>, from first principles, as an independent sanity check.

    Wholly separate from the counting above: it uses the present-day density
    and speed of the neighborhood and never looks at a perihelion. Its own
    uncertainty is large - the local density is known to maybe 30%, and the
    mean speed comes from the radial-velocity subsample, which is biased toward
    brighter and so kinematically colder stars - so it is quoted as a range and
    not used to correct anything.
    """
    rows = list(csv.DictReader(open(RAW / "gaia_30pc.csv")))
    d = sorted(1000.0 / fnum(r, "parallax") for r in rows
               if (fnum(r, "parallax") or 0) > 0)
    dens = {}
    for R in (5, 10, 20, 30):
        n = sum(1 for x in d if x <= R)
        dens[R] = n / (4.0 / 3.0 * math.pi * R ** 3)
    # mean heliocentric speed of the encountering sample
    speeds = []
    for r in csv.DictReader(open(PROC / "orbits.csv")):
        v = [fnum(r, k) for k in ("vx", "vy", "vz")]
        if all(x is not None for x in v) and (fnum(r, "dist_now_pc") or 1e9) <= 30:
            speeds.append(math.sqrt(sum(x * x for x in v)))
    v_mean = sum(speeds) / len(speeds) if speeds else 0.0
    out = {"density_per_pc3": {str(k): round(v, 4) for k, v in dens.items()},
           "mean_speed_km_s": round(v_mean, 2), "gamma_per_myr": {}}
    for n in (0.08, 0.10, 0.12):
        out["gamma_per_myr"][f"n={n}"] = round(
            math.pi * n * v_mean * KM_S_TO_PC_MYR, 1)
    return out


def main() -> None:
    enc = load_encounters()
    C, table = completeness()

    print(f"  uncensored encounters in the catalog: {len(enc)}")
    print("\n  completeness of this catalog against the Gaia 30 pc census")
    print(f"    {'d_now pc':>10}{'census':>8}{'Gaia RV':>9}{'ours':>7}"
          f"{'C(Gaia)':>9}{'C(ours)':>9}")
    for t in table:
        print(f"    {t['bin_pc'][0]:4d}-{t['bin_pc'][1]:<5d}{t['census']:8d}"
              f"{t['gaia_rv']:9d}{t['ours']:7d}"
              f"{t['c_gaia_rv']:9.3f}{t['c_ours']:9.3f}")

    print(f"\n  perihelion epochs are uniform in truth, so this must be flat")
    print(f"    {'|t_ph| Myr':>12}{'N(<3pc)':>9}{'rate/Myr':>10}"
          f"{'median d_now':>14}")
    for lo in range(5):
        sub = [e for e in enc if lo <= e["t"] < lo + 1]
        near = [e for e in sub if e["d"] < 3.0]
        dn = sorted(e["dn"] for e in sub)
        med = dn[len(dn) // 2] if dn else float("nan")
        print(f"    {str(lo) + '-' + str(lo + 1):>12}{len(near):9d}"
              f"{len(near) / 2.0:10.1f}{med:14.1f}")
    print("    it is not: the deficit is stars too far away today to have an RV")

    rate, stat, rows, n_raw = rate_within(enc, C, WINDOW_MYR, FIT_CUTS)
    print(f"\n  inside the complete window |t_ph| < {WINDOW_MYR} Myr")
    print(f"    {'d_ph pc':>9}{'N':>6}{'raw/Myr':>10}{'corrected':>11}{'/d^2':>8}")
    for r in rows:
        print(f"    {r['d_ph_pc']:9.1f}{r['n']:6d}{r['raw_per_myr']:10.1f}"
              f"{r['corrected_per_myr']:11.1f}{r['per_d2']:8.2f}")

    # systematic: how far the answer moves over the fit ranges
    spread = []
    for dmax in (2.0, 3.0, 4.0, 5.0):
        cuts = tuple(c for c in FIT_CUTS if c <= dmax)
        spread.append(rate_within(enc, C, WINDOW_MYR, cuts)[0])
    sysm = (max(spread) - min(spread)) / 2.0

    kin = kinetic_check(enc)

    print(f"\n  {'=' * 60}")
    print(f"  encounter rate within 1 pc = {rate:.1f} "
          f"+/- {stat:.1f} (stat) +/- {sysm:.1f} (sys) per Myr")
    print(f"  published                  = {PUBLISHED[0]} +/- {PUBLISHED[1]} "
          f"per Myr   (Bailer-Jones et al. 2018, Gaia DR2)")
    dev = abs(rate - PUBLISHED[0]) / math.sqrt(stat ** 2 + sysm ** 2
                                               + PUBLISHED[1] ** 2)
    print(f"  agreement                  = {dev:.1f} sigma")
    print(f"  {'=' * 60}")
    print(f"\n  independent check, pi d^2 n <v> with <v> = "
          f"{kin['mean_speed_km_s']} km/s:")
    for k, v in kin["gamma_per_myr"].items():
        print(f"    {k} /pc^3 -> {v} per Myr")
    print("  (that check carries its own ~30% density uncertainty and is not "
          "used to correct anything)")

    out = {
        "window_myr": WINDOW_MYR,
        "rate_within_1pc_per_myr": round(rate, 2),
        "stat_error": round(stat, 2),
        "sys_error": round(sysm, 2),
        "n_encounters_in_fit": n_raw,
        "published": {"value": PUBLISHED[0], "error": PUBLISHED[1],
                      "source": "Bailer-Jones et al. 2018, A&A 616, A37 (Gaia DR2)"},
        "sigma_from_published": round(dev, 2),
        "completeness": table,
        "by_perihelion": rows,
        "kinetic_check": kin,
        "caveat": ("Corrected for this catalog's completeness against the "
                   "Gaia census only. Stars missing from that census - the "
                   "faintest M dwarfs, brown dwarfs, anything lost beside a "
                   "bright neighbor - are not corrected for, so this remains "
                   "a lower bound."),
    }
    (PROC / "encounter_rate.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"\n  -> {(PROC / 'encounter_rate.json').relative_to(BASE)}")

    if "--scan" in sys.argv:
        scan(enc, C)
    if "--json" in sys.argv:
        print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
