"""Straight-line vs integrated Monte Carlo, both computed here.

This script exists to document why stage 5b integrates every draw. The catalog
used to resample astrometry and solve each draw with straight-line motion, which
was justified for encounters a few tens of kyr away - stage 4 measures the
potential's effect there at under 0.005 pc - but not for the "closest ahead"
list, which is dominated by candidates 1-5 Myr out, exactly where stage 4
measures a median shift of 0.53 pc.

It used to read `mc_dmin_med` out of orbits.csv and label it "straight-line".
That was true while the shipped column was straight-line and became quietly
false the moment stage 5b changed, at which point it was comparing an integrated
result with another integrated result and reporting the difference as a model
comparison. Both columns below are now computed here, from the same draws, under
the same resampling rules as stage 5b (zero-crossing parallaxes resampled rather
than replaced, radial-velocity errors floored at 0.1 km/s).
"""

from __future__ import annotations

import csv
import importlib
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))

orb = importlib.import_module("04_orbits")
unc = importlib.import_module("055_uncertainty")
from astro import (icrs_to_galactic_cartesian, linear_perihelion,  # noqa: E402
                   sun_phase_space)

PROC = pathlib.Path(__file__).parent.parent / "data" / "processed"
N_DRAWS = unc.N_DRAWS
RNG = np.random.default_rng(20260730)

TARGETS = [
    "Gliese 710", "HD 7977", "TYC 5584-552-1",
    "Gaia DR3 4591145358813373312", "Gaia DR3 6913732624445112832",
    "2MASS J05250205+0137210",
    "Proxima Centauri", "Ross 248", "LSPM J2146+3813", "Scholz's Star",
]


def draws_for(s: dict, n: int):
    """Resample one star's astrometry exactly as stage 5b does."""
    f = lambda k, d=0.0: float(s[k]) if s.get(k) not in ("", None) else d
    plx, plx_e = f("parallax"), f("parallax_error")
    plx_e = plx_e if plx_e > 0 else 0.02 * abs(plx)
    pmra_e = f("pmra_error") or plx_e
    pmdec_e = f("pmdec_error") or plx_e
    rv_e = max(f("rv_error") or unc.RV_ERR_FALLBACK, 0.1)

    p = unc.truncated_normal_parallax(np.array([plx]), np.array([plx_e]), n).ravel()
    return (np.full(n, f("ra")), np.full(n, f("dec")), p,
            RNG.normal(f("pmra"), pmra_e, n),
            RNG.normal(f("pmdec"), pmdec_e, n),
            RNG.normal(f("rv"), rv_e, n))


def main() -> None:
    stars = {r["name"]: r for r in csv.DictReader(open(PROC / "orbits.csv"))}
    sun_p, sun_v = sun_phase_space()
    orb.SAMPLE_EVERY = 10 ** 9          # no trajectory samples needed
    orb.DT = unc.MC_DT

    print(f"  {N_DRAWS} draws per star; both models computed from the same draws\n")
    print(f"  {'star':30s} {'straight-line':>14s} {'integrated':>13s} {'shift':>9s}"
          f"  {'integrated 68%':>16s}")
    for name in TARGETS:
        s = stars.get(name)
        if s is None:
            print(f"  {name:30s}   not in catalog")
            continue
        ra, dec, plx, pmra, pmdec, rv = draws_for(s, N_DRAWS)
        pos, vel = icrs_to_galactic_cartesian(ra, dec, plx, pmra, pmdec, rv)

        # straight-line: the star and the Sun drift at constant velocity, so the
        # perihelion is closed-form in the relative phase space
        d_lin, _ = linear_perihelion(pos, vel)

        st = np.empty((N_DRAWS + 1, 6))
        st[0, :3], st[0, 3:] = sun_p, sun_v
        st[1:, :3] = sun_p + pos / 1000.0
        st[1:, 3:] = sun_v + vel

        _, _, df, tf = orb.integrate(st, +1)
        _, _, db, tb = orb.integrate(st, -1)
        d_int = np.where(db < df, db, df)

        lin = float(np.median(d_lin))
        lo, med, hi = np.percentile(d_int, [16, 50, 84])
        print(f"  {name[:30]:30s} {lin:11.3f} pc {med:10.3f} pc {med - lin:+8.3f}  "
              f"[{lo:6.3f},{hi:7.3f}]")

    print("\n  The shipped catalog uses the integrated column.")


if __name__ == "__main__":
    main()
