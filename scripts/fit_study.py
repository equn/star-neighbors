"""Diagnose the Chebyshev compression error and pick a fit that is honest.

The number that matters is not the positional residual at the sampled times but
the error in the closest-approach distance recovered from the compressed
trajectory, because that is the quantity the whole project reports.
"""

from __future__ import annotations

import csv
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import importlib
orb = importlib.import_module("04_orbits")
from astro import sun_phase_space  # noqa: E402

PROC = pathlib.Path(__file__).parent.parent / "data" / "processed"


def main() -> None:
    with open(PROC / "named.csv") as fh:
        stars = list(csv.DictReader(fh))
    n = len(stars)
    pos_h = np.array([[float(s["x"]), float(s["y"]), float(s["z"])] for s in stars])
    vel_h = np.array([[float(s["vx"]), float(s["vy"]), float(s["vz"])] for s in stars])
    sp, sv = sun_phase_space()
    st = np.empty((n + 1, 6))
    st[0, :3], st[0, 3:] = sp, sv
    st[1:, :3] = sp + pos_h / 1000.0
    st[1:, 3:] = sv + vel_h

    print("integrating (dense sampling for the fit study) ...")
    orb.SAMPLE_EVERY = 5
    tf, off_f, dmin_f, tmin_f = orb.integrate(st, +1)
    tb, off_b, dmin_b, tmin_b = orb.integrate(st, -1)
    times = np.concatenate([tb[::-1][:-1], tf])
    offs = np.concatenate([off_b[::-1][:-1], off_f], axis=0)
    truth = np.where(dmin_b < dmin_f, dmin_b, dmin_f)      # from full-step tracking
    print(f"  {len(times)} samples per star")

    u = times / orb.T_SPAN
    # dense grid to locate the minimum of each fitted trajectory
    ug = np.linspace(-1, 1, 8001)

    dn = np.array([float(s["dist_now_pc"]) for s in stars])

    for deg in (12, 16, 20, 24):
        coef = np.empty((n, 3, deg + 1))
        for ax in range(3):
            coef[:, ax, :] = np.polynomial.chebyshev.chebfit(u, offs[:, :, ax], deg).T

        for dtype, tag in ((np.float64, "f64"), (np.float32, "f32")):
            c = coef.astype(dtype).astype(np.float64)
            B = np.polynomial.chebyshev.chebvander(ug, deg)
            got = np.empty(n)
            step = 2000
            for a in range(0, n, step):
                b = min(a + step, n)
                p = np.einsum("nk,ijk->nij", B, c[a:b])
                got[a:b] = np.linalg.norm(p, axis=2).min(axis=0)
            err = np.abs(got - truth)
            close = truth < 2.0
            print(f"  deg {deg:2d} {tag}: p50 {np.median(err):.2e}  p99 {np.percentile(err,99):.2e}  "
                  f"max {err.max():.2e} pc | d_min<2pc: max {err[close].max():.2e} pc "
                  f"| payload {n*3*(deg+1)*4/1048576:.2f} MB")

    # is the error concentrated in distant, fast stars?
    coef12 = np.empty((n, 3, 13))
    for ax in range(3):
        coef12[:, ax, :] = np.polynomial.chebyshev.chebfit(u, offs[:, :, ax], 12).T
    B = np.polynomial.chebyshev.chebvander(ug, 12)
    got = np.empty(n)
    for a in range(0, n, 2000):
        b = min(a + 2000, n)
        got[a:b] = np.linalg.norm(np.einsum("nk,ijk->nij", B, coef12[a:b]), axis=2).min(axis=0)
    err = np.abs(got - truth)
    print("\n  degree 12 error vs present distance:")
    for lo, hi in ((0, 30), (30, 100), (100, 200), (200, 1e9)):
        m = (dn >= lo) & (dn < hi)
        if m.sum():
            print(f"    now {lo:>4}-{hi if hi<1e9 else '  +':>4} pc  N={m.sum():5d}  "
                  f"median {np.median(err[m]):.2e}  max {err[m].max():.2e} pc")


if __name__ == "__main__":
    main()
