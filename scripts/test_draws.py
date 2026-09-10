"""Tests for the Monte Carlo's astrometric draws.

`correlated_draws` is the numerical core of the uncertainty stage: get it wrong
and every interval in the catalog is wrong in a way that still looks entirely
plausible. So the properties are asserted rather than eyeballed - that zero
correlation reproduces exactly the independent draws the previous version made,
that a prescribed correlation comes back out of the sample, that an inconsistent
correlation triple degrades to independence instead of raising, and that the
zero-crossing parallax fix does not pile draws up on the mean.

    python3 scripts/test_draws.py
"""

from __future__ import annotations

import importlib
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
unc = importlib.import_module("055_uncertainty")

N = 200_000
TOL_SD, TOL_R = 0.02, 0.01

PLX, PLX_E = np.array([10.0, 50.0]), np.array([0.1, 0.5])
PMA, PMA_E = np.array([100.0, -30.0]), np.array([0.2, 0.4])
PMD, PMD_E = np.array([-40.0, 10.0]), np.array([0.3, 0.15])


def draws(rpa, rpd, rad, n=N, rng=None):
    z = np.full(2, 0.0)
    return unc.correlated_draws(
        PLX, PLX_E, PMA, PMA_E, PMD, PMD_E,
        z + rpa, z + rpd, z + rad, n, rng or np.random.default_rng(1))


def report(label, ok, detail):
    print(f"    [{'PASS' if ok else 'FAIL'}] {label}\n           {detail}")
    return ok


def main() -> None:
    print("  Monte Carlo draw tests:")
    passed = True

    d = draws(0.0, 0.0, 0.0)
    sd = np.stack([d[:, :, k].std(0) for k in range(3)])
    want = np.stack([PLX_E, PMA_E, PMD_E])
    c = np.corrcoef(d[:, 0, :].T)
    off = max(abs(c[0, 1]), abs(c[0, 2]), abs(c[1, 2]))
    passed &= report(
        "zero correlation reproduces independent draws",
        np.allclose(sd, want, rtol=TOL_SD) and off < TOL_R,
        f"sd within {TOL_SD:.0%}, largest spurious correlation {off:.4f}")

    means = d.mean(0)[0]
    passed &= report(
        "draws are centered on the published values",
        np.allclose(means, [PLX[0], PMA[0], PMD[0]], atol=0.02),
        f"{means.round(4)} vs {[PLX[0], PMA[0], PMD[0]]}")

    worst = 0.0
    for rpa, rpd, rad in [(0.5, -0.3, 0.2), (0.9, 0.1, 0.05), (-0.7, 0.6, -0.4)]:
        c = np.corrcoef(draws(rpa, rpd, rad)[:, 1, :].T)
        worst = max(worst, abs(c[0, 1] - rpa), abs(c[0, 2] - rpd), abs(c[1, 2] - rad))
    passed &= report(
        "prescribed correlations come back out of the sample",
        worst < TOL_R, f"largest error over three triples: {worst:.4f}")

    d = draws(0.99, 0.99, -0.99, n=1000)          # not positive definite
    passed &= report(
        "an inconsistent correlation triple degrades instead of raising",
        np.isfinite(d).all() and d.shape == (1000, 2, 3),
        f"shape {d.shape}, all finite, fell back to the diagonal")

    rng = np.random.default_rng(2)
    mean, sd1 = np.array([0.5]), np.array([1.0])
    raw = rng.normal(mean, sd1, (N, 1))
    fixed = unc.resample_nonpositive(raw.copy(), mean, sd1, rng)
    at_mean = float((np.abs(fixed - mean[0]) < 1e-12).mean())
    passed &= report(
        "zero-crossing parallaxes are resampled, not parked on the mean",
        (fixed > 0).all() and at_mean < 1e-4,
        f"{(raw <= 0.01).mean():.1%} of raw draws were non-positive; "
        f"{at_mean:.6f} of the result sits exactly on the mean")

    print(f"\n  {'all draw tests pass' if passed else 'DRAW TESTS FAILED'}")
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
