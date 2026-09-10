"""Stage 9 - how much of each prediction is the Galaxy model rather than the data.

The Monte Carlo answers "how much does this move when the measurements are
resampled". It does not answer "how much does it move if the Milky Way is not
quite the mass model we assumed", and for the candidates that dominate the
closest-approach rankings that second question turns out to be the larger one.
A star 300 pc away whose perihelion is 3 Myr out spends most of that time far
from the Sun, where the disc's pull on the two orbits differs; change the disc
mass by 15% and its closest approach moves further than its whole 68% interval.

So this integrates the entire catalog again under a set of plausible
perturbations - disc and halo mass, disc shape, bulge mass, and the solar
position and motion - and records, per star, the largest shift any of them
produces. That number is written back as `sys_dmin_pc` and shipped, so a reader
can see when a prediction is limited by the model rather than by Gaia.

It is an envelope over hand-chosen variants, not a posterior. It is deliberately
one-sided in its honesty: where it is small the prediction is genuinely
measurement-limited, and where it is large the quoted interval is too narrow.

Cheap - about 12 seconds per variant for 26,567 stars at dt = 5000 yr, which
stage 5b established agrees with a 500 yr reference to 4e-8 pc.
"""

from __future__ import annotations

import csv
import importlib
import json
import pathlib
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import astro  # noqa: E402
from astro import icrs_to_galactic_cartesian  # noqa: E402

orb = importlib.import_module("04_orbits")

PROC = pathlib.Path(__file__).parent.parent / "data" / "processed"
SYS_DT = 0.005

# Each entry perturbs the model by roughly the spread in the literature. These
# are not independent draws and combining them would be wrong; the reported
# figure is the largest single shift, i.e. an envelope.
VARIANTS = {
    "disc mass -15%":        dict(MN_M=0.85),
    "disc mass +15%":        dict(MN_M=1.15),
    "disc scale length -0.5 kpc": dict(MN_A=(3.0 - 0.5) / 3.0),
    "disc scale length +0.5 kpc": dict(MN_A=(3.0 + 0.5) / 3.0),
    "disc scale height 0.20 kpc": dict(MN_B=0.20 / 0.28),
    "disc scale height 0.40 kpc": dict(MN_B=0.40 / 0.28),
    "halo mass -20%":        dict(NFW_M=0.80),
    "halo mass +20%":        dict(NFW_M=1.20),
    "halo scale radius -20%": dict(NFW_RS=0.80),
    "bulge mass -30%":       dict(BULGE_M=0.70),
    "bulge mass +30%":       dict(BULGE_M=1.30),
}
# Solar parameters, given as absolute replacements rather than factors.
SOLAR = {
    "R0 8.091 kpc":   dict(R0_KPC=8.091),
    "R0 8.153 kpc":   dict(R0_KPC=8.153),
    "z_sun 17.8 pc":  dict(Z_SUN_PC=17.8),
    "z_sun 23.8 pc":  dict(Z_SUN_PC=23.8),
    "Vc 224 km/s":    dict(V_CIRC=224.0),
    "Vc 234 km/s":    dict(V_CIRC=234.0),
    "Vc 231.5 (the potential's own)": dict(V_CIRC=astro.circular_velocity(astro.R0_KPC)),
    "U_sun +1 km/s":  dict(U_SUN=astro.U_SUN + 1.0),
    "V_pec +1 km/s":  dict(V_SUN_PEC=astro.V_SUN_PEC + 1.0),
    "W_sun +1 km/s":  dict(W_SUN=astro.W_SUN + 1.0),
}


def col(stars, key):
    out = np.empty(len(stars))
    for i, s in enumerate(stars):
        try:
            out[i] = float(s[key])
        except (KeyError, TypeError, ValueError):
            out[i] = np.nan
    return out


# Gaia DR3 parallaxes are known to be slightly too small. Lindegren et al.
# (2021) give a magnitude-, color- and position-dependent correction Z, with
# corrected parallax = parallax - Z; Z is typically -0.017 to -0.05 mas. That
# function needs five extra Gaia columns and a published coefficient table, so it
# is not applied here. What is measured instead is how much it would matter: a
# global offset over the plausible range, which brackets the true correction for
# all but the extreme tails.
ZERO_POINTS = {
    "parallax zero-point Z = -0.017 mas": -0.017,
    "parallax zero-point Z = -0.030 mas": -0.030,
    "parallax zero-point Z = -0.050 mas": -0.050,
    "parallax zero-point Z = +0.010 mas": +0.010,
}


def perihelia(pos_pc, vel_kms) -> np.ndarray:
    """Closest approach of every star under the currently patched model."""
    sun_p, sun_v = astro.sun_phase_space()
    st = np.empty((pos_pc.shape[0] + 1, 6))
    st[0, :3], st[0, 3:] = sun_p, sun_v
    st[1:, :3] = sun_p + pos_pc / 1000.0
    st[1:, 3:] = sun_v + vel_kms
    _, _, df, _ = orb.integrate(st, +1)
    _, _, db, _ = orb.integrate(st, -1)
    return np.where(db < df, db, df)


def main() -> None:
    with open(PROC / "orbits.csv") as fh:
        stars = list(csv.DictReader(fh))
    n = len(stars)
    pos, vel = icrs_to_galactic_cartesian(
        col(stars, "ra"), col(stars, "dec"), col(stars, "parallax"),
        col(stars, "pmra"), col(stars, "pmdec"), col(stars, "rv"))

    orb.SAMPLE_EVERY = 10 ** 9
    orb.DT = SYS_DT
    baseline_defaults = {k: getattr(astro, k) for k in
                         ("MN_M", "MN_A", "MN_B", "BULGE_M", "NUC_M", "NFW_M",
                          "NFW_RS", "R0_KPC", "Z_SUN_PC", "V_CIRC",
                          "U_SUN", "V_SUN_PEC", "W_SUN")}

    print(f"  {n} stars, {len(VARIANTS) + len(SOLAR)} model variants "
          f"at dt={SYS_DT*1e6:.0f} yr")
    t0 = time.time()
    base = perihelia(pos, vel)

    worst = np.zeros(n)
    which = [""] * n
    summary = []
    for label, patch in list(VARIANTS.items()) + list(SOLAR.items()):
        scale = label in VARIANTS
        old = {k: getattr(astro, k) for k in patch}
        for k, v in patch.items():
            setattr(astro, k, baseline_defaults[k] * v if scale else v)
        d = perihelia(pos, vel)
        for k, v in old.items():
            setattr(astro, k, v)

        shift = np.abs(d - base)
        bigger = shift > worst
        worst = np.where(bigger, shift, worst)
        for i in np.flatnonzero(bigger):
            which[i] = label
        close = base < 2.0
        summary.append({
            "variant": label,
            "median_shift_pc": round(float(np.median(shift)), 6),
            "p99_shift_pc": round(float(np.percentile(shift, 99)), 6),
            "median_shift_within_2pc": round(float(np.median(shift[close])), 6),
            "max_shift_within_2pc": round(float(np.max(shift[close])), 6),
        })
        print(f"    {label:34s} median {np.median(shift):9.5f} pc   "
              f"p99 {np.percentile(shift, 99):9.4f} pc")

    # ---- parallax zero-point, kept separate from the Galaxy model -----------
    # It is an astrometric systematic, not a dynamical one, and folding the two
    # into one number would hide which is which.
    ra, dec = col(stars, "ra"), col(stars, "dec")
    plx = col(stars, "parallax")
    pmra, pmdec = col(stars, "pmra"), col(stars, "pmdec")
    rv = col(stars, "rv")
    zp_worst = np.zeros(n)
    for label, z in ZERO_POINTS.items():
        p2, v2 = icrs_to_galactic_cartesian(ra, dec, np.maximum(plx - z, 0.01),
                                            pmra, pmdec, rv)
        d = perihelia(p2, v2)
        shift = np.abs(d - base)
        zp_worst = np.maximum(zp_worst, shift)
        summary.append({
            "variant": label,
            "median_shift_pc": round(float(np.median(shift)), 6),
            "p99_shift_pc": round(float(np.percentile(shift, 99)), 6),
            "median_shift_within_2pc": round(float(np.median(shift[base < 2.0])), 6),
            "max_shift_within_2pc": round(float(np.max(shift[base < 2.0])), 6),
        })
        print(f"    {label:34s} median {np.median(shift):9.5f} pc   "
              f"p99 {np.percentile(shift, 99):9.4f} pc")

    for i, s in enumerate(stars):
        s["sys_dmin_pc"] = f"{worst[i]:.5f}"
        s["sys_dominant"] = which[i]
        s["zp_dmin_pc"] = f"{zp_worst[i]:.5f}"
    with open(PROC / "orbits.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(stars[0].keys()))
        w.writeheader()
        w.writerows(stars)

    # How often does the model matter more than the measurements?
    width = col(stars, "mc_dmin_hi") - col(stars, "mc_dmin_lo")
    med = col(stars, "mc_dmin_med")
    interesting = np.isfinite(width) & (med < 5.0) & (col(stars, "mc_censored") <= 0.5)
    dominated = interesting & (worst > width)
    now = col(stars, "dist_now_pc")
    report = {
        "note": "largest single-variant shift in the integrated closest approach; "
                "an envelope over hand-chosen model perturbations, not a posterior",
        "dt_myr": SYS_DT,
        "n_variants": len(VARIANTS) + len(SOLAR),
        "variants": summary,
        "n_within_5pc": int(interesting.sum()),
        "n_zeropoint_dominated": int((interesting & (zp_worst > width)).sum()),
        "median_zeropoint_shift_within_5pc_pc":
            round(float(np.median(zp_worst[interesting])), 5),
        "n_model_dominated": int(dominated.sum()),
        "median_shift_within_5pc_pc": round(float(np.median(worst[interesting])), 5),
        "median_shift_near_25pc_pc":
            round(float(np.median(worst[interesting & (now < 25)])), 5),
        "median_shift_beyond_100pc_pc":
            round(float(np.median(worst[interesting & (now >= 100)])), 5),
    }
    (PROC / "sensitivity.json").write_text(json.dumps(report, indent=2))

    print(f"\n  total {time.time()-t0:.0f}s")
    print(f"  of {interesting.sum()} stars approaching within 5 pc, "
          f"{dominated.sum()} have a model shift larger than their 68% interval")
    print(f"  median shift: {report['median_shift_near_25pc_pc']:.5f} pc within 25 pc, "
          f"{report['median_shift_beyond_100pc_pc']:.5f} pc beyond 100 pc")
    print(f"  parallax zero-point over its plausible range moves "
          f"{int((interesting & (zp_worst > width)).sum())} of them further than "
          f"their own interval (median {np.median(zp_worst[interesting]):.5f} pc)")
    print("  -> sys_dmin_pc, zp_dmin_pc written to data/processed/orbits.csv")
    print("  -> data/processed/sensitivity.json")


if __name__ == "__main__":
    main()
