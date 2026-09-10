"""Stage 4 - integrate every star and the Sun through the Galactic potential.

Why not straight lines
----------------------
Both reference charts assume the Sun and each star drift in straight lines. Over
the +/-80 kyr they cover that is fine. Over 5 Myr it is not: the Sun completes
~2% of a Galactic orbit, the disc's vertical restoring force bends trajectories
back toward the midplane on a ~70 Myr period, and shear from differential
rotation separates stars that a linear model keeps together. So both the Sun and
every star are integrated as test particles in a Milky Way potential and the
separation is measured between the integrated orbits.

Compression
-----------
Storing sampled trajectories for 12k stars would be tens of megabytes. Instead
each star's heliocentric offset is fitted with a Chebyshev polynomial in time.
The paths are smooth, so a low-order fit reproduces them to a small fraction of
a parsec, the whole catalog collapses to ~1 MB, and a polynomial is cheap
enough to evaluate per-star per-frame on the GPU.
"""

from __future__ import annotations

import csv
import json
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from astro import acceleration, sun_phase_space  # noqa: E402

BASE = pathlib.Path(__file__).parent.parent
PROC = BASE / "data" / "processed"

# 1 km/s expressed in kpc/Myr. The same factor converts (km/s)^2/kpc into
# (km/s)/Myr, so it serves for both halves of the state derivative.
KMS_TO_KPC_PER_MYR = 1.0227121650537077e-3

T_SPAN = 5.0            # Myr, integrated both directions
DT = 0.002              # Myr per step (2000 yr); orbital period is ~230 Myr
SAMPLE_EVERY = 10       # store a sample every 20 kyr for the polynomial fit
CHEB_DEG = 12


def derivative(state: np.ndarray) -> np.ndarray:
    """d/dt of (pos_kpc, vel_kms) stacked as (N, 6), time in Myr."""
    out = np.empty_like(state)
    out[:, :3] = state[:, 3:] * KMS_TO_KPC_PER_MYR
    out[:, 3:] = acceleration(state[:, :3]) * KMS_TO_KPC_PER_MYR
    return out


def rk4(state: np.ndarray, dt: float) -> np.ndarray:
    k1 = derivative(state)
    k2 = derivative(state + 0.5 * dt * k1)
    k3 = derivative(state + 0.5 * dt * k2)
    k4 = derivative(state + dt * k3)
    return state + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def integrate(state0: np.ndarray, direction: int):
    """Integrate one direction. Returns (times, samples, dmin, tmin).

    samples is (n_samples, N, 3) heliocentric offsets in pc, star minus Sun,
    where index 0 of the body axis is the Sun itself.

    The closest approach is refined between steps rather than taken as the
    smallest sampled separation. That matters more than it sounds: the
    catalog contains stars moving at hundreds of km/s, and one hypervelocity
    star covers 1.65 pc per 2000-year step, so simply picking the smallest
    sampled distance overestimated its perihelion by 0.2 pc (1.58 against a
    converged 1.38). Over one step the relative motion is effectively straight,
    which makes the squared separation an exact parabola in time, so fitting a
    parabola through the three samples bracketing the minimum recovers the true
    perihelion almost exactly - and far more cheaply than shrinking the step.
    """
    n_steps = int(round(T_SPAN / DT))
    state = state0.copy()
    n_bodies = state.shape[0]

    times = [0.0]
    offs = [(state[1:, :3] - state[0, :3]) * 1000.0]

    d2 = np.sum(offs[0] ** 2, axis=1)
    best2 = d2.copy()                      # raw sampled minimum, as a floor
    tmin = np.zeros(n_bodies - 1)
    prev2, prev1 = None, d2               # squared separation, two steps back

    for step in range(1, n_steps + 1):
        state = rk4(state, direction * DT)
        t = direction * step * DT
        rel = (state[1:, :3] - state[0, :3]) * 1000.0
        d2 = np.sum(rel ** 2, axis=1)

        closer = d2 < best2
        best2 = np.where(closer, d2, best2)
        tmin = np.where(closer, t, tmin)

        if prev2 is not None:
            # was the previous step a local minimum? if so, interpolate it
            is_min = (prev1 <= prev2) & (prev1 <= d2)
            denom = prev2 - 2.0 * prev1 + d2
            ok = is_min & (denom > 0)
            delta = np.where(ok, 0.5 * (prev2 - d2) / np.where(ok, denom, 1.0), 0.0)
            refined = prev1 - 0.25 * (prev2 - d2) * delta
            t_ref = (t - direction * DT) + delta * direction * DT
            upd = ok & (refined >= 0) & (refined < best2)
            best2 = np.where(upd, refined, best2)
            tmin = np.where(upd, t_ref, tmin)

        prev2, prev1 = prev1, d2

        if step % SAMPLE_EVERY == 0:
            times.append(t)
            offs.append(rel)

    return (np.array(times), np.array(offs, dtype=np.float64),
            np.sqrt(np.maximum(best2, 0.0)), tmin)


def main() -> None:
    src = PROC / "named.csv"
    if not src.exists():
        src = PROC / "merged.csv"
    with open(src) as fh:
        stars = list(csv.DictReader(fh))
    n = len(stars)
    print(f"  {n} stars from {src.name}")

    # heliocentric Galactic Cartesian (pc, km/s) -> Galactocentric (kpc, km/s)
    pos_helio = np.array([[float(s["x"]), float(s["y"]), float(s["z"])] for s in stars])
    vel_helio = np.array([[float(s["vx"]), float(s["vy"]), float(s["vz"])] for s in stars])
    sun_pos, sun_vel = sun_phase_space()

    state0 = np.empty((n + 1, 6))
    state0[0, :3] = sun_pos
    state0[0, 3:] = sun_vel
    state0[1:, :3] = sun_pos + pos_helio / 1000.0
    state0[1:, 3:] = sun_vel + vel_helio

    print(f"  integrating +/-{T_SPAN} Myr at dt={DT*1e6:.0f} yr "
          f"({int(T_SPAN/DT)} steps each way, {n+1} bodies) ...")
    t_f, off_f, dmin_f, tmin_f = integrate(state0, +1)
    print("    forward done")
    t_b, off_b, dmin_b, tmin_b = integrate(state0, -1)
    print("    backward done")

    # stitch the two halves into one increasing time axis
    times = np.concatenate([t_b[::-1][:-1], t_f])
    offs = np.concatenate([off_b[::-1][:-1], off_f], axis=0)

    better_back = dmin_b < dmin_f
    dmin = np.where(better_back, dmin_b, dmin_f)
    tmin = np.where(better_back, tmin_b, tmin_f)

    # ---- Chebyshev compression ---------------------------------------------
    # (dmin/tmin above are the refined perihelia from the integration)
    # Fit each coordinate as a function of time over the whole window.
    u = times / T_SPAN                       # map to [-1, 1]
    coeffs = np.empty((n, 3, CHEB_DEG + 1), dtype=np.float64)
    for axis in range(3):
        c = np.polynomial.chebyshev.chebfit(u, offs[:, :, axis], CHEB_DEG)
        coeffs[:, axis, :] = c.T

    # honest error measurement against the integrated samples
    recon = np.empty_like(offs)
    for axis in range(3):
        # chebval with 2-D coefficients returns (n_series, n_points)
        recon[:, :, axis] = np.polynomial.chebyshev.chebval(
            u, coeffs[:, axis, :].T
        ).T
    err = np.linalg.norm(recon - offs, axis=2)
    print(f"  Chebyshev deg {CHEB_DEG}: median err {np.median(err):.3e} pc, "
          f"p99 {np.percentile(err, 99):.3e} pc, max {err.max():.3e} pc")

    # The residual at the sampled times understates what matters. The quantity
    # this project reports is the closest approach, so measure the error in that
    # directly: minimize the stored polynomial on a fine grid and compare with
    # the refined perihelion from the integration.
    # The grid minimum must be refined the same way the orbit's was, or the
    # measurement is limited by grid spacing rather than by the fit: a fast star
    # crosses a parsec between grid points and the comparison then reports the
    # sampling error of the yardstick.
    ug = np.linspace(-1.0, 1.0, 8001)
    Bg = np.polynomial.chebyshev.chebvander(ug, CHEB_DEG)
    c32 = coeffs.astype(np.float32).astype(np.float64)      # exactly what ships
    fit_dmin = np.empty(n)
    for a in range(0, n, 1500):
        b = min(a + 1500, n)
        d2 = np.sum(np.einsum("nk,ijk->nij", Bg, c32[a:b]) ** 2, axis=2)
        k = d2.argmin(axis=0)
        idx = np.arange(b - a)
        f0 = d2[k, idx]
        kk = np.clip(k, 1, len(ug) - 2)
        fm, fp = d2[kk - 1, idx], d2[kk + 1, idx]
        den = fm - 2 * f0 + fp
        dl = np.where(den > 0, 0.5 * (fm - fp) / np.where(den > 0, den, 1.0), 0.0)
        # trust the parabola only when its vertex lies inside the bracket
        ok = (k > 0) & (k < len(ug) - 1) & (den > 0) & (np.abs(dl) <= 0.5)
        ref = np.where(ok, f0 - 0.25 * (fm - fp) * dl, f0)
        fit_dmin[a:b] = np.sqrt(np.maximum(np.minimum(ref, f0), 0.0))
    perr = np.abs(fit_dmin - dmin)
    close = dmin < 2.0
    print(f"  perihelion error of the shipped float32 fit: median {np.median(perr):.2e} pc, "
          f"p99 {np.percentile(perr, 99):.2e} pc, max {perr.max():.2e} pc")
    if close.any():
        print(f"    among the {close.sum()} stars closing to within 2 pc: "
              f"max {perr[close].max():.2e} pc")

    # ---- outputs ------------------------------------------------------------
    np.save(PROC / "cheb_coeffs.npy", coeffs.astype(np.float32))

    lin_d = np.array([float(s["lin_dmin_pc"]) for s in stars])
    lin_t = np.array([float(s["lin_tmin_myr"]) for s in stars])
    for i, s in enumerate(stars):
        s["gal_dmin_pc"] = f"{dmin[i]:.6f}"
        s["gal_tmin_myr"] = f"{tmin[i]:.6f}"
        s["dmin_shift_pc"] = f"{dmin[i] - lin_d[i]:.6f}"
        s["tmin_shift_myr"] = f"{tmin[i] - lin_t[i]:.6f}"
    with open(PROC / "orbits.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(stars[0].keys()))
        w.writeheader()
        w.writerows(stars)

    # how much did the Galactic potential actually matter?
    interesting = lin_d < 5.0
    if interesting.any():
        shift = np.abs(dmin - lin_d)[interesting]
        tsh = np.abs(tmin - lin_t)[interesting]
        print(f"  linear vs integrated, stars with d_lin < 5 pc (N={interesting.sum()}):")
        print(f"    |delta d_min|: median {np.median(shift):.4f} pc, p90 {np.percentile(shift, 90):.4f} pc, max {shift.max():.4f} pc")
        print(f"    |delta t_min|: median {np.median(tsh)*1000:.1f} kyr, max {tsh.max()*1000:.1f} kyr")

    meta = {
        "n_stars": n, "t_span_myr": T_SPAN, "dt_myr": DT,
        "cheb_degree": CHEB_DEG, "n_samples": len(times),
        "cheb_pos_err_median_pc": float(np.median(err)),
        "cheb_pos_err_max_pc": float(err.max()),
        "perihelion_err_median_pc": float(np.median(perr)),
        "perihelion_err_p99_pc": float(np.percentile(perr, 99)),
        "perihelion_err_max_pc": float(perr.max()),
        "perihelion_err_max_close_pc": float(perr[dmin < 2.0].max()) if (dmin < 2.0).any() else 0.0,
    }
    (PROC / "orbit_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"  -> cheb_coeffs.npy {coeffs.astype(np.float32).nbytes/1024:.0f} KB (float32), orbits.csv")


if __name__ == "__main__":
    main()
