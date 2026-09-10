"""Stage 5b - Monte Carlo perihelion uncertainty, integrated through the potential.

A closest-approach distance quoted without an error bar is close to meaningless
here. The catalog contains stars presently hundreds of parsecs away whose
solution brings them within a parsec of the Sun; at that range a 5% parallax
error moves the whole trajectory by tens of parsecs. The formal quality grade
does not capture this, because the astrometry itself can be excellent while the
extrapolation is not.

Every draw is integrated through the same Milky Way potential as the nominal
orbit, and distance *and* epoch are both taken from that one draw population.

An earlier version solved the draws with straight-line motion, on the grounds
that stage 4 measured the potential's effect at under 0.005 pc for encounters a
few tens of kyr away. That reasoning was sound for those stars and wrong for
everything else: the same stage measured a median shift of 0.53 pc at 2.5-4.5
Myr, and the closest-approach rankings are dominated by multi-Myr candidates.
Integrating the draws moves HD 7977 from 0.150 to 0.042 pc, TYC 5584-552-1 from
0.901 to 1.627 pc and 2MASS J05250205+0137210 from 5.659 to 1.582 pc, while
leaving Proxima Centauri, Ross 248 and Scholz's star unchanged to three
decimals. Mixing a straight-line distance with an integrated epoch, as that
version did, was not one coherent estimate of anything.

Gaia solves parallax and both proper-motion components in one fit, so their
errors are correlated - the median |r| across this catalog is about 0.15 and
reaches 0.91. Drawing them independently samples an axis-aligned box instead of
the real error ellipsoid, which for a fifth of the close candidates changes the
interval width by more than 10%. The draws below use the published correlations
where Gaia supplies them, and fall back to independence for the Hipparcos and
SIMBAD rows, which publish none.

What is still missing, and why the intervals are called measurement-sensitivity
ranges rather than confidence intervals:
  * The Gaia parallax zero-point correction is not applied. Stage 9 measures how
    much that matters rather than leaving it as a bare caveat.
  * Radial-velocity errors are assumed where a catalog does not publish one.
  * Uncertainty in the Galactic potential and in the solar parameters is not
    folded in here; stage 9 reports it separately as sys_dmin_pc, and for a
    tenth of the close candidates it is the larger term.
"""

from __future__ import annotations

import csv
import importlib
import pathlib
import sys
import time

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
orb = importlib.import_module("04_orbits")
from astro import icrs_to_galactic_cartesian, sun_phase_space  # noqa: E402

PROC = pathlib.Path(__file__).parent.parent / "data" / "processed"

N_DRAWS = 256
STARS_PER_BATCH = 200
# Integrating the draws at 5000 yr agrees with a 500 yr reference to 4e-8 pc,
# including for an 885 km/s star reaching 0.022 pc, because the perihelion is
# refined between steps rather than sampled. That is far below every other error
# here and halves a ~2 hour run.
MC_DT = 0.005
RV_ERR_FALLBACK = 2.0
RNG = np.random.default_rng(20260730)


def col(stars, key, default=np.nan):
    out = np.empty(len(stars))
    for i, s in enumerate(stars):
        v = s.get(key, "")
        try:
            x = float(v)
        except (TypeError, ValueError):
            x = default
        out[i] = default if x != x else x
    return out


def correlated_draws(plx, plx_e, pmra, pmra_e, pmdec, pmdec_e,
                     r_pa, r_pd, r_ad, n, rng):
    """Draw (parallax, pmRA*, pmDec) from Gaia's actual error ellipsoid.

    Returns (n, m, 3). Built as a Cholesky factor per star rather than by
    calling multivariate_normal m times, which would be minutes rather than
    milliseconds at this catalog size. Where a correlation is missing - the
    Hipparcos and SIMBAD rows - it is zero, and the result is exactly the
    independent draw those stars had before.

    A published correlation triple need not be a positive-definite matrix once
    rounded, so a failed factorisation falls back to the diagonal rather than
    producing nonsense.
    """
    m = len(plx)
    sd = np.stack([plx_e, pmra_e, pmdec_e], axis=-1)              # (m, 3)
    corr = np.zeros((m, 3, 3))
    corr[:, 0, 0] = corr[:, 1, 1] = corr[:, 2, 2] = 1.0
    corr[:, 0, 1] = corr[:, 1, 0] = r_pa
    corr[:, 0, 2] = corr[:, 2, 0] = r_pd
    corr[:, 1, 2] = corr[:, 2, 1] = r_ad
    cov = corr * sd[:, :, None] * sd[:, None, :]
    try:
        L = np.linalg.cholesky(cov)
    except np.linalg.LinAlgError:
        L = np.zeros_like(cov)
        for i in range(m):
            try:
                L[i] = np.linalg.cholesky(cov[i])
            except np.linalg.LinAlgError:
                L[i] = np.diag(sd[i])                 # treat as independent
    z = rng.standard_normal((n, m, 3))
    mean = np.stack([plx, pmra, pmdec], axis=-1)
    return mean[None, :, :] + np.einsum('mij,nmj->nmi', L, z)


def resample_nonpositive(plx_draws, mean, sigma, rng):
    """Redraw parallaxes that crossed zero, keeping the rest untouched.

    Substituting the mean for a bad draw quietly narrows the distribution and
    biases it toward the nominal value. Resampling keeps the draws honest;
    anything still non-positive after a few passes is clipped to a small
    positive parallax, which puts the star far away rather than behind the
    observer. The correlated companions of a redrawn parallax are left alone -
    the bias this removes is far larger than the correlation it breaks, and
    only 1 draw in ~10^4 is affected.
    """
    out = plx_draws
    for _ in range(6):
        bad = out <= 0.01
        if not bad.any():
            break
        out = np.where(bad, rng.normal(np.broadcast_to(mean, out.shape),
                                       np.broadcast_to(sigma, out.shape)), out)
    return np.maximum(out, 0.01)


def truncated_normal_parallax(mean, sigma, n):
    """Draw parallaxes, resampling the tail that crosses zero.

    Substituting the mean for a bad draw (the previous behavior) quietly
    narrows the distribution and biases it toward the nominal value. Resampling
    keeps the draws honest; anything still non-positive after a few passes is
    clipped to a small positive parallax, which puts the star far away rather
    than behind the observer.
    """
    out = RNG.normal(mean, sigma, (n, len(mean)))
    for _ in range(6):
        bad = out <= 0.01
        if not bad.any():
            break
        out[bad] = RNG.normal(
            np.broadcast_to(mean, out.shape)[bad],
            np.broadcast_to(sigma, out.shape)[bad],
        )
    return np.maximum(out, 0.01)


def main() -> None:
    with open(PROC / "orbits.csv") as fh:
        stars = list(csv.DictReader(fh))
    n = len(stars)

    ra, dec = col(stars, "ra"), col(stars, "dec")
    plx, plx_e = col(stars, "parallax"), col(stars, "parallax_error")
    pmra, pmdec = col(stars, "pmra"), col(stars, "pmdec")
    pmra_e, pmdec_e = col(stars, "pmra_error"), col(stars, "pmdec_error")
    rv, rv_e = col(stars, "rv"), col(stars, "rv_error")
    # Gaia's astrometric correlations; 0 where a catalog publishes none.
    r_pa = col(stars, "plx_pmra_corr", 0.0)
    r_pd = col(stars, "plx_pmdec_corr", 0.0)
    r_ad = col(stars, "pmra_pmdec_corr", 0.0)
    for a_ in (r_pa, r_pd, r_ad):
        np.nan_to_num(a_, copy=False)
    n_corr = int(((r_pa != 0) | (r_pd != 0) | (r_ad != 0)).sum())

    plx_e = np.where(np.isnan(plx_e) | (plx_e <= 0), 0.02 * np.abs(plx), plx_e)
    pmra_e = np.where(np.isnan(pmra_e) | (pmra_e <= 0), plx_e, pmra_e)
    pmdec_e = np.where(np.isnan(pmdec_e) | (pmdec_e <= 0), plx_e, pmdec_e)
    # A published error of 1e-4 km/s on a literature radial velocity is not
    # credible and would produce an absurdly tight interval; floor it.
    rv_e = np.where(np.isnan(rv_e) | (rv_e <= 0), RV_ERR_FALLBACK, rv_e)
    n_floored = int((rv_e < 0.1).sum())
    rv_e = np.maximum(rv_e, 0.1)

    sun_p, sun_v = sun_phase_space()
    orb.SAMPLE_EVERY = 10 ** 9          # trajectories are not needed here
    orb.DT = MC_DT

    keys = ("dmin_med", "dmin_lo", "dmin_hi", "tmin_med", "tmin_lo", "tmin_hi",
            "censored")
    res = {k: np.zeros(n) for k in keys}

    print(f"  {N_DRAWS} draws x {n} stars, integrated at dt={MC_DT*1e6:.0f} yr")
    print(f"  floored {n_floored} implausibly small radial-velocity errors to 0.1 km/s")
    print(f"  drawing from the published astrometric error ellipsoid for {n_corr} "
          f"stars; the remaining {n - n_corr} publish no correlations")
    t0 = time.time()
    for a in range(0, n, STARS_PER_BATCH):
        b = min(a + STARS_PER_BATCH, n)
        m = b - a
        sh = (N_DRAWS, m)

        astro_s = correlated_draws(
            plx[a:b], plx_e[a:b], pmra[a:b], pmra_e[a:b], pmdec[a:b], pmdec_e[a:b],
            r_pa[a:b], r_pd[a:b], r_ad[a:b], N_DRAWS, RNG)
        p_s = resample_nonpositive(astro_s[:, :, 0], plx[a:b], plx_e[a:b], RNG)
        a_s = astro_s[:, :, 1]
        d_s = astro_s[:, :, 2]
        r_s = RNG.normal(rv[a:b], rv_e[a:b], sh)
        ra_s = np.broadcast_to(ra[a:b], sh)
        dec_s = np.broadcast_to(dec[a:b], sh)

        pos, vel = icrs_to_galactic_cartesian(
            ra_s.ravel(), dec_s.ravel(), p_s.ravel(),
            a_s.ravel(), d_s.ravel(), r_s.ravel())

        k = pos.shape[0]
        st = np.empty((k + 1, 6))
        st[0, :3], st[0, 3:] = sun_p, sun_v
        st[1:, :3] = sun_p + pos / 1000.0
        st[1:, 3:] = sun_v + vel

        _, _, df, tf = orb.integrate(st, +1)
        _, _, db, tb = orb.integrate(st, -1)
        back = db < df
        dmin = np.where(back, db, df).reshape(sh)
        tmin = np.where(back, tb, tf).reshape(sh)

        res["dmin_lo"][a:b], res["dmin_med"][a:b], res["dmin_hi"][a:b] = \
            np.percentile(dmin, [16, 50, 84], axis=0)
        res["tmin_lo"][a:b], res["tmin_med"][a:b], res["tmin_hi"][a:b] = \
            np.percentile(tmin, [16, 50, 84], axis=0)
        # a draw whose minimum sits on the window edge has not had an encounter
        res["censored"][a:b] = np.mean(np.abs(tmin) >= orb.T_SPAN - 1e-6, axis=0)

        if (a // STARS_PER_BATCH) % 10 == 0:
            el = time.time() - t0
            done = b / n
            print(f"    {b}/{n}  {el/60:.1f} min elapsed"
                  f"{f', ~{el/done*(1-done)/60:.0f} min left' if done > 0.02 else ''}")

    for i, s in enumerate(stars):
        s["mc_dmin_med"] = f"{res['dmin_med'][i]:.5f}"
        s["mc_dmin_lo"] = f"{res['dmin_lo'][i]:.5f}"
        s["mc_dmin_hi"] = f"{res['dmin_hi'][i]:.5f}"
        s["mc_tmin_med"] = f"{res['tmin_med'][i]:.5f}"
        s["mc_tmin_lo"] = f"{res['tmin_lo'][i]:.5f}"
        s["mc_tmin_hi"] = f"{res['tmin_hi'][i]:.5f}"
        s["mc_censored"] = f"{res['censored'][i]:.3f}"
        s["mc_model"] = "integrated"

    with open(PROC / "orbits.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(stars[0].keys()))
        w.writeheader()
        w.writerows(stars)

    gal = col(stars, "gal_dmin_pc")
    width = res["dmin_hi"] - res["dmin_lo"]
    close = gal < 2.0
    dist_now = col(stars, "dist_now_pc")
    print(f"\n  total {(time.time()-t0)/60:.1f} min")
    print(f"  stars with integrated perihelion < 2 pc: {close.sum()}")
    for lo, hi, lab in ((0, 25, 'now < 25 pc'), (25, 100, '25-100 pc'), (100, 1e9, '> 100 pc')):
        k = close & (dist_now >= lo) & (dist_now < hi)
        if k.sum():
            print(f"    {lab:12s} N={k.sum():4d}  median 68% width = {np.median(width[k]):7.3f} pc"
                  f"   median predicted d_min = {np.median(gal[k]):6.3f} pc")
    print(f"  trajectories censored at the window edge in >50% of draws: "
          f"{int((res['censored'] > 0.5).sum())}")
    print("  -> integrated uncertainty columns written to data/processed/orbits.csv")


if __name__ == "__main__":
    main()
