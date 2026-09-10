"""Stage 6 - pack the catalog for the web client.

Emits three files into web/data/:
  cheb.bin    float32 Chebyshev coefficients, [n_stars][3 axes][deg+1]
  stars.json  per-star metadata as parallel arrays (compact, gzips well)
  meta.json   provenance, model parameters, validation summary

Also groups gravitationally bound multiples into systems. This matters for
honesty as much as for tidiness: catalogs list each component's instantaneous
velocity, which in a close pair includes orbital motion about the barycenter.
Alpha Centauri A and B are cataloged 7.3 km/s apart in radial velocity, which
alone moves their computed closest approach by 0.23 pc, so the components must
be presented as one system with a spread rather than as independent stars.
"""

from __future__ import annotations

import csv
import json
import pathlib
import re

import numpy as np

import sys
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from gates import require_complete  # noqa: E402
from systems import dense_ids, group_systems  # noqa: E402

BASE = pathlib.Path(__file__).parent.parent
PROC = BASE / "data" / "processed"
OUT = BASE / "web" / "data"
OUT.mkdir(parents=True, exist_ok=True)

LY_PER_PC = 3.261563777

# Effective temperature by spectral class: luminosity class V anchored on
# Pecaut & Mamajek (2013), class III on the standard giant scale. Each entry is
# the temperature at subclass 0 and at subclass 10, interpolated linearly.
#
# This is consulted only when Gaia supplies neither a temperature nor a BP-RP
# color, which is precisely the case for the bright stars that saturate Gaia
# and reach the catalog through Hipparcos. Without it they all took a flat
# 5000 K default and rendered as the same cream dot: Sirius, Vega, Arcturus,
# Altair, Procyon, Capella, Pollux, Fomalhaut and both Alpha Centauri
# components — that is, essentially every star a viewer would recognize well
# enough to check the color of. 99 of the 235 stars on the default have a
# usable spectral type; the rest keep 5000 K.
SPT_TEFF_V = {"O": (42000, 30000), "B": (30000, 9700), "A": (9700, 7200),
              "F": (7200, 5920), "G": (5920, 5280), "K": (5280, 3850),
              "M": (3850, 2400), "L": (2250, 1400), "T": (1300, 800)}
SPT_TEFF_III = {"G": (5850, 4900), "K": (4750, 3800), "M": (3800, 3200)}
SPT_RE = re.compile(r"([OBAFGKMLT])\s*(\d+(?:\.\d+)?)?\s*([IV]+)?")


def sptype_to_teff(sp: str | None) -> float | None:
    """Effective temperature from an MK spectral type, or None if unparseable.

    Luminosity class matters most for the cool end: Arcturus is K1.5 III at
    4290 K, where the same subclass on the main sequence would be near 5100 K.
    Class IV is left on the dwarf scale, which puts Procyon (F5 IV-V) at
    6560 K against a measured 6530 K.
    """
    m = SPT_RE.match((sp or "").strip())
    if not m:
        return None
    cls, sub, lum = m.group(1), m.group(2), m.group(3)
    table = SPT_TEFF_III if lum in ("I", "II", "III") and cls in SPT_TEFF_III \
        else SPT_TEFF_V
    hot, cool = table[cls]
    frac = min(float(sub), 9.9) / 10.0 if sub else 0.5
    return hot + (cool - hot) * frac


def teff_to_rgb(teff: float | None, bp_rp: float | None,
                sp_type: str | None = None) -> tuple[int, int, int]:
    """Approximate blackbody color, so O stars look blue and M dwarfs red."""
    if teff is None and bp_rp is None:
        teff = sptype_to_teff(sp_type)
    if teff is None and bp_rp is not None:
        # rough color-temperature relation for main-sequence Gaia colors
        teff = 4600.0 * (1.0 / (0.92 * bp_rp + 1.7) + 1.0 / (0.92 * bp_rp + 0.62))
    if teff is None:
        teff = 5000.0
    t = max(1000.0, min(40000.0, teff)) / 100.0
    if t <= 66:
        r = 255.0
        g = 99.4708025861 * np.log(t) - 161.1195681661
    else:
        r = 329.698727446 * (t - 60) ** -0.1332047592
        g = 288.1221695283 * (t - 60) ** -0.0755148492
    if t >= 66:
        b = 255.0
    elif t <= 19:
        b = 0.0
    else:
        b = 138.5177312231 * np.log(t - 10) - 305.0447927307
    return tuple(int(max(0, min(255, v))) for v in (r, g, b))


DISPLAY_DMIN_PC = 10.0      # what the interface calls an encounter
DISPLAY_NOW_PC = 30.0


def astrometry_provenance(rows: list[dict]) -> str:
    """Where the positions and motions came from, counted from the data.

    Written by hand this said "Gaia DR3; Hipparcos-2 for Gaia-saturated bright
    stars", which quietly dropped the one star that has neither: Scholz's star
    has DR3 photometry but no astrometric solution, being a very red L/T binary,
    so its discovery-paper astrometry is used instead. The Methods panel renders
    this string, so the omission was visible in the application. Counting the
    rows means the sentence cannot drift from the catalog again.
    """
    labels = {
        "gaiadr3": "Gaia DR3",
        "hipparcos2": "Hipparcos-2 (van Leeuwen 2007), for bright stars that "
                      "saturate Gaia and are absent from DR3",
        "simbad": "SIMBAD literature astrometry, for stars with no astrometric "
                  "solution in either catalog",
    }
    counts: dict[str, int] = {}
    for r in rows:
        counts[r.get("astrom_source") or "unknown"] = \
            counts.get(r.get("astrom_source") or "unknown", 0) + 1
    parts = [f"{n} from {labels.get(k, k)}"
             for k, n in sorted(counts.items(), key=lambda kv: -kv[1])]
    return "; ".join(parts)


def main() -> None:
    with open(PROC / "orbits.csv") as fh:
        all_stars = list(csv.DictReader(fh))
    all_coeffs = np.load(PROC / "cheb_coeffs.npy")       # (n, 3, deg+1) float32
    meta_in = json.loads((PROC / "orbit_meta.json").read_text())
    assert all_coeffs.shape[0] == len(all_stars), "coefficient/catalog mismatch"
    require_complete(all_stars, "export")

    # The catalog is screened at 20 pc so that Galactic deflection cannot push
    # a real close pass outside it. Interest is decided here, after integration,
    # on the integrated and resampled results rather than the straight-line value
    # used to find the star. The full table stays in data/star_encounters.csv.
    def num(s, k, d=1e9):
        try:
            return float(s.get(k, ""))
        except (TypeError, ValueError):
            return d

    keep = [
        i for i, s in enumerate(all_stars)
        if num(s, "dist_now_pc") < DISPLAY_NOW_PC
        or num(s, "gal_dmin_pc") < DISPLAY_DMIN_PC
        or num(s, "mc_dmin_med") < DISPLAY_DMIN_PC
        or s.get("featured")
    ]
    stars = [all_stars[i] for i in keep]
    coeffs = all_coeffs[keep]
    print(f"  shipping {len(stars)} of {len(all_stars)} stars "
          f"(within {DISPLAY_NOW_PC:.0f} pc now, or approaching within "
          f"{DISPLAY_DMIN_PC:.0f} pc, or curated)")

    def fnum(s, k, default=0.0):
        v = s.get(k, "")
        try:
            x = float(v)
        except (TypeError, ValueError):
            return default
        return default if x != x else x

    # stage 2 already grouped the systems; reuse its ids so the two stages
    # cannot drift apart
    sys_id = [int(fnum(s, "system", -1)) for s in stars]
    if all(v < 0 for v in sys_id):
        sys_id = dense_ids(group_systems(
            np.array([fnum(s, "ra") for s in stars]),
            np.array([fnum(s, "dec") for s in stars]),
            np.array([fnum(s, "parallax", 1.0) for s in stars]),
            np.array([fnum(s, "pmra") for s in stars]),
            np.array([fnum(s, "pmdec") for s in stars])))
    # Stage 2 grouped the whole 20 pc catalog, but only part of it ships. A
    # pair whose companion was not selected would otherwise arrive here as a
    # one-member "system", which is not a system at all.
    counts: dict[int, int] = {}
    for v in sys_id:
        if v >= 0:
            counts[v] = counts.get(v, 0) + 1
    remap: dict[int, int] = {}
    for k, v in enumerate(sys_id):
        if v >= 0 and counts[v] < 2:
            sys_id[k] = -1
    for k, v in enumerate(sys_id):
        if v >= 0:
            sys_id[k] = remap.setdefault(v, len(remap))
    nxt = len(remap)
    print(f"  {nxt} candidate co-moving systems covering "
          f"{sum(1 for s in sys_id if s >= 0)} shipped stars")

    # ---- per-system radial-velocity spread (the binary-contamination flag) ---
    rvs: dict[int, list[float]] = {}
    for i, sid in enumerate(sys_id):
        if sid >= 0:
            rvs.setdefault(sid, []).append(fnum(stars[i], "rv"))
    sys_rv_spread = {k: (max(v) - min(v)) for k, v in rvs.items()}
    contaminated = sum(1 for v in sys_rv_spread.values() if v > 2.0)
    print(f"  systems whose components disagree in RV by >2 km/s: {contaminated}")

    # ---- assemble ------------------------------------------------------------
    names, dist_now, dmin, tmin, grades, colors = [], [], [], [], [], []
    lin_dmin, lin_tmin, absmag, rv, flags, sptype = [], [], [], [], [], []
    mc_lo, mc_hi, mct_lo, mct_hi, mc_med, mct_med = [], [], [], [], [], []
    censored, sysd, zpd = [], [], []
    for i, s in enumerate(stars):
        names.append(s["name"])
        dist_now.append(round(fnum(s, "dist_now_pc"), 4))
        dmin.append(round(fnum(s, "gal_dmin_pc"), 4))
        tmin.append(round(fnum(s, "gal_tmin_myr"), 5))
        lin_dmin.append(round(fnum(s, "lin_dmin_pc"), 4))
        lin_tmin.append(round(fnum(s, "lin_tmin_myr"), 5))
        grades.append(s["grade"])
        rv.append(round(fnum(s, "rv"), 2))
        sptype.append(s.get("sp_type", "")[:12])
        teff = fnum(s, "teff", 0.0) or None
        bp = s.get("bp_rp", "")
        bp = float(bp) if bp not in ("", None) else None
        colors.append(teff_to_rgb(teff, bp, s.get("sp_type", "")))
        absmag.append(round(fnum(s, "abs_g", 99.0), 2))
        mc_med.append(round(fnum(s, "mc_dmin_med"), 4))
        mct_med.append(round(fnum(s, "mc_tmin_med"), 5))
        mc_lo.append(round(fnum(s, "mc_dmin_lo"), 4))
        mc_hi.append(round(fnum(s, "mc_dmin_hi"), 4))
        mct_lo.append(round(fnum(s, "mc_tmin_lo"), 5))
        mct_hi.append(round(fnum(s, "mc_tmin_hi"), 5))
        censored.append(round(fnum(s, "mc_censored"), 3))
        sysd.append(round(fnum(s, "sys_dmin_pc"), 4))
        zpd.append(round(fnum(s, "zp_dmin_pc"), 4))
        f = 0
        if s.get("is_wd") == "1":
            f |= 1
        if s.get("rv_reliability") in ("suspect_wd", "extreme", "low_snr"):
            f |= 2
        if "high_ruwe" in (s.get("astrom_flag") or ""):
            f |= 4
        if s.get("rv_source") == "simbad":
            f |= 8
        if s.get("featured"):
            f |= 16
        if sys_id[i] >= 0 and sys_rv_spread.get(sys_id[i], 0) > 2.0:
            f |= 32
        if s.get("rv_source") == "simbad_system":
            f |= 64
        # Spectroscopic binaries have no single radial velocity: the cataloged
        # number is one orbital phase. Sirius is the clearest case here - its
        # 8.17 ly / 46 kyr against the reference chart's 7.86 ly / 64 kyr
        # corresponds to an RV of -9.3 km/s rather than the cataloged -5.5.
        ot = (s.get("otype") or "")
        if ot.startswith("SB") or ot in ("EB*", "**", "SB*"):
            f |= 128
        # More than half the draws never reach a minimum inside +/-5 Myr: the
        # star is still closing at the boundary, so there is no perihelion to
        # report, only a bound.
        if fnum(s, "mc_censored") > 0.5:
            f |= 256
        # Stage 9 measured how far this prediction moves when the Galaxy model
        # is perturbed. Where that exceeds the whole resampling interval, the
        # quoted range is describing the wrong uncertainty.
        if fnum(s, "sys_dmin_pc") > (fnum(s, "mc_dmin_hi") - fnum(s, "mc_dmin_lo")):
            f |= 512
        if fnum(s, "zp_dmin_pc") > (fnum(s, "mc_dmin_hi") - fnum(s, "mc_dmin_lo")):
            f |= 1024
        flags.append(f)

    payload = {
        "names": names,
        "dist_now_pc": dist_now,
        "dmin_pc": dmin,
        "tmin_myr": tmin,
        "lin_dmin_pc": lin_dmin,
        "lin_tmin_myr": lin_tmin,
        "grade": grades,
        "rv_kms": rv,
        "abs_g": absmag,
        "sp_type": sptype,
        "color": [c for rgb in colors for c in rgb],
        "system": sys_id,
        "flags": flags,
        "mc_dmin_med": mc_med, "mc_tmin_med": mct_med,
        "mc_dmin_lo": mc_lo, "mc_dmin_hi": mc_hi,
        "mc_tmin_lo": mct_lo, "mc_tmin_hi": mct_hi,
        "mc_censored": censored,
        "sys_dmin_pc": sysd,
        "zp_dmin_pc": zpd,
    }
    (OUT / "stars.json").write_text(json.dumps(payload, separators=(",", ":")))

    coeffs.astype("<f4").tofile(OUT / "cheb.bin")

    val_path = PROC / "validation.json"
    import datetime
    raw = BASE / "data" / "raw" / "gaia_encounters.csv"
    retrieved = (datetime.date.fromtimestamp(raw.stat().st_mtime).isoformat()
                 if raw.exists() else "unknown")
    meta = {
        "catalog": "Gaia DR3",
        "retrieved": retrieved,
        "n_stars_shipped": len(stars),
        "n_stars_catalog": len(all_stars),
        "display_cut": f"within {DISPLAY_NOW_PC:.0f} pc now, or closest approach "
                       f"under {DISPLAY_DMIN_PC:.0f} pc, or curated",
        "n_stars": len(stars),
        "t_span_myr": meta_in["t_span_myr"],
        "cheb_degree": meta_in["cheb_degree"],
        "cheb_pos_err_median_pc": meta_in["cheb_pos_err_median_pc"],
        "cheb_pos_err_max_pc": meta_in["cheb_pos_err_max_pc"],
        "perihelion_err_median_pc": meta_in["perihelion_err_median_pc"],
        "perihelion_err_p99_pc": meta_in["perihelion_err_p99_pc"],
        "perihelion_err_max_pc": meta_in["perihelion_err_max_pc"],
        "perihelion_err_max_close_pc": meta_in["perihelion_err_max_close_pc"],
        "n_candidate_comoving_systems": nxt,
        "n_systems": nxt,
        "ly_per_pc": LY_PER_PC,
        "monte_carlo": {
            "draws": 256,
            "model": "every draw integrated through the Galactic potential; "
                     "distance and epoch come from the same draw population",
            "interval_meaning": "measurement-sensitivity range, not a full "
                     "confidence interval: Gaia astrometric covariance, the "
                     "parallax zero-point, and potential/solar-parameter "
                     "uncertainty are not propagated",
            # An earlier version of this note called the nominal perihelion a
            # biased estimator. That was an artifact of comparing an integrated
            # nominal orbit with straight-line draws; under one model the two
            # agree to a median 0.2%. Ranking still uses the median, because it
            # is the estimate the interval belongs to.
            "note": "mc_dmin_med/lo/hi are the 50th/16th/84th percentiles of the "
                    "closest-approach distance under resampled astrometry. Rank "
                    "on the median: it is the point estimate the quoted interval "
                    "describes, and it degrades gracefully when the parallax "
                    "error is large.",
        },
        "flag_bits": {
            "1": "white dwarf", "2": "unreliable radial velocity",
            "4": "high RUWE (probable unresolved binary)",
            "8": "radial velocity from literature, not Gaia RVS",
            "16": "curated / notable", "32": "component RVs disagree (binary contamination)",
            "64": "barycentric system velocity adopted in place of the component's",
            "128": "spectroscopic binary - cataloged RV is one orbital phase, not the barycenter",
            "512": "the Galactic mass model moves this closest approach further "
                   "than the measurement resampling does; sys_dmin_pc gives the size",
            "1024": "the Gaia parallax zero-point, over its plausible range, moves "
                    "this closest approach further than the measurement resampling "
                    "does; zp_dmin_pc gives the size",
            "256": "closest approach lies outside the +/-5 Myr window, so the quoted "
                   "distance is a bound. Which side: mc_tmin_med pinned to +5 Myr "
                   "means still closing, perihelion later; pinned to -5 Myr means "
                   "already receding, perihelion earlier; an interior mc_tmin_med "
                   "means the draws reach both ends and the side is unresolved",
        },
        "sources": {
            "astrometry": astrometry_provenance(all_stars),
            "radial_velocity": "Gaia DR3 RVS, else literature via SIMBAD",
            "potential": "Miyamoto-Nagai disc + Hernquist bulge and nucleus + NFW halo "
                         "(gala MilkyWayPotential, Price-Whelan 2017)",
            "solar_motion": "R0=8.122 kpc (GRAVITY 2018), z=20.8 pc (Bennett & Bovy 2019), "
                            "Vc=229 km/s (Eilers 2019), peculiar (11.1, 12.24, 7.25) km/s "
                            "(Schoenrich, Binney & Dehnen 2010)",
        },
    }
    # How many close candidates the Galaxy model, not Gaia, is the limit on.
    # This was quoted as a literal in the interface and again in the README,
    # and after a rebuild the two disagreed with each other and with the data
    # (168 of 1,691 and 165 of 1,684, against a true 165 of 1,697). Counted
    # here, from the same test that sets flag 512, so it cannot drift again.
    close = [s for s in all_stars
             if (fnum(s, "mc_dmin_med") or 9e9) < 5.0]
    limited = [s for s in close
               if (fnum(s, "sys_dmin_pc") or 0.0)
               > ((fnum(s, "mc_dmin_hi") or 0.0) - (fnum(s, "mc_dmin_lo") or 0.0))]
    meta["model_limited"] = {"within_pc": 5.0,
                             "candidates": len(close),
                             "limited": len(limited)}

    if val_path.exists():
        meta["validation"] = json.loads(val_path.read_text())
    # The measured encounter rate is the one derived number the project
    # produces, so the interface quotes it from here rather than from a figure
    # typed into the page, which could then drift from the data behind it.
    rate_path = PROC / "encounter_rate.json"
    if rate_path.exists():
        r = json.loads(rate_path.read_text())
        meta["encounter_rate"] = {
            k: r[k] for k in ("rate_within_1pc_per_myr", "stat_error",
                              "sys_error", "window_myr", "published",
                              "sigma_from_published", "caveat")
        }
    (OUT / "meta.json").write_text(json.dumps(meta, indent=1))

    kb = lambda p: p.stat().st_size / 1024
    print(f"  -> web/data/cheb.bin    {kb(OUT/'cheb.bin'):8.0f} KB")
    print(f"  -> web/data/stars.json  {kb(OUT/'stars.json'):8.0f} KB")
    print(f"  -> web/data/meta.json   {kb(OUT/'meta.json'):8.0f} KB")


if __name__ == "__main__":
    main()
