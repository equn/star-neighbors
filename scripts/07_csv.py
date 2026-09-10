"""Stage 7 - write the flat catalog deliverable.

This existed only as an ad-hoc snippet, which meant the published CSV could not
be regenerated from the repository. Column meanings live here and are written
out alongside it so the file documents itself.
"""

from __future__ import annotations

import csv
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from gates import require_complete  # noqa: E402

BASE = pathlib.Path(__file__).parent.parent
PROC = BASE / "data" / "processed"
OUT = BASE / "data"

COLUMNS = [
    ("name", "common name"),
    ("simbad_id", "SIMBAD identifier"),
    ("source_id", "Gaia DR3 source_id"),
    ("hip", "Hipparcos number"),
    ("sp_type", "spectral type"),
    ("system", "candidate co-moving system id, -1 if none"),
    ("dist_now_pc", "present distance, pc"),
    ("ra", "RA deg (ICRS)"), ("dec", "Dec deg (ICRS)"),
    ("parallax", "parallax mas"), ("parallax_error", "parallax error mas"),
    ("parallax_catalog", "the Gaia/Hipparcos parallax, where a better published "
                           "value replaced it"),
    ("parallax_source", "source of the replacement, blank if the catalog value stands"),
    ("pmra", "pmRA* mas/yr"), ("pmdec", "pmDec mas/yr"),
    ("pmra_error", "pmRA* error mas/yr"), ("pmdec_error", "pmDec error mas/yr"),
    ("rv", "radial velocity used, km/s"),
    ("rv_error", "RV error km/s, as published. The Monte Carlo substitutes "
                 "2.0 where this is missing or non-positive and floors it at "
                 "0.1; this column is the untouched catalog value, so the "
                 "substitution is reconstructible"),
    ("rv_source", "gaia | simbad | simbad_system (barycentric)"),
    ("rv_component", "component RV, where a barycentric value replaced it"),
    ("system_parent", "SIMBAD parent object supplying the barycentric RV"),
    ("x", "heliocentric Galactic X pc (toward the Galactic center)"),
    ("y", "Y pc (toward l=90)"), ("z", "Z pc (toward the North Galactic Pole)"),
    ("vx", "U km/s"), ("vy", "V km/s"), ("vz", "W km/s"),
    ("gal_dmin_pc", "closest approach of the NOMINAL orbit, integrated, pc"),
    ("gal_tmin_myr", "epoch of that closest approach, Myr"),
    ("lin_dmin_pc", "closest approach, straight-line model, pc"),
    ("lin_tmin_myr", "epoch, straight-line model, Myr"),
    ("mc_dmin_med", "closest approach, median of integrated draws, pc  [PREFER THIS]"),
    ("mc_dmin_lo", "16th percentile of integrated draws, pc"),
    ("mc_dmin_hi", "84th percentile of integrated draws, pc"),
    ("mc_tmin_med", "epoch, median of the SAME draws, Myr  [PREFER THIS]"),
    ("mc_tmin_lo", "16th percentile epoch, Myr"),
    ("mc_tmin_hi", "84th percentile epoch, Myr"),
    ("mc_censored", "fraction of draws whose minimum separation is clipped to "
                    "EITHER end of the +/-5 Myr window. Above 0.5 there is no "
                    "perihelion inside the window and mc_dmin_med is an upper "
                    "bound, not an encounter. Which side: mc_tmin_med at +5 Myr "
                    "means the star is still closing and the encounter is later; "
                    "at -5 Myr it was already receding and the encounter is "
                    "earlier; an interior mc_tmin_med with mc_censored > 0.5 "
                    "means the draws reach both ends and the side is unresolved"),
    ("mc_model", "uncertainty model used ('integrated')"),
    ("sys_dmin_pc", "largest shift in the closest approach produced by any one "
                    "of 21 perturbations of the Galactic mass model and the solar "
                    "position/motion. Compare with (mc_dmin_hi - mc_dmin_lo): where "
                    "it is larger, the prediction is limited by the Galaxy model "
                    "rather than by Gaia. An envelope over chosen variants, not a "
                    "posterior"),
    ("sys_dominant", "which perturbation produced that largest shift"),
    ("zp_dmin_pc", "largest shift produced by a Gaia parallax zero-point in the "
                   "range Z = -0.05 to +0.01 mas. The Lindegren et al. (2021) "
                   "correction is NOT applied to the catalog; this is how much "
                   "applying it could move each result"),
    ("grade", "data quality A/B/C - about the measurements, not the prediction"),
    ("rv_reliability", "RV reliability flag"),
    ("astrom_flag", "astrometry flags"),
    ("is_wd", "1 if white dwarf"),
    ("featured", "curated label if notable"),
]

HEADER_NOTE = """# Column reference for `star_encounters.csv`

Sorted by the median closest approach of the integrated Monte Carlo draws
(nearest first). Distances in parsec (1 pc = 3.2616 light years); times in Myr
from today, negative = past.

Closest approaches are bounded by the +/-5 Myr integration window. Where
`mc_censored` exceeds 0.5 the encounter lies outside the window and
`mc_dmin_med` is an upper bound, not an encounter. That happens in three ways,
distinguished by `mc_tmin_med`: pinned to +5 Myr the star is still closing and
its perihelion is later; pinned to -5 Myr it was already receding when the
window opened and its perihelion is earlier; and for a handful of stars the
median sits inside the window because the draws reach both ends, meaning the
astrometry cannot say which side the encounter is on.

`sys_dmin_pc` and `zp_dmin_pc` are separate, non-statistical figures: how far the closest
approach moves when the Galactic mass model or the solar parameters are
perturbed within the range the literature allows, and how far it moves under a
Gaia parallax zero-point across its plausible range. Neither is part of the
`mc_*` interval and neither should be added to it in quadrature.

`mc_*` columns come from resampling the published astrometry within its stated
errors and integrating every draw through the Galactic potential. They are
measurement-sensitivity ranges, not full confidence intervals: Gaia's
astrometric covariance and the parallax zero-point are not propagated, and
uncertainty in the potential and the solar parameters is reported separately in
`sys_dmin_pc` rather than folded in. See README.md.

"""


def main() -> None:
    with open(PROC / "orbits.csv") as fh:
        rows = list(csv.DictReader(fh))

    def key(r):
        for k in ("mc_dmin_med", "gal_dmin_pc"):
            try:
                v = float(r.get(k, ""))
                if v == v:
                    return v
            except (TypeError, ValueError):
                pass
        return 1e9

    # A published CSV with silently blank science columns is worse than no CSV:
    # the header still promises them and nothing downstream can tell. The same
    # gate stage 6 uses, so the two exporters cannot disagree about what counts
    # as a finished run.
    require_complete(rows, "catalog")

    names = [c for c, _ in COLUMNS]
    missing = [c for c in names if c not in rows[0]]
    if missing:
        raise SystemExit(f"refusing to write the catalog: orbits.csv has no "
                         f"{missing}. Run the missing pipeline stage first.")

    out = OUT / "star_encounters.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(names)
        for r in sorted(rows, key=key):
            w.writerow([r.get(c, "") for c in names])

    (OUT / "COLUMNS.md").write_text(
        HEADER_NOTE + "\n".join(f"- `{c}` — {d}" for c, d in COLUMNS) + "\n")
    print(f"  -> {out.relative_to(BASE)}  ({out.stat().st_size/1024:.0f} KB, {len(rows)} rows)")
    print(f"  -> {(OUT / 'COLUMNS.md').relative_to(BASE)}")


if __name__ == "__main__":
    main()
