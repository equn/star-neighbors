"""Stage 1 - acquire raw astrometry.

Three complementary sources, because no single catalog covers the problem:

  A. Gaia DR3 "encounter candidates" - every star with 6D astrometry whose
     straight-line closest approach to the Sun is < 20 pc within +/-6 Myr.
     The screen is deliberately twice the distance anyone would call an
     encounter: Galactic deflection moves a multi-Myr perihelion by of order a
     parsec, and in one measured case by 4 pc, so screening at the same 10 pc
     used for scientific interest would let real close passes fall outside the
     candidate list entirely.
     The perihelion solution is evaluated server-side in ADQL so the 33M-star
     radial-velocity catalog is filtered by ESA rather than downloaded.
  B. Gaia DR3 "present-day neighborhood" - everything within 30 pc, including
     stars with no radial velocity (they still belong on the map of today).
  C. Bright stars Gaia cannot measure. Sirius, Alpha/Beta Centauri et al.
     saturate Gaia's detectors, so their DR3 astrometry is absent or unusable.
     Those come from Hipparcos-2 (van Leeuwen 2007) + SIMBAD radial velocities.

Outputs raw CSVs into data/raw/.
"""

from __future__ import annotations

import csv
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from tap import SIMBAD_TAP, async_query, sync_query  # noqa: E402

RAW = pathlib.Path(__file__).parent.parent / "data" / "raw"
RAW.mkdir(parents=True, exist_ok=True)

# 4.740470446 km/s = 1 AU/yr; converts mas/yr at a given parallax into km/s.
KAPPA = 4.740470446
# 1 pc / (km/s) expressed in Myr.
PC_KMS_TO_MYR = 0.977792

GAIA_COLS = """source_id, ra, dec, parallax, parallax_error, pmra, pmra_error,
    pmdec, pmdec_error, radial_velocity, radial_velocity_error,
    phot_g_mean_mag, bp_rp, teff_gspphot, ruwe,
    rv_expected_sig_to_noise, rv_nb_transits, rv_template_teff,
    astrometric_excess_noise, non_single_star,
    parallax_pmra_corr, parallax_pmdec_corr, pmra_pmdec_corr"""


def write(name: str, rows: list[dict]) -> None:
    path = RAW / name
    if not rows:
        print(f"  !! {name}: no rows")
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"  -> {path.name}: {len(rows)} rows")


def fetch_encounters() -> None:
    """A. Stars whose straight-line closest approach is < 20 pc within +/-6 Myr."""
    print("[A] Gaia DR3 encounter candidates (server-side perihelion filter)")
    q = f"""
    SELECT {GAIA_COLS}
    FROM (
      SELECT {GAIA_COLS},
             1000.0/parallax AS dist_pc,
             {KAPPA} * SQRT(pmra*pmra + pmdec*pmdec) / parallax AS vt
      FROM gaiadr3.gaia_source
      WHERE parallax > 1
        AND radial_velocity IS NOT NULL
        AND parallax_over_error > 5
    ) AS s
    WHERE s.dist_pc * s.vt / SQRT(s.vt*s.vt + s.radial_velocity*s.radial_velocity) < 20
      AND ABS({PC_KMS_TO_MYR} * s.dist_pc * s.radial_velocity
              / (s.vt*s.vt + s.radial_velocity*s.radial_velocity)) < 6
    """
    write("gaia_encounters.csv", async_query(q))


def fetch_neighborhood() -> None:
    """B. Everything Gaia sees within 30 pc, radial velocity or not."""
    print("[B] Gaia DR3 present-day neighborhood (< 30 pc)")
    q = f"""
    SELECT {GAIA_COLS}
    FROM gaiadr3.gaia_source
    WHERE parallax > 33.333333
      AND parallax_over_error > 10
    """
    write("gaia_30pc.csv", async_query(q))


def fetch_hipparcos_bright() -> None:
    """C. Bright stars saturated in Gaia, from the Hipparcos-2 reduction."""
    print("[C] Hipparcos-2 bright stars (Gaia-saturated)")
    q = """
    SELECT hip, ra, dec, plx, e_plx, pm_ra AS pmra, pm_de AS pmdec,
           e_pm_ra, e_pm_de, hp_mag, b_v
    FROM public.hipparcos_newreduction
    WHERE plx > 5 AND hp_mag < 6.5 AND plx/e_plx > 5
    """
    write("hipparcos_bright.csv", async_query(q))


def fetch_simbad_rv() -> None:
    """C2. SIMBAD radial velocities + identities for bright nearby stars."""
    print("[C2] SIMBAD radial velocities for bright nearby stars")
    q = """
    SELECT b.main_id, b.ra, b.dec, b.plx_value, b.plx_err,
           b.pmra, b.pmdec, b.rvz_radvel, b.rvz_err, b.rvz_bibcode,
           b.sp_type, b.otype, id.id AS hip_id
    FROM basic AS b
    JOIN ident AS id ON id.oidref = b.oid
    WHERE b.plx_value > 5
      AND b.rvz_radvel IS NOT NULL
      AND b.rvz_type = 'v'
      AND id.id LIKE 'HIP %'
    """
    write("simbad_bright.csv", async_query(q, base=SIMBAD_TAP))


def fetch_simbad_gaia_rv() -> None:
    """D. SIMBAD radial velocities keyed by Gaia DR3 source_id.

    Most nearby M dwarfs have no Gaia RV (too faint/red for RVS) but do have a
    published RV. Without this, the actual nearest stars drop out of the
    time-domain analysis entirely.
    """
    print("[D] SIMBAD radial velocities keyed to Gaia DR3 ids (< 50 pc)")
    q = """
    SELECT b.main_id, b.plx_value, b.rvz_radvel, b.rvz_err, b.rvz_bibcode,
           b.sp_type, b.otype, id.id AS gaia_id
    FROM basic AS b
    JOIN ident AS id ON id.oidref = b.oid
    WHERE b.plx_value > 20
      AND b.rvz_radvel IS NOT NULL
      AND b.rvz_type = 'v'
      AND id.id LIKE 'Gaia DR3 %'
    """
    write("simbad_gaia_rv.csv", async_query(q, base=SIMBAD_TAP))


def fetch_featured_astrometry() -> None:
    """E. Full SIMBAD astrometry for the curated featured stars.

    A few important objects have no usable Gaia solution: Scholz's star, the
    closest known past flyby, is a very red L/T binary with DR3 photometry but
    no astrometry at all. SIMBAD carries the discovery-paper values.
    """
    print("[E] SIMBAD astrometry for featured stars")
    import featured as _f
    quoted = ",".join("'" + n.replace("'", "''") + "'" for n in _f.FEATURED)
    q = f"""
    SELECT i.id AS query_id, b.main_id, b.ra, b.dec, b.plx_value, b.plx_err,
           b.pmra, b.pmdec, b.rvz_radvel, b.rvz_err, b.rvz_type, b.rvz_bibcode,
           b.sp_type, b.otype
    FROM ident AS i JOIN basic AS b ON b.oid = i.oidref
    WHERE i.id IN ({quoted})
    """
    write("featured_astrometry.csv", sync_query(q, base=SIMBAD_TAP))


if __name__ == "__main__":
    which = sys.argv[1:] or ["a", "b", "c", "c2", "d", "e"]
    if "a" in which:
        fetch_encounters()
    if "b" in which:
        fetch_neighborhood()
    if "c" in which:
        fetch_hipparcos_bright()
    if "c2" in which:
        fetch_simbad_rv()
    if "d" in which:
        fetch_simbad_gaia_rv()
    if "e" in which:
        fetch_featured_astrometry()
