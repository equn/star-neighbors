"""Stage 2 - merge the catalogs into one 6D phase-space table.

Resolution rules
----------------
* Gaia DR3 is the default source of astrometry.
* A star's radial velocity comes from Gaia RVS if measured, otherwise from the
  literature via SIMBAD. Nearby M dwarfs are usually too faint/red for RVS, so
  without the SIMBAD fallback the genuinely nearest stars would be missing.
* Stars bright enough to saturate Gaia (Sirius, Alpha/Beta Centauri, ...) are
  absent from DR3 entirely; those come from Hipparcos-2 + a SIMBAD RV, matched
  against Gaia both by catalog cross-match and by position so nothing is
  duplicated.

Every star also gets a data-quality grade. This matters more than it sounds:
a spurious parallax or radial velocity turns an ordinary distant star into a
spectacular fake close encounter, and those fakes are precisely what rise to the
top of any "closest approach" ranking.

Kept if within 30 pc today, or straight-line closest approach under 20 pc within
+/-5.5 Myr, or on the curated featured list. That 20 pc is a deliberately
generous screen, not a claim of interest - see the keep rule below.
"""

from __future__ import annotations

import csv
import math
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from astro import icrs_to_galactic_cartesian, linear_perihelion  # noqa: E402
from systems import (RV_DISAGREEMENT_KMS, dense_ids, group_systems,  # noqa: E402
                     rv_spread)
from tap import SIMBAD_TAP, async_query, sync_query  # noqa: E402

BASE = pathlib.Path(__file__).parent.parent
RAW = BASE / "data" / "raw"
PROC = BASE / "data" / "processed"
PROC.mkdir(parents=True, exist_ok=True)

HIP_EPOCH, GAIA_EPOCH = 1991.25, 2016.0

# Literature radial velocities are scraped from heterogeneous papers and some are
# simply wrong: SIMBAD lists EZ Aquarii (an M dwarf 3.4 pc away) at 6824 km/s and
# HD 163318 at -830 km/s. Left unchecked each becomes a spectacular fake close
# encounter, because a huge inward RV is exactly what the perihelion filter
# rewards. Genuine hypervelocity stars are vanishingly rare locally, so bound the
# fallback. Gaia's own RVs are kept but graded.
RV_SANITY_KMS = 500.0

# Systems where a published parallax is demonstrably better than the catalog
# value this pipeline would otherwise use. Deliberately tiny and deliberately
# sourced - this is not a place to tune results.
#
# Alpha Centauri A and B saturate Gaia and are absent from DR3, so they arrive
# from Hipparcos-2, which gives them parallaxes that disagree with each other by
# 5.6%: 754.81 and 796.92 mas. That put the two components 0.228 ly apart in
# present distance, for a pair that orbits at 23 AU - a physical separation of
# 0.00036 ly, so the catalog was wrong by a factor of 630 and the error showed
# up directly in their closest approaches (2.88 vs 3.13 ly). Akeson et al. 2021
# measured the system to 750.81 +/- 0.38 mas with millimeter astrometry
# (doi:10.3847/1538-3881/abfaff), which is both far more precise and, being one
# measurement of one system, self-consistent between the components.
PUBLISHED_PARALLAX = {
    "71683": (750.81, 0.38, "Akeson et al. 2021, mm astrometry of alpha Cen AB"),
    "71681": (750.81, 0.38, "Akeson et al. 2021, mm astrometry of alpha Cen AB"),
}


def load(name: str) -> list[dict]:
    with open(RAW / name) as fh:
        return list(csv.DictReader(fh))


def num(row: dict, key: str):
    v = row.get(key)
    if v in (None, "", "NaN", "nan", "null"):
        return None
    try:
        x = float(v)
    except ValueError:
        return None
    return None if x != x else x


def fetch_xmatch() -> dict[str, str]:
    """HIP number -> Gaia DR3 source_id, from the official cross-match."""
    path = RAW / "hip_gaia_xmatch.csv"
    if not path.exists():
        print("  fetching Hipparcos<->Gaia cross-match ...")
        rows = async_query(
            "SELECT x.source_id, x.original_ext_source_id AS hip "
            "FROM gaiadr3.hipparcos2_best_neighbor AS x"
        )
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["source_id", "hip"])
            w.writeheader()
            w.writerows(rows)
    return {r["hip"].strip(): r["source_id"] for r in load("hip_gaia_xmatch.csv")}


def _ident_pairs(pairs: list[tuple[str, str]]) -> set[tuple[str, str]]:
    """Ask SIMBAD which (HIP, Gaia source_id) pairs are one object.

    Used only for the handful of positional matches the parallax agreement test
    rejects. Relaxing that test globally would merge genuinely different stars
    that happen to sit close together; asking the identifier database resolves
    the specific cases on evidence instead. rho Ori (HIP 24331 = Gaia DR3
    3235349837027282560) is the reason this exists: Hipparcos gives 9.32 +/- 0.94
    mas and Gaia 13.329 +/- 0.031 mas, a 30.1% disagreement, so the star arrived
    in the catalog twice, once at 75 pc and once at 107 pc.

    The answer is cached, so a rebuild needs no network.
    """
    path = RAW / "hip_gaia_ident.csv"
    known: set[tuple[str, str]] = set()
    if path.exists():
        known = {(r["hip"], r["source_id"]) for r in load("hip_gaia_ident.csv")}
    todo = [p for p in pairs if p not in known]
    if todo:
        ids = ",".join(f"'Gaia DR3 {sid}'" for _, sid in todo)
        rows = sync_query(
            "SELECT i1.id AS gaia_id, i2.id AS hip_id FROM ident AS i1 "
            "JOIN ident AS i2 ON i1.oidref = i2.oidref "
            f"WHERE i1.id IN ({ids}) AND i2.id LIKE 'HIP %'",
            base=SIMBAD_TAP)
        same = {(r["hip_id"].replace("HIP", "").strip(),
                 r["gaia_id"].replace("Gaia DR3", "").strip()) for r in rows}
        confirmed = [p for p in todo if p in same]
        known |= set(confirmed)
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["hip", "source_id"])
            w.writeheader()
            w.writerows({"hip": h, "source_id": s} for h, s in sorted(known))
        print(f"    SIMBAD identifiers resolved {len(confirmed)} of "
              f"{len(todo)} ambiguous positional matches")
    return known


def _positional_dedup(gaia_stars: list[dict], hip_rows: list[dict]) -> set[str]:
    """HIP numbers already present in the Gaia list, matched on sky + parallax.

    gaiadr3.hipparcos2_best_neighbor covers only ~99.5k of ~118k Hipparcos
    entries and the gaps cluster among bright stars, so this backs it up.
    Hipparcos positions are epoch 1991.25 and Gaia's 2016.0; a high-proper-motion
    star such as eps Eri moves ~24 arcsec between them, so epochs must be
    aligned first or every nearby star looks like a non-match.

    Position alone is not enough - two stars can sit 10 arcsec apart - so the
    parallaxes must also agree. Where they do not, the pair goes to SIMBAD's
    identifier tables rather than being assumed distinct.
    """
    from scipy.spatial import cKDTree

    def unit(ra_deg, dec_deg):
        ra, dec = np.radians(ra_deg), np.radians(dec_deg)
        return np.stack([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra),
                         np.sin(dec)], axis=-1)

    g_ra = np.array([s["ra"] for s in gaia_stars])
    g_dec = np.array([s["dec"] for s in gaia_stars])
    g_plx = np.array([s["parallax"] for s in gaia_stars])
    tree = cKDTree(unit(g_ra, g_dec))

    dt = GAIA_EPOCH - HIP_EPOCH
    matched: set[str] = set()
    ambiguous: list[tuple[str, str]] = []
    tol_rad = np.radians(10.0 / 3600.0)
    for r in hip_rows:
        ra, dec = num(r, "ra"), num(r, "dec")
        pmra, pmdec, plx = num(r, "pmra"), num(r, "pmdec"), num(r, "plx")
        if None in (ra, dec, pmra, pmdec, plx) or plx <= 0:
            continue
        hip = r["hip"].strip()
        dec2 = dec + pmdec * dt / 3.6e6
        ra2 = ra + (pmra * dt / 3.6e6) / max(np.cos(np.radians(dec)), 1e-6)
        near = tree.query_ball_point(unit(ra2, dec2), r=2 * np.sin(tol_rad / 2))
        if not near:
            continue
        if any(abs(g_plx[i] - plx) / plx < 0.30 for i in near):
            matched.add(hip)
            continue
        # position agrees, parallax does not: let the identifier tables decide
        best = min(near, key=lambda i: (g_ra[i] - ra2) ** 2 + (g_dec[i] - dec2) ** 2)
        ambiguous.append((hip, gaia_stars[best]["source_id"]))

    if ambiguous:
        print(f"  {len(ambiguous)} positional matches with disagreeing parallaxes; "
              f"checking SIMBAD identifiers")
        same = _ident_pairs(ambiguous)
        for pair in ambiguous:
            if pair in same:
                matched.add(pair[0])
    return matched


def grade(star: dict) -> tuple[str, str, str]:
    """Return (rv_reliability, astrometry_flag, overall grade A/B/C)."""
    plx, plx_err = star["parallax"], star["parallax_error"]
    gmag, bp_rp = star["gmag"], star["bp_rp"]
    ruwe, nss = star["ruwe"], star["nss"]

    # White dwarfs: absolute magnitude far below the main sequence for their color.
    m_g = None
    if gmag is not None and plx and plx > 0:
        m_g = gmag + 5 * math.log10(plx) - 10
    is_wd = bool(
        (m_g is not None and bp_rp is not None and m_g > 10 and bp_rp < 1.5)
        or "WD" in (star.get("otype") or "")
    )

    # Gaia's RVS cross-correlates against normal-star templates. On a degenerate
    # atmosphere that fails outright: within 50 pc only 3 white dwarfs have any
    # DR3 radial velocity and their median |RV| is 374 km/s at S/N ~ 6, against
    # 17.7 km/s at S/N ~ 367 for comparable main-sequence stars. Treat as junk.
    rv_snr = star.get("rv_snr")
    if is_wd and star["rv_source"] == "gaia":
        rv_rel = "suspect_wd"
    elif abs(star["rv"]) > 400:
        rv_rel = "extreme"       # locally almost always an artifact, not a real runaway
    elif rv_snr is not None and rv_snr < 5:
        rv_rel = "low_snr"
    elif star["rv_source"] == "simbad_system":
        # A barycentric velocity is the right quantity, but it is still a
        # literature value and often an old one (Alpha Cen's is from 1979), so
        # it should not outrank a Gaia RVS measurement.
        rv_rel = "literature_system"
    elif star["rv_source"] == "simbad":
        rv_rel = "literature"
    else:
        rv_rel = "good"

    astro = ""
    if ruwe is not None and ruwe > 1.4:
        astro = "high_ruwe"          # usually an unresolved binary
    if nss:
        astro = "non_single" if not astro else astro + ",non_single"
    if plx and plx_err and plx / plx_err < 10:
        astro = "weak_parallax" if not astro else astro + ",weak_parallax"

    if rv_rel in ("suspect_wd", "low_snr", "extreme") or "weak_parallax" in astro:
        g = "C"
    elif rv_rel in ("literature", "literature_system") or astro:
        g = "B"
    else:
        g = "A"
    return rv_rel, astro, g, is_wd, m_g


def apply_system_rv(stars: list[dict]) -> int:
    """Give bound components their system's barycentric radial velocity.

    Only systems whose components actually disagree are touched, and only where
    SIMBAD's parent object carries a real (type 'v') velocity. Everything else
    is left exactly as measured.
    """
    ra = np.array([s["ra"] for s in stars])
    dec = np.array([s["dec"] for s in stars])
    plx = np.array([s["parallax"] for s in stars])
    pmra = np.array([s["pmra"] for s in stars])
    pmdec = np.array([s["pmdec"] for s in stars])
    sys_id = dense_ids(group_systems(ra, dec, plx, pmra, pmdec))
    for s, sid in zip(stars, sys_id):
        s["system"] = sid

    spread = rv_spread(sys_id, [s["rv"] for s in stars])
    bad = {k for k, v in spread.items() if v > RV_DISAGREEMENT_KMS}
    members = [i for i, sid in enumerate(sys_id) if sid in bad]
    print(f"  {len(bad)} systems disagree internally by >{RV_DISAGREEMENT_KMS} km/s "
          f"({len(members)} components); looking up barycentric velocities")
    if not members:
        return 0

    # child identifier -> star index
    key_of: dict[str, int] = {}
    for i in members:
        st = stars[i]
        if st["source_id"]:
            key_of[f"Gaia DR3 {st['source_id']}"] = i
        elif st["hip"]:
            key_of[f"HIP {st['hip']}"] = i
    keys = list(key_of)

    # SIMBAD's h_link table gives each component's parent object
    parent_rv: dict[int, tuple] = {}
    CH = 300
    for a in range(0, len(keys), CH):
        part = keys[a:a + CH]
        quoted = ",".join("'" + k.replace("'", "''") + "'" for k in part)
        q = f"""
        SELECT c.id AS child_id, p.main_id AS parent_id, p.otype AS parent_otype,
               p.rvz_radvel, p.rvz_err, p.rvz_type
        FROM h_link AS h
        JOIN ident AS c ON c.oidref = h.child
        JOIN basic AS p ON p.oid = h.parent
        WHERE c.id IN ({quoted})
        """
        try:
            rows = sync_query(q, base=SIMBAD_TAP)
        except Exception as exc:                       # keep the pipeline usable
            print(f"    !! SIMBAD h_link lookup failed: {exc}")
            return 0
        for r in rows:
            idx = key_of.get(" ".join(r["child_id"].split()))
            rv = num(r, "rvz_radvel")
            if idx is None or rv is None:
                continue
            if (r.get("rvz_type") or "").strip() != "v" or abs(rv) > RV_SANITY_KMS:
                continue
            parent_rv[idx] = (rv, num(r, "rvz_err"), r["parent_id"].strip())

    n_fixed = 0
    for i, got in parent_rv.items():
        rv, rv_err, parent = got
        st = stars[i]
        st["rv_component"] = st["rv"]
        st["rv"] = rv
        st["rv_error"] = rv_err if rv_err else st["rv_error"]
        st["rv_source"] = "simbad_system"
        st["system_parent"] = parent
        n_fixed += 1
    print(f"  adopted a barycentric RV for {n_fixed} components")
    return n_fixed


def main() -> None:
    # ---- featured stars ----------------------------------------------------
    feat = load("featured_resolved.csv")
    feat_gaia = {r["gaia_source_id"]: r["label"] for r in feat if r["gaia_source_id"]}
    feat_hip = {r["hip"]: r["label"] for r in feat if r["hip"]}

    # ---- radial-velocity lookups -------------------------------------------
    rv_by_gaia: dict[str, tuple] = {}
    simbad_by_gaia: dict[str, dict] = {}
    n_rej = 0
    for r in load("simbad_gaia_rv.csv"):
        sid = r["gaia_id"].replace("Gaia DR3", "").strip()
        simbad_by_gaia[sid] = r
        rv = num(r, "rvz_radvel")
        if rv is None:
            continue
        if abs(rv) > RV_SANITY_KMS:
            n_rej += 1
            continue
        rv_by_gaia[sid] = (rv, num(r, "rvz_err"))

    rv_by_hip: dict[str, tuple] = {}
    for r in load("simbad_bright.csv"):
        hip = r["hip_id"].replace("HIP", "").strip()
        rv = num(r, "rvz_radvel")
        if rv is None:
            continue
        if abs(rv) > RV_SANITY_KMS:
            n_rej += 1
            continue
        rv_by_hip[hip] = (rv, num(r, "rvz_err"), r["main_id"].strip(),
                          r.get("sp_type", ""), r.get("otype", ""))
    print(f"  literature RVs rejected as unphysical (>{RV_SANITY_KMS:.0f} km/s): {n_rej}")

    # ---- Gaia stars ---------------------------------------------------------
    gaia: dict[str, dict] = {}
    for fname in ("gaia_encounters.csv", "gaia_30pc.csv"):
        for r in load(fname):
            gaia.setdefault(r["source_id"], r)

    stars: list[dict] = []
    n_gaia_rv = n_simbad_rv = n_norv = 0
    for sid, r in gaia.items():
        plx = num(r, "parallax")
        if plx is None or plx <= 0:
            continue
        rv, rv_err, rv_src = num(r, "radial_velocity"), num(r, "radial_velocity_error"), "gaia"
        if rv is None:
            got = rv_by_gaia.get(sid)
            if not got:
                n_norv += 1
                continue
            rv, rv_err, rv_src = got[0], got[1], "simbad"
            n_simbad_rv += 1
        else:
            n_gaia_rv += 1
        sb = simbad_by_gaia.get(sid, {})
        stars.append({
            "id": f"gaia{sid}", "source_id": sid, "hip": "",
            "ra": num(r, "ra"), "dec": num(r, "dec"),
            "parallax": plx, "parallax_error": num(r, "parallax_error"),
            "pmra": num(r, "pmra"), "pmdec": num(r, "pmdec"),
            "pmra_error": num(r, "pmra_error"), "pmdec_error": num(r, "pmdec_error"),
            "rv": rv, "rv_error": rv_err, "rv_source": rv_src,
            "rv_snr": num(r, "rv_expected_sig_to_noise"),
            # Gaia solves parallax and both proper-motion components together,
            # so their errors are correlated - median |r| ~ 0.15 across this
            # catalog and up to 0.91. Stage 5b needs these to draw from the
            # actual error ellipsoid rather than an axis-aligned box.
            "plx_pmra_corr": num(r, "parallax_pmra_corr"),
            "plx_pmdec_corr": num(r, "parallax_pmdec_corr"),
            "pmra_pmdec_corr": num(r, "pmra_pmdec_corr"),
            "gmag": num(r, "phot_g_mean_mag"), "bp_rp": num(r, "bp_rp"),
            "teff": num(r, "teff_gspphot"), "ruwe": num(r, "ruwe"),
            "nss": int(num(r, "non_single_star") or 0),
            "astrom_source": "gaiadr3",
            "simbad_id": sb.get("main_id", "").strip(),
            "sp_type": sb.get("sp_type", "").strip(),
            "otype": sb.get("otype", "").strip(),
            "featured": feat_gaia.get(sid, ""),
        })
    print(f"  Gaia stars usable: {len(stars)}  "
          f"(RV: {n_gaia_rv} Gaia + {n_simbad_rv} SIMBAD; dropped {n_norv} with no RV)")

    # ---- Hipparcos bright stars Gaia cannot see -----------------------------
    hip_to_gaia = fetch_xmatch()
    have_gaia = {s["source_id"] for s in stars}
    hip_rows = load("hipparcos_bright.csv")
    by_id = sum(1 for r in hip_rows if hip_to_gaia.get(r["hip"].strip()) in have_gaia)
    dedup = _positional_dedup(stars, hip_rows)
    print(f"  cross-match links: {by_id} by catalog id, {len(dedup)} by position")

    n_hip = 0
    for r in hip_rows:
        hip = r["hip"].strip()
        sid = hip_to_gaia.get(hip)
        if (sid and sid in have_gaia) or hip in dedup:
            continue
        got = rv_by_hip.get(hip)
        plx = num(r, "plx")
        if not got or plx is None or plx <= 0:
            continue
        rv, rv_err, main_id, sp, otype = got
        stars.append({
            "id": f"hip{hip}", "source_id": "", "hip": hip,
            "ra": num(r, "ra"), "dec": num(r, "dec"),
            "parallax": plx, "parallax_error": num(r, "e_plx"),
            "pmra": num(r, "pmra"), "pmdec": num(r, "pmdec"),
            "pmra_error": num(r, "e_pm_ra"), "pmdec_error": num(r, "e_pm_de"),
            "rv": rv, "rv_error": rv_err, "rv_source": "simbad", "rv_snr": None,
            "gmag": num(r, "hp_mag"), "bp_rp": None, "teff": None,
            "ruwe": None, "nss": 0,
            "astrom_source": "hipparcos2",
            "simbad_id": main_id, "sp_type": sp, "otype": otype,
            "featured": feat_hip.get(hip, ""),
        })
        n_hip += 1
    print(f"  Hipparcos-only bright stars added: {n_hip}")

    # ---- better published parallaxes, where they exist ----------------------
    n_plx = 0
    for st in stars:
        got = PUBLISHED_PARALLAX.get((st.get("hip") or "").strip())
        if not got:
            continue
        plx, err, src = got
        st["parallax_catalog"] = st["parallax"]
        st["parallax"], st["parallax_error"] = plx, err
        st["parallax_source"] = src
        n_plx += 1
    if n_plx:
        print(f"  adopted a published parallax for {n_plx} components "
              f"(see PUBLISHED_PARALLAX)")

    # ---- featured stars neither catalog can supply ------------------------
    # Scholz's star (the closest known past flyby) has DR3 photometry but no
    # astrometric solution, being a very red L/T binary. Fall back to SIMBAD's
    # discovery-paper astrometry rather than lose the object.
    present = {s["featured"] for s in stars if s["featured"]}
    feat_label = {}
    for r in feat:
        feat_label[r["query_id"]] = r["label"]
    # A featured entry may name a system whose components are already in the
    # catalog under different names (Luyten 726-8 is UV Cet + BL Cet; GJ 71 is
    # Tau Ceti). Adding it would duplicate them, so also check by position.
    # Genuine binaries such as 61 Cygni A/B stay - they are separate stars that
    # arrive on their own Gaia rows, not system-level aliases.
    from scipy.spatial import cKDTree

    def unit(ra_deg, dec_deg):
        ra, dec = np.radians(ra_deg), np.radians(dec_deg)
        return np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])

    existing = cKDTree(np.array([unit(s["ra"], s["dec"]) for s in stars]))
    chord = 2 * np.sin(np.radians(60.0 / 3600.0) / 2)     # 60 arcsec

    def already_there(ra, dec, plx, pmra, pmdec) -> bool:
        # SIMBAD positions are J2000, Gaia's are 2016.0. Luyten 726-8 moves
        # 3.3 arcsec per year, so by the Gaia epoch it has shifted ~54 arcsec -
        # far enough to miss its own components unless the epoch is aligned.
        dt = GAIA_EPOCH - 2000.0
        dec2 = dec + pmdec * dt / 3.6e6
        ra2 = ra + (pmra * dt / 3.6e6) / max(np.cos(np.radians(dec)), 1e-6)
        for i in existing.query_ball_point(unit(ra2, dec2), r=chord):
            p = stars[i]["parallax"]
            if p and abs(p - plx) / plx < 0.15:
                return True
        return False

    n_simbad_only = 0
    for r in load("featured_astrometry.csv"):
        label = feat_label.get(r["query_id"].strip(), "")
        if not label or label in present:
            continue
        plx, rv = num(r, "plx_value"), num(r, "rvz_radvel")
        pmra, pmdec = num(r, "pmra"), num(r, "pmdec")
        if (r.get("rvz_type") or "").strip() != "v" or rv is None or abs(rv) > RV_SANITY_KMS:
            continue
        if None in (plx, pmra, pmdec) or plx <= 0:
            continue
        if already_there(num(r, "ra"), num(r, "dec"), plx, pmra, pmdec):
            print(f"    skipping featured '{label}': already present as its components")
            continue
        present.add(label)
        stars.append({
            "id": "simbad" + label.replace(" ", "_"), "source_id": "", "hip": "",
            "ra": num(r, "ra"), "dec": num(r, "dec"),
            "parallax": plx, "parallax_error": num(r, "plx_err"),
            "pmra": pmra, "pmdec": pmdec,
            "rv": rv, "rv_error": num(r, "rvz_err"),
            "pmra_error": None, "pmdec_error": None,
            "rv_source": "simbad", "rv_snr": None,
            "gmag": None, "bp_rp": None, "teff": None, "ruwe": None, "nss": 0,
            "astrom_source": "simbad",
            "simbad_id": r["main_id"].strip(), "sp_type": r.get("sp_type", "").strip(),
            "otype": r.get("otype", "").strip(), "featured": label,
        })
        n_simbad_only += 1
    print(f"  SIMBAD-only featured stars added: {n_simbad_only}")

    # ---- bound systems: use the barycenter's velocity, not a component's ----
    for st in stars:
        st.setdefault("parallax_catalog", "")
        st.setdefault("parallax_source", "")
        st.setdefault("rv_component", "")
        st.setdefault("system_parent", "")
        st.setdefault("system", -1)
        # Hipparcos and SIMBAD rows publish no correlations; blank means
        # "assume independent", which is all that was ever available for them.
        for k in ("plx_pmra_corr", "plx_pmdec_corr", "pmra_pmdec_corr"):
            st.setdefault(k, "")
    apply_system_rv(stars)

    # ---- phase space + straight-line encounter solution ---------------------
    keys = ("ra", "dec", "parallax", "pmra", "pmdec", "rv")
    stars = [s for s in stars if all(s[k] is not None for k in keys)]
    arr = {k: np.array([s[k] for s in stars], dtype=float) for k in keys}
    pos, vel = icrs_to_galactic_cartesian(
        arr["ra"], arr["dec"], arr["parallax"], arr["pmra"], arr["pmdec"], arr["rv"]
    )
    d_now = np.linalg.norm(pos, axis=1)
    d_ph, t_ph = linear_perihelion(pos, vel)

    is_feat = np.array([bool(s["featured"]) for s in stars])
    # Keep the full 20 pc screen through integration. Whether a star is
    # scientifically interesting is decided after its orbit is integrated, not
    # from the straight-line value used to find it: deflection has been measured
    # moving a perihelion by 4 pc, so a 10 pc cut here would discard real close
    # passes before they were ever computed. Stage 6 narrows this for display.
    keep = (d_now < 30.0) | ((d_ph < 20.0) & (np.abs(t_ph) < 5.5)) | is_feat
    print(f"  after selection cut: {keep.sum()} of {len(stars)}")

    out, grades = [], {"A": 0, "B": 0, "C": 0}
    for i, s in enumerate(stars):
        if not keep[i]:
            continue
        rv_rel, astro, g, is_wd, m_g = grade(s)
        grades[g] += 1
        s = dict(s)
        s.update({
            "x": pos[i, 0], "y": pos[i, 1], "z": pos[i, 2],
            "vx": vel[i, 0], "vy": vel[i, 1], "vz": vel[i, 2],
            "dist_now_pc": d_now[i],
            "lin_dmin_pc": d_ph[i], "lin_tmin_myr": t_ph[i],
            "abs_g": m_g, "is_wd": int(is_wd),
            "rv_reliability": rv_rel, "astrom_flag": astro, "grade": g,
        })
        out.append(s)

    out.sort(key=lambda s: s["dist_now_pc"])
    with open(PROC / "merged.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    print(f"  quality grades: A={grades['A']}  B={grades['B']}  C={grades['C']}")
    print(f"  featured present: {sum(1 for s in out if s['featured'])} of {len(feat)}")
    print(f"  -> data/processed/merged.csv: {len(out)} stars")


if __name__ == "__main__":
    main()
