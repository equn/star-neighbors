"""Bound multiple systems: grouping, and barycentric radial velocities.

Two components of one bound pair must not produce two different encounter
solutions, but that is exactly what happens if each is propagated on its own
cataloged velocity. Alpha Centauri A and B are listed 7.3 km/s apart in radial
velocity - not because the system is doing anything unusual, but because each
measurement catches the component at a different phase of an 80-year orbit
about the barycenter. Propagating those separately moved their computed closest
approach apart by 0.23 pc.

What actually travels through the Galaxy is the barycenter, and for many systems
SIMBAD carries a system-level entry with exactly that velocity: `* alf Cen` is
listed at -22.3 km/s, against components at -15.25 and -22.59. So bound groups
are identified first, and where the parent object has a usable velocity it is
adopted for every component.
"""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

# Angular separation below which two co-moving stars are treated as one system.
# 120 arcsec is wide enough for the resolved nearby binaries (61 Cygni is 30",
# Alpha Cen ~8") and far tighter than the ~1 pc typical field-star spacing.
SYSTEM_LINK_ARCSEC = 120.0
RV_DISAGREEMENT_KMS = 2.0


def group_systems(ra, dec, plx, pmra, pmdec, verbose: bool = True) -> list[int]:
    """Union-find over common-proper-motion pairs -> system index per star.

    Grouping on 3-D separation looks natural but fails on exactly the pairs that
    matter. Hipparcos gives Alpha Centauri A and B parallaxes that disagree by
    5.6%, which puts them 0.07 pc apart in Cartesian space although they orbit
    at 23 AU - the orbital motion corrupts each component's astrometric
    solution. The reliable test is the classic one: close on the sky, consistent
    parallax, and shared proper motion.
    """
    ra_r, dec_r = np.radians(ra), np.radians(dec)
    unit = np.stack([np.cos(dec_r) * np.cos(ra_r),
                     np.cos(dec_r) * np.sin(ra_r), np.sin(dec_r)], axis=-1)
    tree = cKDTree(unit)
    chord = 2 * np.sin(np.radians(SYSTEM_LINK_ARCSEC / 3600.0) / 2)

    parent = list(range(len(ra)))

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    pm = np.hypot(pmra, pmdec)
    linked = 0
    for i, j in tree.query_pairs(r=chord):
        if abs(plx[i] - plx[j]) / max(plx[i], plx[j]) > 0.25:
            continue                                   # different distances
        dpm = np.hypot(pmra[i] - pmra[j], pmdec[i] - pmdec[j])
        if dpm > 0.30 * max(pm[i], pm[j], 1e-6):
            continue                                   # not co-moving
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[ri] = rj
            linked += 1
    if verbose:
        print(f"  linked {linked} common-proper-motion pairs")
    return [find(i) for i in range(len(ra))]


def dense_ids(raw: list[int]) -> list[int]:
    """Renumber group roots, giving -1 to stars that are on their own."""
    counts: dict[int, int] = {}
    for r in raw:
        counts[r] = counts.get(r, 0) + 1
    remap: dict[int, int] = {}
    out = []
    for r in raw:
        if counts[r] < 2:
            out.append(-1)
            continue
        if r not in remap:
            remap[r] = len(remap)
        out.append(remap[r])
    return out


def rv_spread(sys_id: list[int], rv: list[float]) -> dict[int, float]:
    """Largest radial-velocity disagreement within each system, km/s."""
    by: dict[int, list[float]] = {}
    for s, v in zip(sys_id, rv):
        if s >= 0 and v is not None:
            by.setdefault(s, []).append(v)
    return {k: (max(v) - min(v)) for k, v in by.items() if len(v) > 1}
