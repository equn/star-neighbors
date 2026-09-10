"""Stage 3 - attach human-readable names from SIMBAD.

Gaia source_ids are useless in a UI. SIMBAD holds every alias each object has;
this picks the most recognizable one by catalog priority, so a star shows up
as "Ross 248" or "Barnard's Star" rather than "gaia1926461164913660160".
"""

from __future__ import annotations

import csv
import json
import pathlib
import re
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from tap import SIMBAD_TAP, sync_query  # noqa: E402

BASE = pathlib.Path(__file__).parent.parent
PROC = BASE / "data" / "processed"

CHUNK = 400

GREEK = {
    "alf": "α", "bet": "β", "gam": "γ", "del": "δ",
    "eps": "ε", "zet": "ζ", "eta": "η", "tet": "θ",
    "iot": "ι", "kap": "κ", "lam": "λ", "mu.": "μ",
    "nu.": "ν", "ksi": "ξ", "omi": "ο", "pi.": "π",
    "rho": "ρ", "sig": "σ", "tau": "τ", "ups": "υ",
    "phi": "φ", "chi": "χ", "psi": "ψ", "ome": "ω",
}

# Lower rank = preferred. Proper names first, then the catalogs people
# actually recognize, then survey designations, with Gaia last.
#
# The split between ranks 5-7 matters. A proper-motion catalog number says
# something about the object - LSPM covers stars moving faster than 0.15 arcsec
# per year, which is most of why a star turns up in an encounter list at all -
# whereas TYC/UCAC/2MASS are positional survey identifiers that every point
# source has. Ties are broken by length, so leaving them in one tier let
# "TYC 3135-52-1" beat "LSPM J1935+3746" purely for being shorter.
#
# Survey names are also coordinate snapshots, which ages badly for exactly these
# stars: 2MASS J21462206+3813047 encodes where LSPM J2146+3813 sat around 1999,
# and at 0.2 arcsec/yr it has since moved several arcseconds off that position.
PRIORITY = [
    (0, re.compile(r"^NAME\s+", re.I)),
    (1, re.compile(r"^\*\s")),            # Bayer / Flamsteed
    (2, re.compile(r"^V\*\s")),           # variable-star designation
    (3, re.compile(r"^(GJ|Gl|Wolf|Ross|LHS|LAL|LFT|LTT|LP|L)\s", re.I)),
    (4, re.compile(r"^(HD|HIP|HR|BD|CD|CPD|SAO)\s", re.I)),
    (5, re.compile(r"^(LSPM|UPM|EGGR|GD|WD|G)\s", re.I)),      # proper-motion / nearby-star
    (6, re.compile(r"^(PM|StKM)\s", re.I)),                    # other PM designations
    (7, re.compile(r"^(TYC|UCAC|2MASS|WISE|USNO|LAMOST|SDSS|WDS|SCR|RX)\s", re.I)),
    # Running numbers last: CNS5, Karmn and TIC identify an object but carry no
    # position and mean nothing to a reader, so a coordinate-based survey name
    # is more useful even though it is longer. "** ABC 12A" is a double-star
    # discoverer code, equally opaque - the WDS coordinate name beats it.
    # Gaia is deliberately absent so it falls to the default rank, last of all.
    (8, re.compile(r"^(CNS5|Karmn|TIC)\s", re.I)),
    (8, re.compile(r"^\*\*\s")),
]


def squash(s: str) -> str:
    return " ".join(s.split())


def prettify(ident: str) -> str:
    """Turn a SIMBAD identifier into something readable."""
    s = squash(ident)
    if s.upper().startswith("NAME "):
        return s[5:]
    if s.startswith("* "):
        body = s[2:]
        parts = body.split()
        if parts and parts[0].lower() in GREEK:
            parts[0] = GREEK[parts[0].lower()]
        return " ".join(parts)
    if s.startswith("V* "):
        return s[3:]
    if s.startswith("** "):
        return s[3:]
    return s


def rank(ident: str) -> int:
    s = squash(ident)
    for r, pat in PRIORITY:
        if pat.match(s):
            return r
    return 9


def chunks(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def main() -> None:
    with open(PROC / "merged.csv") as fh:
        stars = list(csv.DictReader(fh))

    # query key -> star index
    want: dict[str, int] = {}
    for i, s in enumerate(stars):
        if s["source_id"]:
            want[f"Gaia DR3 {s['source_id']}"] = i
        elif s["hip"]:
            want[f"HIP {s['hip']}"] = i
    keys = list(want)
    cache_path = PROC / "simbad_cache.json"
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else None
    if cache:
        print(f"  using cached SIMBAD lookups ({len(cache['aliases'])} objects)")
        oid_by_key = cache["oid_by_key"]
        aliases = {int(k): v for k, v in cache["aliases"].items()}
        oid_of = {i: oid_by_key[k] for k, i in want.items() if k in oid_by_key}
        return finish(stars, oid_of, aliases)
    print(f"  resolving {len(keys)} identifiers in chunks of {CHUNK}")

    # ---- pass 1: query key -> SIMBAD oid ------------------------------------
    oid_of: dict[int, int] = {}          # star index -> oid
    for n, part in enumerate(chunks(keys, CHUNK), 1):
        quoted = ",".join("'" + k.replace("'", "''") + "'" for k in part)
        rows = sync_query(
            f"SELECT i.id, i.oidref FROM ident AS i WHERE i.id IN ({quoted})",
            base=SIMBAD_TAP,
        )
        for r in rows:
            idx = want.get(squash(r["id"]))
            if idx is None:                       # spacing differences
                idx = want.get(squash(r["id"]).replace("Gaia DR3 ", "Gaia DR3 "))
            if idx is not None:
                oid_of[idx] = int(r["oidref"])
        print(f"    chunk {n}/{-(-len(keys)//CHUNK)}: {len(oid_of)} matched so far")

    # ---- pass 2: all aliases for those objects ------------------------------
    oids = sorted(set(oid_of.values()))
    aliases: dict[int, list[str]] = {}
    for n, part in enumerate(chunks(oids, CHUNK), 1):
        lst = ",".join(str(o) for o in part)
        rows = sync_query(
            f"SELECT i.oidref, i.id FROM ident AS i WHERE i.oidref IN ({lst})",
            base=SIMBAD_TAP,
        )
        for r in rows:
            aliases.setdefault(int(r["oidref"]), []).append(r["id"])
        if n % 5 == 0 or n == 1:
            print(f"    alias chunk {n}/{-(-len(oids)//CHUNK)}")

    cache_path.write_text(json.dumps({
        "oid_by_key": {k: oid_of[i] for k, i in want.items() if i in oid_of},
        "aliases": {str(k): v for k, v in aliases.items()},
    }))
    return finish(stars, oid_of, aliases)


def finish(stars, oid_of, aliases) -> None:
    # ---- choose display names -----------------------------------------------
    named = 0
    for i, s in enumerate(stars):
        alist = aliases.get(oid_of.get(i, -1), [])
        if alist:
            best = min(alist, key=lambda a: (rank(a), len(squash(a))))
            s["name"] = prettify(best)
            s["aliases"] = "|".join(sorted({squash(a) for a in alist}))
            named += 1
        else:
            s["name"] = ""
            s["aliases"] = ""
        # the curated label always wins if present
        if s["featured"]:
            s["name"] = s["featured"]
        if not s["name"]:
            s["name"] = f"Gaia DR3 {s['source_id']}" if s["source_id"] else f"HIP {s['hip']}"

    out = PROC / "named.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(stars[0].keys()))
        w.writeheader()
        w.writerows(stars)
    print(f"  named {named} of {len(stars)} from SIMBAD")
    print(f"  -> {out.relative_to(BASE)}")


if __name__ == "__main__":
    main()
