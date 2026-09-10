"""Resolve the curated featured-star names to catalog identifiers via SIMBAD."""

from __future__ import annotations

import csv
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from featured import FEATURED  # noqa: E402
from tap import SIMBAD_TAP, sync_query  # noqa: E402

RAW = pathlib.Path(__file__).parent.parent / "data" / "raw"


def main() -> None:
    names = list(FEATURED)
    quoted = ",".join("'" + n.replace("'", "''") + "'" for n in names)
    q = f"""
    SELECT i.id AS query_id, b.main_id, b.otype, b.sp_type, b.plx_value,
           g.id AS gaia_id, h.id AS hip_id
    FROM ident AS i
    JOIN basic AS b ON b.oid = i.oidref
    LEFT OUTER JOIN ident AS g ON (g.oidref = b.oid AND g.id LIKE 'Gaia DR3 %')
    LEFT OUTER JOIN ident AS h ON (h.oidref = b.oid AND h.id LIKE 'HIP %')
    WHERE i.id IN ({quoted})
    """
    rows = sync_query(q, base=SIMBAD_TAP)

    # SIMBAD stores identifiers with its own spacing ("Ross  248", "HD   7977"),
    # so compare on a whitespace- and case-normalized key rather than verbatim.
    def norm(s: str) -> str:
        return " ".join(s.split()).upper()

    by_norm = {norm(n): lbl for n, lbl in FEATURED.items()}

    out, seen = [], set()
    for r in rows:
        qid = r["query_id"].strip()
        key = norm(qid)
        if key in seen:
            continue
        seen.add(key)
        out.append({
            "query_id": qid,
            "label": by_norm.get(key, r["main_id"].strip()),
            "main_id": r["main_id"].strip(),
            "otype": r["otype"].strip(),
            "sp_type": r["sp_type"].strip(),
            "gaia_source_id": (r["gaia_id"] or "").replace("Gaia DR3", "").strip(),
            "hip": (r["hip_id"] or "").replace("HIP", "").strip(),
        })

    missing = [n for n in names if norm(n) not in seen]
    with open(RAW / "featured_resolved.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    print(f"  -> featured_resolved.csv: {len(out)} resolved")
    if missing:
        print(f"  !! unresolved ({len(missing)}): {', '.join(missing)}")


if __name__ == "__main__":
    main()
