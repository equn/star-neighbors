"""Completeness gates shared by the two exporters.

Both `06_export.py` and `07_csv.py` publish from `orbits.csv`, and both used to
tolerate a half-finished science run in their own way: one turned missing Monte
Carlo values into zeros, the other wrote blank cells under a header that
promised them. Neither failure is detectable downstream - the bundle loads, the
CSV parses, and every number is wrong.

The rule lives here so the two cannot drift apart: a row is publishable only if
every required column exists, is non-empty, parses as a finite number where a
number is expected, and carries the integrated uncertainty model.
"""

from __future__ import annotations

import math

# Columns without which an export would be a lie: the shipped metadata states
# the integrated Monte Carlo unconditionally, so the data had better contain it.
NUMERIC = (
    "dist_now_pc", "gal_dmin_pc", "gal_tmin_myr",
    "mc_dmin_med", "mc_dmin_lo", "mc_dmin_hi",
    "mc_tmin_med", "mc_tmin_lo", "mc_tmin_hi",
    "mc_censored", "sys_dmin_pc", "zp_dmin_pc",
)
TEXT = ("name", "grade", "mc_model")
REQUIRED = NUMERIC + TEXT


def require_complete(rows: list[dict], what: str) -> None:
    """Raise SystemExit unless every row is publishable. `what` names the caller."""
    if not rows:
        raise SystemExit(f"{what} refused: orbits.csv is empty.")

    missing = [c for c in REQUIRED if c not in rows[0]]
    if missing:
        raise SystemExit(f"{what} refused: orbits.csv has no {missing}. "
                         f"Run the missing pipeline stage first.")

    bad_model = sum(1 for r in rows if r.get("mc_model") != "integrated")
    if bad_model:
        raise SystemExit(f"{what} refused: {bad_model} of {len(rows)} rows are not "
                         f"mc_model='integrated'. The uncertainty stage did not finish.")

    blank, nonfinite, example = 0, 0, ""
    for r in rows:
        for c in REQUIRED:
            v = r.get(c, "")
            if v is None or str(v).strip() == "":
                blank += 1
                example = example or f"{r.get('name', '?')}.{c} is empty"
                break
            if c in NUMERIC:
                try:
                    x = float(v)
                except (TypeError, ValueError):
                    x = float("nan")
                if not math.isfinite(x):
                    nonfinite += 1
                    example = example or f"{r.get('name', '?')}.{c} = {v!r}"
                    break
    if blank or nonfinite:
        raise SystemExit(
            f"{what} refused: {blank} rows with an empty required field and "
            f"{nonfinite} with a non-finite number (e.g. {example}).")
