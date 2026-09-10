"""Tests for the export completeness gate.

The gate exists to stop a half-finished science run reaching the published
artifacts, which means the interesting behavior is the *rejections* - and a
gate that silently stopped rejecting would look exactly like a gate that was
working. These cases were run by hand once; keeping them here means they run
every time.

    python3 scripts/test_gates.py            # synthetic cases only
    python3 scripts/test_gates.py --catalog  # also gate the real orbits.csv

No test framework: the pipeline's only dependencies are numpy and scipy, and a
test that needs an install is a test that does not get run.
"""

from __future__ import annotations

import csv
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from gates import NUMERIC, REQUIRED, require_complete  # noqa: E402

BASE = pathlib.Path(__file__).parent.parent


def good_row(i: int = 0) -> dict:
    row = {c: "1.0" for c in NUMERIC}
    row.update(name=f"Test {i}", grade="A", mc_model="integrated")
    return row


def rows(n: int = 5) -> list[dict]:
    return [good_row(i) for i in range(n)]


def edit(n: int, col: str, value) -> list[dict]:
    """A clean set with one cell replaced - one defect, one row."""
    out = rows()
    if value is None:
        del out[n][col]
    else:
        out[n][col] = value
    return out


def drop_column(col: str) -> list[dict]:
    return [{k: v for k, v in r.items() if k != col} for r in rows()]


# (label, rows, must_be_rejected)
CASES = [
    ("clean run", rows(), False),
    ("single clean row", rows(1), False),
    ("negative epoch is legal", edit(0, "mc_tmin_med", "-4.2"), False),
    ("zero epoch is legal", edit(1, "mc_tmin_med", "0"), False),
    ("zero distance is legal", edit(2, "mc_dmin_lo", "0.0"), False),
    ("zero censored fraction is legal", edit(3, "mc_censored", "0"), False),
    ("scientific notation is legal", edit(0, "mc_dmin_lo", "1.5e-3"), False),

    ("empty input", [], True),
    ("blank mc_dmin_med", edit(2, "mc_dmin_med", ""), True),
    ("whitespace-only grade", edit(1, "grade", "   "), True),
    ("missing key on one row", edit(3, "mc_censored", None), True),
    ("column absent from every row", drop_column("mc_dmin_hi"), True),
    ("non-numeric where a number belongs", edit(0, "mc_tmin_lo", "n/a"), True),
    ("literal nan", edit(4, "mc_tmin_med", "nan"), True),
    ("positive inf", edit(0, "mc_dmin_hi", "inf"), True),
    ("negative inf", edit(2, "mc_dmin_lo", "-inf"), True),
    ("straight-line model on one row", edit(1, "mc_model", "straight_line"), True),
    ("empty model on one row", edit(0, "mc_model", ""), True),
    # deleting from row 0 trips the header check rather than the per-row one,
    # which is the point: the two paths report different things
    ("key deleted from the first row", edit(0, "name", None), True),
]


def check(label: str, data: list[dict], expect_reject: bool) -> bool:
    try:
        require_complete(data, "test")
        rejected, detail = False, ""
    except SystemExit as exc:
        rejected, detail = True, str(exc)
    ok = rejected == expect_reject
    verb = "rejected" if rejected else "accepted"
    print(f"    [{'PASS' if ok else 'FAIL'}] {label}: {verb}"
          + (f"\n           {detail}" if rejected and ok else "")
          + ("" if ok else f"  (expected {'reject' if expect_reject else 'accept'})"))
    return ok


def main() -> None:
    bad = sum(1 for _, _, reject in CASES if reject)
    print(f"  {len(CASES)} gate cases over {len(REQUIRED)} required columns "
          f"({bad} defects that must be caught, {len(CASES) - bad} valid inputs "
          f"that must not be):")
    passed = all(check(*c) for c in CASES)

    if "--catalog" in sys.argv:
        path = BASE / "data" / "processed" / "orbits.csv"
        if path.exists():
            with open(path) as fh:
                real = list(csv.DictReader(fh))
            print(f"\n  Real catalog ({len(real)} rows):")
            passed &= check("the shipped catalog passes", real, False)
        else:
            print(f"\n  {path} absent; skipping the catalog case")

    print(f"\n  {'all gate cases pass' if passed else 'GATE TESTS FAILED'}")
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
