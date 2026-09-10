"""Stage 8 - verify that the published artifacts agree with each other.

Stage 5 validates the science. This validates the *build*: that the CSV, the web
bundle, the coefficient blob and the metadata all describe the same catalog.
Those are four files written by three scripts, and nothing else in the pipeline
would notice if one of them were stale - a half-finished export leaves a bundle
that loads perfectly and shows the wrong numbers.

Runs last, after stages 6 and 7, and writes data/processed/build_report.json.
Exits non-zero on any failure.
"""

from __future__ import annotations

import csv
import json
import pathlib

import numpy as np

BASE = pathlib.Path(__file__).parent.parent
PROC = BASE / "data" / "processed"
WEB = BASE / "web" / "data"


def main() -> None:
    orbits = list(csv.DictReader(open(PROC / "orbits.csv")))
    coeffs = np.load(PROC / "cheb_coeffs.npy")
    omit = json.loads((PROC / "orbit_meta.json").read_text())
    stars = json.loads((WEB / "stars.json").read_text())
    meta = json.loads((WEB / "meta.json").read_text())
    pub = list(csv.DictReader(open(BASE / "data" / "star_encounters.csv")))
    cheb_bytes = (WEB / "cheb.bin").stat().st_size

    n_ship = len(stars["names"])
    deg = meta["cheb_degree"]
    by_name = {s["name"]: s for s in orbits}

    def fnum(row, key):
        try:
            return float(row[key])
        except (KeyError, TypeError, ValueError):
            return float("nan")

    checks: list[tuple[str, bool, str]] = []

    def chk(label, ok, detail=""):
        checks.append((label, bool(ok), detail))

    chk("orbits.csv matches the integration metadata",
        len(orbits) == omit["n_stars"],
        f"{len(orbits)} rows vs n_stars={omit['n_stars']}")
    chk("coefficient array matches the catalog",
        coeffs.shape[0] == len(orbits),
        f"{coeffs.shape} vs {len(orbits)} rows")
    chk("public CSV covers the whole catalog",
        len(pub) == len(orbits), f"{len(pub)} rows vs {len(orbits)}")

    order = [fnum(r, "mc_dmin_med") for r in pub]
    chk("public CSV sorted by the integrated median",
        all(a <= b for a, b in zip(order, order[1:])),
        "ascending closest approach")

    lens = {k: len(v) for k, v in stars.items() if isinstance(v, list)}
    expect = {k: (n_ship * 3 if k == "color" else n_ship) for k in lens}
    chk("every stars.json array has the shipped length",
        lens == expect,
        f"{n_ship} stars; mismatched: "
        f"{ {k: lens[k] for k in lens if lens[k] != expect[k]} or 'none'}")

    want = n_ship * 3 * (deg + 1) * 4
    chk("cheb.bin is exactly the declared size",
        cheb_bytes == want, f"{cheb_bytes} bytes vs {want} expected")
    chk("metadata star count matches the bundle",
        meta.get("n_stars_shipped") == n_ship and meta.get("n_stars") == n_ship,
        f"meta {meta.get('n_stars_shipped')} vs bundle {n_ship}")
    chk("catalog count in metadata matches orbits.csv",
        meta.get("n_stars_catalog") == len(orbits),
        f"{meta.get('n_stars_catalog')} vs {len(orbits)}")

    # every shipped star is a real row, with the same numbers
    missing, drift = [], []
    for i, nm in enumerate(stars["names"]):
        row = by_name.get(nm)
        if row is None:
            missing.append(nm)
            continue
        if abs(fnum(row, "mc_dmin_med") - stars["mc_dmin_med"][i]) > 5e-4:
            drift.append(nm)
    chk("every shipped star exists in the catalog", not missing,
        f"{len(missing)} unknown, e.g. {missing[:3]}" if missing else "all resolved")
    chk("shipped closest approaches match the catalog", not drift,
        f"{len(drift)} disagree, e.g. {drift[:3]}" if drift else
        "agree to the exported rounding")

    # the censored flag bit must be derived from the censored fraction
    bad = [stars["names"][i] for i in range(n_ship)
           if bool(stars["flags"][i] & 256) != (stars["mc_censored"][i] > 0.5)]
    chk("censored flag bit agrees with the censored fraction", not bad,
        f"{len(bad)} disagree, e.g. {bad[:3]}" if bad else
        f"{sum(1 for f in stars['flags'] if f & 256)} flagged")

    # The rate is quoted in the interface, so the figure the page shows has to
    # be the figure the stage computed - not one left behind by an earlier run.
    rate_path = PROC / "encounter_rate.json"
    shipped = meta.get("encounter_rate", {})
    if rate_path.exists():
        src = json.loads(rate_path.read_text())
        same = all(shipped.get(k) == src.get(k) for k in
                   ("rate_within_1pc_per_myr", "stat_error", "sys_error",
                    "sigma_from_published"))
        chk("shipped encounter rate matches the stage that computed it", same,
            f"{shipped.get('rate_within_1pc_per_myr')} shipped vs "
            f"{src.get('rate_within_1pc_per_myr')} computed"
            if not same else
            f"{src['rate_within_1pc_per_myr']} +/- {src['stat_error']} per Myr, "
            f"{src['sigma_from_published']} sigma from published")
    else:
        chk("shipped encounter rate matches the stage that computed it",
            not shipped, "encounter_rate.json absent and none shipped"
            if not shipped else "metadata quotes a rate with no stage output")

    val = meta.get("validation", {})
    chk("embedded validation reports a pass",
        val.get("status") == "pass",
        f"status={val.get('status', 'absent')}, "
        f"{val.get('n_failed', '?')} of {val.get('n_checks', '?')} failed")

    print("  Build verification:")
    for label, ok, detail in checks:
        print(f"    [{'PASS' if ok else 'FAIL'}] {label}"
              + (f"\n           {detail}" if detail else ""))
    passed = all(ok for _, ok, _ in checks)

    (PROC / "build_report.json").write_text(json.dumps({
        "status": "pass" if passed else "fail",
        "n_stars_catalog": len(orbits),
        "n_stars_shipped": n_ship,
        "checks": [{"check": c, "passed": ok, "detail": d} for c, ok, d in checks],
    }, indent=2))
    print(f"\n  -> data/processed/build_report.json  "
          f"({'pass' if passed else 'FAIL'})")
    if not passed:
        raise SystemExit("build verification FAILED")


if __name__ == "__main__":
    main()
