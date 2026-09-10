"""Compare this catalog against the peer-reviewed Bailer-Jones encounter catalogs.

`compare_published.py` checks 47 stars against a hand-transcribed compilation.
This checks thousands against the refereed source that compilation draws on,
joined on Gaia source_id rather than on sky position, so there is no matching
step to get wrong.

Two reference catalogs, deliberately kept apart:

  BJ2022  Bailer-Jones (2022), ApJL 935, L9 - VizieR J/ApJ/935/L9/table12.
          61 stars with a median perihelion inside 1 pc, computed from Gaia
          DR3. Same data release as this project, keyed by DR3 source_id, and
          covering exactly the regime the project is about. This is the
          primary comparison.

  BJ2018  Bailer-Jones et al. (2018), A&A 616, A37 - VizieR J/A+A/616/A37/table23.
          3379 stars with a median perihelion inside 10 pc, computed from Gaia
          DR2, keyed by DR2 source_id and translated through
          gaiadr3.dr2_neighbourhood. Two data releases apart, so a
          disagreement here conflates method with DR2-vs-DR3 astrometry. Used
          for bulk distribution statistics, not for adjudicating single stars.

Both quote a median perihelion distance and epoch from their own Monte Carlo
sampling, so the comparison is median-to-median with an interval-overlap test,
the same convention `compare_published.py` uses.

Stars whose encounter is clipped by this project's +/-5 Myr window are reported
separately rather than counted as disagreements: the two calculations are not
answering the same question there.

    python3 scripts/compare_bailer_jones.py [--csv out.csv] [--refresh]
"""

from __future__ import annotations

import csv
import math
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import tap  # noqa: E402

BASE = pathlib.Path(__file__).parent.parent
RAW = BASE / "data" / "raw"
PROC = BASE / "data" / "processed"

BJ2022_TABLE = "J/ApJ/935/L9/table12"
BJ2018_TABLE = "J/A+A/616/A37/table23"

# The window this project integrates. Beyond it our perihelion is a bound, not
# an encounter, and comparing it to a reference that integrated further is
# comparing two different quantities.
T_SPAN_KYR = 5000.0

# dr2_neighbourhood returns several DR3 candidates for a DR2 source when the
# solution was split or merged between releases. Accept a translation only if
# one candidate is clearly the same star.
XMATCH_ARCSEC = 1.0
XMATCH_DMAG = 0.5


def cached(name: str, fetch, refresh: bool = False) -> list[dict]:
    """Fetch once and keep it on disk; these tables are static publications."""
    path = RAW / name
    if path.exists() and not refresh:
        with open(path) as fh:
            return list(csv.DictReader(fh))
    print(f"  fetching {name} ...")
    rows = fetch()
    if not rows:
        raise SystemExit(f"{name}: query returned nothing")
    RAW.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"    -> {path.relative_to(BASE)}  ({len(rows)} rows)")
    return rows


def fetch_bj2022(refresh: bool = False) -> list[dict]:
    return cached(
        "bj2022_dr3_encounters.csv",
        lambda: tap.sync_query(f'SELECT * FROM "{BJ2022_TABLE}"',
                               base=tap.VIZIER_TAP, timeout=300),
        refresh,
    )


def fetch_bj2018(refresh: bool = False) -> list[dict]:
    return cached(
        "bj2018_dr2_encounters.csv",
        lambda: tap.sync_query(f'SELECT * FROM "{BJ2018_TABLE}"',
                               base=tap.VIZIER_TAP, timeout=600),
        refresh,
    )


def fetch_dr2_to_dr3(dr2_ids: list[str], refresh: bool = False) -> dict[str, str]:
    """Translate DR2 source_ids to DR3, keeping only unambiguous matches."""

    def go() -> list[dict]:
        out: list[dict] = []
        # An IN list of 3379 ids overflows what the sync endpoint will accept,
        # and an upload job needs a VOTable; chunking is the plain way.
        for i in range(0, len(dr2_ids), 500):
            chunk = dr2_ids[i:i + 500]
            q = ("SELECT dr2_source_id, dr3_source_id, angular_distance, "
                 "magnitude_difference FROM gaiadr3.dr2_neighbourhood "
                 "WHERE dr2_source_id IN (%s)" % ",".join(chunk))
            out.extend(tap.sync_query(q, timeout=600))
            print(f"    translated {min(i + 500, len(dr2_ids))}/{len(dr2_ids)}")
        return out

    rows = cached("bj2018_dr2_to_dr3.csv", go, refresh)

    by_dr2: dict[str, list[dict]] = {}
    for r in rows:
        by_dr2.setdefault(r["dr2_source_id"], []).append(r)

    out: dict[str, str] = {}
    for dr2, cands in by_dr2.items():
        def score(c):
            return (float(c["angular_distance"] or 9e9),
                    abs(float(c["magnitude_difference"] or 9e9)))
        cands = sorted(cands, key=score)
        sep, dmag = score(cands[0])
        if sep <= XMATCH_ARCSEC and dmag <= XMATCH_DMAG:
            out[dr2] = cands[0]["dr3_source_id"]
    return out


def load_ours() -> dict[str, dict]:
    with open(PROC / "orbits.csv") as fh:
        return {r["source_id"]: r for r in csv.DictReader(fh) if r.get("source_id")}


def fnum(row: dict, key: str):
    try:
        v = float(row[key])
        return v if v == v else None
    except (TypeError, ValueError, KeyError):
        return None


def compare(ref: list[dict], ours: dict[str, dict], label: str, note: str,
            id_of, dmin_of, tmin_of, lo_of, hi_of) -> list[dict]:
    """Median-to-median comparison of one reference catalog against ours."""
    rows, unmatched, censored = [], 0, []
    for r in ref:
        sid = id_of(r)
        s = ours.get(sid) if sid else None
        if s is None:
            unmatched += 1
            continue
        ref_d, ref_t = dmin_of(r), tmin_of(r)
        our_d = fnum(s, "mc_dmin_med")
        our_t = fnum(s, "mc_tmin_med")
        if ref_d is None or our_d is None or our_d <= 0 or ref_d <= 0:
            continue
        our_t_kyr = our_t * 1000.0 if our_t is not None else None
        cens = fnum(s, "mc_censored") or 0.0

        rec = {
            "name": s["name"], "source_id": sid,
            "ref_dmin_pc": round(ref_d, 4), "our_dmin_pc": round(our_d, 4),
            "ref_lo_pc": lo_of(r), "ref_hi_pc": hi_of(r),
            "our_lo_pc": fnum(s, "mc_dmin_lo"), "our_hi_pc": fnum(s, "mc_dmin_hi"),
            "d_ratio": round(our_d / ref_d, 3),
            "ref_tmin_kyr": ref_t, "our_tmin_kyr": round(our_t_kyr, 1)
                if our_t_kyr is not None else None,
            "mc_censored": round(cens, 3),
            "grade": s.get("grade", ""),
            "sys_dmin_pc": fnum(s, "sys_dmin_pc"),
        }
        # Outside our window the two calculations answer different questions.
        rec["window_limited"] = bool(
            cens > 0.5 or (ref_t is not None and abs(ref_t) > T_SPAN_KYR))
        (censored if rec["window_limited"] else rows).append(rec)

    print(f"\n{'=' * 74}\n{label}\n{note}\n{'=' * 74}")
    print(f"  reference rows           {len(ref)}")
    print(f"  joined to our catalog  {len(rows) + len(censored)}"
          f"   (unmatched {unmatched})")
    print(f"  comparable               {len(rows)}"
          f"   (outside our +/-5 Myr window: {len(censored)})")
    if not rows:
        return rows

    ratios = sorted(r["d_ratio"] for r in rows)
    n = len(ratios)

    def pct(p):
        return ratios[min(n - 1, int(p * n))]

    overlap = sum(
        1 for r in rows
        if (r["our_lo_pc"] is not None and r["our_hi_pc"] is not None
            and r["our_lo_pc"] - 1e-9 <= r["ref_dmin_pc"] <= r["our_hi_pc"] + 1e-9)
        or (r["ref_lo_pc"] is not None and r["ref_hi_pc"] is not None
            and r["ref_lo_pc"] - 1e-9 <= r["our_dmin_pc"] <= r["ref_hi_pc"] + 1e-9))

    ep = [r for r in rows if r["ref_tmin_kyr"] is not None
          and r["our_tmin_kyr"] is not None]
    ep_ok = sum(1 for r in ep
                if abs(r["our_tmin_kyr"] - r["ref_tmin_kyr"])
                <= max(30.0, 0.05 * abs(r["ref_tmin_kyr"])))
    within10 = sum(1 for x in ratios if 0.9 <= x <= 1.1)

    print(f"\n  perihelion distance, ours / reference")
    print(f"    median   {pct(0.5):.3f}")
    print(f"    16-84%   {pct(0.16):.3f} - {pct(0.84):.3f}")
    print(f"    5-95%    {pct(0.05):.3f} - {pct(0.95):.3f}")
    print(f"    within 10%          {within10}/{n}  ({100 * within10 / n:.1f}%)")
    print(f"    intervals overlap   {overlap}/{n}  ({100 * overlap / n:.1f}%)")
    if ep:
        print(f"  epoch agrees to 5% or 30 kyr  {ep_ok}/{len(ep)}"
              f"  ({100 * ep_ok / len(ep):.1f}%)")

    worst = sorted(rows, key=lambda r: -abs(math.log(r["d_ratio"] or 1)))[:12]
    print(f"\n  largest disagreements")
    print(f"    {'star':<30}{'ref':>8}{'ours':>8}{'ratio':>7}"
          f"{'ref t':>9}{'our t':>9}  gr  Galaxy+/-")
    for r in worst:
        sysv = f"{r['sys_dmin_pc']:.3f}" if r["sys_dmin_pc"] is not None else "  -  "
        rt = f"{r['ref_tmin_kyr']:.0f}" if r["ref_tmin_kyr"] is not None else "-"
        ot = f"{r['our_tmin_kyr']:.0f}" if r["our_tmin_kyr"] is not None else "-"
        print(f"    {r['name'][:29]:<30}{r['ref_dmin_pc']:8.4f}"
              f"{r['our_dmin_pc']:8.4f}{r['d_ratio']:7.2f}"
              f"{rt:>9}{ot:>9}   {r['grade']}  {sysv}")

    if censored:
        print(f"\n  outside our window, reported not scored ({len(censored)})")
        for r in sorted(censored, key=lambda x: x["ref_dmin_pc"])[:8]:
            rt = f"{r['ref_tmin_kyr']:.0f}" if r["ref_tmin_kyr"] is not None else "-"
            print(f"    {r['name'][:29]:<30}ref {r['ref_dmin_pc']:7.4f} pc "
                  f"at {rt:>8} kyr   censored {r['mc_censored']:.2f}")
    return rows + censored


def gravity_budget(ref: list[dict], ours: dict[str, dict], id_of, dmin_of) -> None:
    """Ask whether a disagreement is even dynamically reachable.

    For every star this project ships both a straight-line perihelion and an
    integrated one, so |integrated - straight line| is the entire budget the
    Galactic potential has to work with over the +/-5 Myr window. A reference
    value further from our straight line than that budget cannot be explained
    by any choice of potential: switching gravity off completely would not
    close the gap, so the difference has to live in the astrometry or in the
    definition of the quantity, not in the dynamics.
    """
    tested = covered = reachable = 0
    beyond = []
    for r in ref:
        s = ours.get(id_of(r))
        if s is None:
            continue
        ref_d = dmin_of(r)
        lin, gal = fnum(s, "lin_dmin_pc"), fnum(s, "gal_dmin_pc")
        lo, hi = fnum(s, "mc_dmin_lo"), fnum(s, "mc_dmin_hi")
        sysd = fnum(s, "sys_dmin_pc") or 0.0
        cens = fnum(s, "mc_censored") or 0.0
        if not (ref_d and lin and gal) or lin <= 0 or cens > 0.5 or lo is None:
            continue
        tested += 1
        # First ask the ordinary question: is the reference inside our stated
        # uncertainty, widened by the Galaxy-model envelope we already publish?
        if lo - sysd - 1e-9 <= ref_d <= hi + sysd + 1e-9:
            covered += 1
            continue
        # It is not. Gravity's entire contribution over the window is
        # |gal - lin|; allow it to act in either direction, and allow the
        # measurement interval on top. Anything outside that is unreachable by
        # any potential, because switching gravity off does not close it either.
        budget = abs(gal - lin)
        half = (hi - lo) / 2.0
        span_lo = min(lin, gal) - budget - half - sysd
        span_hi = max(lin, gal) + budget + half + sysd
        if span_lo <= ref_d <= span_hi:
            reachable += 1
        else:
            beyond.append((abs(ref_d - gal), s["name"], ref_d, lin, gal,
                           budget, half + sysd))

    print(f"\n  is the disagreement dynamically reachable?")
    print(f"    stars tested                                     {tested}")
    print(f"    reference inside our interval + Galaxy envelope  {covered}")
    print(f"    outside, but within the gravity budget           {reachable}")
    print(f"    beyond what ANY potential could produce          {len(beyond)}")
    if beyond:
        print(f"\n    {'star':<30}{'ref':>8}{'our lin':>9}{'our gal':>9}"
              f"{'gravity':>9}{'meas':>8}{'gap':>8}")
        for gap, nm, rd, lin, gal, budget, meas in sorted(beyond, reverse=True)[:10]:
            print(f"    {nm[:29]:<30}{rd:8.4f}{lin:9.4f}{gal:9.4f}"
                  f"{budget:9.4f}{meas:8.4f}{gap:8.4f}")


def main() -> None:
    refresh = "--refresh" in sys.argv
    ours = load_ours()
    print(f"  our catalog: {len(ours)} stars carrying a Gaia source_id")

    all_rows = []

    bj22 = fetch_bj2022(refresh)
    all_rows += [dict(r, catalog="BJ2022") for r in compare(
        bj22, ours,
        "BJ2022 - Bailer-Jones (2022), ApJL 935, L9",
        "Gaia DR3, 61 stars inside 1 pc. Same data release, exact source_id join.",
        id_of=lambda r: r.get("GaiaDR3"),
        dmin_of=lambda r: fnum(r, "dphmed"),
        tmin_of=lambda r: fnum(r, "tphmed"),
        lo_of=lambda r: fnum(r, "dph5"), hi_of=lambda r: fnum(r, "dph95"))]
    gravity_budget(bj22, ours, lambda r: r.get("GaiaDR3"),
                   lambda r: fnum(r, "dphmed"))

    bj18 = fetch_bj2018(refresh)
    xm = fetch_dr2_to_dr3([r["ID"] for r in bj18], refresh)
    print(f"\n  DR2->DR3: {len(xm)}/{len(bj18)} unambiguous translations")
    all_rows += [dict(r, catalog="BJ2018") for r in compare(
        bj18, ours,
        "BJ2018 - Bailer-Jones et al. (2018), A&A 616, A37",
        "Gaia DR2, 3379 stars inside 10 pc, ids translated through "
        "dr2_neighbourhood.\nA different data release: read the spread, not "
        "individual stars.",
        id_of=lambda r: xm.get(r["ID"]),
        dmin_of=lambda r: fnum(r, "dphmed"),
        tmin_of=lambda r: fnum(r, "tphmed"),
        lo_of=lambda r: fnum(r, "b_dph"), hi_of=lambda r: fnum(r, "B_dph"))]

    if "--csv" in sys.argv:
        out = pathlib.Path(sys.argv[sys.argv.index("--csv") + 1])
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(all_rows[0].keys()))
            w.writeheader()
            w.writerows(all_rows)
        print(f"\n  -> {out}")


if __name__ == "__main__":
    main()
