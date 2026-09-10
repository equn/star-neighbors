"""Compare this catalog against the published list of sub-5-light-year encounters.

The reference is the compilation of close stellar encounters maintained on
Wikipedia's "List of nearest stars" (drawing on Bailer-Jones et al. 2018/2022 and
Bobylev & Bajkova), which quotes a perihelion distance, an epoch and a present
distance for every star known to pass within 5 ly of the Sun.

This is the strongest validation available to the project. Everything else it
checks is internal - that the artifacts agree with each other, that the numbers
are finite and ordered - or compares a handful of stars against two reference
charts. This compares the *result* against an independent, peer-reviewed
compilation, star by star, on 47 objects.

Matching is by sky position, not by name: the reference abbreviates 2MASS
designations ("2MASS J0628+1845") while the catalog carries the full form, so
names would produce spurious misses.

Position alone is not enough either. The reference quotes positions at a mix of
epochs and the catalog is at Gaia's 2016.0, and these are by construction
high-proper-motion stars: Barnard's Star moves 10.3 arcsec a year, so sixteen
years of it is 165 arcsec and a naive match silently "loses" the most famous
entries in the table. Each catalog star is therefore also propagated back to
J2000 and the better of the two positions is used, with a present-distance
agreement test as a guard against matching the wrong object.

    python3 scripts/compare_published.py [--csv out.csv]
"""

from __future__ import annotations

import csv
import math
import pathlib
import re
import sys

BASE = pathlib.Path(__file__).parent.parent
PROC = BASE / "data" / "processed"
LY_PER_PC = 3.261563777
MATCH_ARCSEC = 30.0

# name, d_min ly, +err, -err, epoch kyr, present distance ly, RA, Dec
REFERENCE = [
    ("Gliese 710",            0.167, 0.012, 0.012,   1296.0,  62.248, "18 19 50.843", "-01 56 18.98"),
    ("HD 7977",               0.478, 0.104, 0.078,  -2764.0, 246.740, "01 20 31.597", "+61 52 57.08"),
    ("Scholz's Star",         0.82,  0.37,  0.22,     -78.5,  22.200, "07 20 03.20",  "-08 46 51.2"),
    ("2MASS J0628+1845",      1.61,  0.28,  0.24,    1720.0, 272.280, "06 28 11.593", "+18 45 12.91"),
    ("2MASS J0805+4624",      1.610, 0.099, 0.092,   -363.0, 238.100, "08 05 29.038", "+46 24 51.78"),
    ("CD-69 2001",            1.616, 0.070, 0.068,  -1907.0, 332.610, "21 40 31.514", "-69 25 14.58"),
    ("HD 49995",              1.70,  0.23,  0.20,   -4034.0, 439.740, "06 50 20.810", "-18 37 30.58"),
    ("2MASS J0621-0101",      1.71,  0.46,  0.39,   -3206.0, 428.800, "06 21 34.807", "-01 01 55.01"),
    ("LSPM J2146+3813",       1.8557, 0.0048, 0.0048,  84.59, 22.9858, "21 46 22.285", "+38 13 03.12"),
    ("2MASS J0455+1144",      1.94,  0.16,  0.15,    1702.0, 349.500, "04 55 21.427", "+11 44 41.25"),
    ("2MASS J0734-0637",      1.950, 0.021, 0.021,   -554.6, 130.660, "07 34 39.097", "-06 37 12.21"),
    ("2MASS J1151-0313",      1.98,  0.20,  0.18,    1017.0, 125.880, "11 51 37.434", "-03 13 45.24"),
    ("UCAC4,076-006432",      2.042, 0.034, 0.033,   -893.8, 212.410, "06 34 29.385", "-74 49 47.12"),
    ("2MASS J0120+4739",      2.25,  0.17,  0.15,     473.0, 237.560, "01 20 04.561", "+47 39 46.56"),
    ("TYC 6760-1510-1",       2.46,  0.19,  0.18,   -1708.0, 102.890, "15 00 09.536", "-29 05 27.67"),
    ("UCAC2 15719371",        2.46,  0.10,  0.10,   -4282.0, 280.800, "09 44 09.884", "-37 45 31.09"),
    ("TYC 1662-1962-1",       2.637, 0.055, 0.054,  -1536.6, 286.510, "21 14 32.911", "+21 53 32.76"),
    ("HD 179939",             2.65,  0.17,  0.17,    3020.0, 334.320, "19 14 10.043", "+07 45 50.72"),
    ("BD-21 1529",            2.701, 0.059, 0.058,  -1660.1, 368.480, "06 37 48.004", "-21 22 21.94"),
    ("2MASS J1310-1307",      2.79,  0.59,  0.47,   -1520.0, 433.000, "13 10 30.804", "-13 07 33.55"),
    ("UPM J1121-5549",        2.803, 0.020, 0.020,   -282.5,  72.498, "11 21 18.136", "-55 49 17.77"),
    ("UCAC4,464-006057",      2.812, 0.052, 0.051,    932.0, 101.570, "04 09 02.050", "+02 45 38.32"),
    ("UCAC4,213-008644",      2.91,  0.13,  0.12,    -306.0,  80.987, "06 21 54.714", "-47 25 31.33"),
    ("GJ 3649",               3.016, 0.024, 0.024,   -520.4,  54.435, "11 12 38.97",  "+18 56 05.4"),
    ("Ross 248",              3.0446, 0.0077, 0.0077,  38.500, 10.3057, "23 41 54.99", "+44 10 40.8"),
    ("2MASS J1921-1244",      3.08,  0.21,  0.19,   -3490.0, 376.460, "19 21 58.124", "-12 43 58.61"),
    ("Proxima Centauri",      3.123, 0.015, 0.015,     28.65,  4.24646, "14 29 42.949", "-62 40 46.14"),
    ("TYC 9387-2515-1",       3.220, 0.081, 0.079,  -1509.1, 401.960, "06 18 54.643", "-80 19 16.54"),
    ("Alpha Centauri AB",     3.242, 0.060, 0.060,     29.63,  4.321, "14 39 36.495", "-60 50 02.31"),
    ("Gliese 445",            3.3400, 0.0051, 0.0051,  46.341, 17.1368, "11 47 41.377", "+78 41 28.18"),
    ("2MASS J1638-6355",      3.37,  0.29,  0.28,   -1428.0, 468.500, "16 38 21.759", "-63 55 13.16"),
    ("2MASS J0542+3217",      3.43,  0.75,  0.71,    5823.0, 884.600, "05 42 38.349", "+32 17 29.85"),
    ("2MASS J0625-2408",      3.700, 0.082, 0.080,  -1874.0, 534.880, "06 25 42.744", "-24 08 35.02"),
    ("Barnard's Star",        3.7682, 0.0031, 0.0031, 11.735,  5.96290, "17 57 48.498", "+04 41 36.25"),
    ("BD+05 1792",            3.965, 0.040, 0.040,   -962.7, 239.730, "07 48 07.037", "+05 27 22.51"),
    ("2MASS J2241-2759",      4.05,  0.16,  0.16,   -2810.0, 411.060, "22 41 50.996", "-27 59 47.04"),
    ("2MASS J1724-0522",      4.15,  0.26,  0.25,    3058.0, 489.500, "17 24 55.056", "-05 22 11.45"),
    ("StKM 1-554",            4.217, 0.036, 0.035,   -549.9, 151.970, "05 14 01.871", "+05 22 56.26"),
    ("GJ 3379",               4.227, 0.024, 0.024,   -157.43, 16.9861, "06 00 03.824", "+02 42 22.97"),
    ("2MASS J1936+3627",      4.23,  0.62,  0.57,    3830.0, 671.600, "19 36 57.294", "+36 27 57.71"),
    ("2MASS J0710+5228",      4.303, 0.039, 0.039,    507.6,  90.949, "07 10 52.167", "+52 28 18.49"),
    ("HD 146248",             4.341, 0.040, 0.039,  -1141.5, 334.870, "16 19 27.875", "-64 50 34.38"),
    ("2MASS J1724+0355",      4.37,  0.12,  0.12,    1991.0, 254.990, "17 24 34.633", "+03 55 26.75"),
    ("StKM 1-1456",           4.396, 0.043, 0.043,   1240.2, 144.934, "17 17 31.118", "+15 34 55.35"),
    ("Zeta Leporis",          4.43,  0.33,  0.30,    -878.0,  72.810, "05 46 57.341", "-14 49 19.02"),
    ("Lalande 21185",         4.6807, 0.0055, 0.0055, 21.973,  8.30437, "11 03 20.194", "+35 58 11.55"),
    ("HD 68814",              4.724, 0.090, 0.089,  -2242.0, 259.850, "08 13 57.112", "-04 03 12.56"),
    ("2MASS J1941-4602",      4.814, 0.050, 0.049,   -456.5,  66.848, "19 41 53.18",  "-46 02 31.4"),
]


WIKI_RAW = ("https://en.wikipedia.org/w/index.php"
            "?title=List_of_nearest_stars&action=raw")
WIKI_CACHE = BASE / "data" / "raw" / "wikipedia_nearest_stars.wikitext"


def verify_source(refresh: bool = False) -> int:
    """Check REFERENCE above still matches the page it was transcribed from.

    A hardcoded reference table is a silent liability: the source is edited,
    the comparison keeps passing against a stale copy, and the agreement it
    reports slowly stops meaning anything. This re-reads the live wikitext and
    diffs it, so drift is something the repository can detect rather than
    something to be assumed away.
    """
    import urllib.request

    if refresh or not WIKI_CACHE.exists():
        req = urllib.request.Request(
            WIKI_RAW, headers={"User-Agent": "star-neighbors/1.0 (research)"})
        WIKI_CACHE.write_text(
            urllib.request.urlopen(req, timeout=120).read().decode("utf-8", "replace"))
        print(f"  fetched {WIKI_CACHE.relative_to(BASE)}")
    raw = WIKI_CACHE.read_text()

    i = raw.find("Distant future and past encounters")
    seg = raw[i:]
    tbl = seg[seg.find("{|"):seg.find("\n|}", seg.find("{|"))]

    # Whole rows can be commented out - the page currently hides WD 0810-353
    # for an implausible radial velocity - and a comment can also sit on the
    # row separator without hiding anything, so strip comments before splitting.
    omitted = [re.search(r"\|\s*\[?\[?([^\n|\]]+)", m.group(1))
               for m in re.finditer(r"<!--(.*?)-->", tbl, re.S) if "\n|" in m.group(1)]
    body = re.sub(r"<!--.*?-->", "", tbl, flags=re.S)

    def val(c):
        m = re.search(r"\{\{val\|([-−\d.]+)", c)
        return float(m.group(1).replace("−", "-")) if m else None

    def name_of(c):
        c = re.sub(r"\[\[[^|\]]*\|([^\]]*)\]\]", r"\1", c)
        c = re.sub(r"\[\[([^\]]*)\]\]", r"\1", c)
        return re.sub(r"<.*?>", "", c).replace("–", "-").strip()

    live = {}
    for blk in body.split("\n|-")[1:]:
        cells = blk.split("\n|")[1:]
        if len(cells) < 4:
            continue
        nm, d, ep, dn = name_of(cells[0]), val(cells[1]), val(cells[2]), val(cells[3])
        if nm and d is not None:
            live[nm] = (d, ep, dn)

    hard = {n.replace("–", "-"): (d, epoch, dnow)
            for (n, d, _p, _m, epoch, dnow, _ra, _de) in REFERENCE}
    # The page labels one row "Scholz's Star and companion brown dwarf".
    alias = {"Scholz's Star and companion brown dwarf": "Scholz's Star"}
    live = {alias.get(k, k): v for k, v in live.items()}

    drift = [(n, hard[n], live[n]) for n in set(hard) & set(live)
             if abs(hard[n][0] - live[n][0]) > 1e-6
             or abs((hard[n][1] or 0) - (live[n][1] or 0)) > 0.5
             or abs((hard[n][2] or 0) - (live[n][2] or 0)) > 0.02]

    print(f"\n  source check against {WIKI_RAW.split('?')[0]}")
    print(f"    rows on the page      {len(live)}")
    print(f"    rows in REFERENCE     {len(hard)}")
    print(f"    names matched         {len(set(hard) & set(live))}")
    print(f"    values drifted        {len(drift)}")
    for n in sorted(set(live) - set(hard)):
        print(f"      only on the page:  {n}")
    for n in sorted(set(hard) - set(live)):
        print(f"      only in REFERENCE: {n}")
    for n, h, l in drift:
        print(f"      {n}: hardcoded {h} vs live {l}")
    if omitted:
        print("    commented out on the page (excluded there, kept here): "
              + ", ".join(m.group(1).strip() for m in omitted if m))
    ok = not drift and not (set(hard) ^ set(live))
    print(f"    -> {'REFERENCE is current' if ok else 'REFERENCE HAS DRIFTED'}")
    return 0 if ok else 1


def hms(s: str) -> float:
    h, m, sec = (float(x) for x in s.split())
    return (h + m / 60 + sec / 3600) * 15.0


def dms(s: str) -> float:
    s = s.replace("−", "-")
    sign = -1.0 if s.lstrip().startswith("-") else 1.0
    d, m, sec = (abs(float(x)) for x in s.replace("+", "").replace("-", "").split())
    return sign * (d + m / 60 + sec / 3600)


def sep_arcsec(ra1, de1, ra2, de2) -> float:
    r1, d1, r2, d2 = map(math.radians, (ra1, de1, ra2, de2))
    c = (math.sin(d1) * math.sin(d2)
         + math.cos(d1) * math.cos(d2) * math.cos(r1 - r2))
    return math.degrees(math.acos(max(-1.0, min(1.0, c)))) * 3600.0


def fnum(r, k, d=float("nan")):
    try:
        return float(r[k])
    except (KeyError, TypeError, ValueError):
        return d


def main() -> None:
    with open(PROC / "orbits.csv") as fh:
        stars = list(csv.DictReader(fh))
    ras = [fnum(s, "ra") for s in stars]
    des = [fnum(s, "dec") for s in stars]
    # the same positions wound back to J2000, for reference rows quoted there
    DT = -16.0                                   # Gaia epoch 2016.0 -> 2000.0
    ras0, des0 = [], []
    for s, r, d in zip(stars, ras, des):
        pmr, pmd = fnum(s, "pmra", 0.0), fnum(s, "pmdec", 0.0)
        cd = max(math.cos(math.radians(d)), 1e-6)
        ras0.append(r + (pmr * DT / 3.6e6) / cd)
        des0.append(d + pmd * DT / 3.6e6)

    rows, matched, missing = [], 0, []
    for (name, dmin, ep, em, epoch, dnow, ra_s, dec_s) in REFERENCE:
        ra, dec = hms(ra_s), dms(dec_s)
        best, best_sep, best_epoch = -1, 1e9, ""
        cd = math.cos(math.radians(dec))
        for label, rr, dd in (("2016", ras, des), ("J2000", ras0, des0)):
            for i, (r, d) in enumerate(zip(rr, dd)):
                if abs(d - dec) > 0.05 or abs((r - ra) * cd) > 0.05:
                    continue                  # cheap box before the exact sep
                sep = sep_arcsec(ra, dec, r, d)
                if sep < best_sep:
                    best, best_sep, best_epoch = i, sep, label
        # guard against grabbing a neighbor: present distances must agree
        if best >= 0:
            ours_now_ly = fnum(stars[best], "dist_now_pc") * LY_PER_PC
            if dnow > 0 and abs(ours_now_ly - dnow) / dnow > 0.10:
                best, best_sep = -1, 1e9
        if best < 0 or best_sep > MATCH_ARCSEC:
            missing.append(name)
            continue
        matched += 1
        s = stars[best]
        ours_d = fnum(s, "mc_dmin_med") * LY_PER_PC
        ours_lo = fnum(s, "mc_dmin_lo") * LY_PER_PC
        ours_hi = fnum(s, "mc_dmin_hi") * LY_PER_PC
        ours_t = fnum(s, "mc_tmin_med") * 1000.0
        ours_now = fnum(s, "dist_now_pc") * LY_PER_PC
        # does either interval reach the other's point estimate?
        agree = (ours_lo - 1e-9 <= dmin <= ours_hi + 1e-9) or \
                (dmin - em <= ours_d <= dmin + ep)
        rows.append({
            "star": name, "sep_arcsec": round(best_sep, 2), "matched_at": best_epoch,
            "ours_name": s["name"], "grade": s["grade"],
            "ref_dmin_ly": dmin, "our_dmin_ly": round(ours_d, 4),
            "our_lo_ly": round(ours_lo, 4), "our_hi_ly": round(ours_hi, 4),
            "d_ratio": round(ours_d / dmin, 3) if dmin else None,
            "ref_epoch_kyr": epoch, "our_epoch_kyr": round(ours_t, 1),
            "ref_now_ly": dnow, "our_now_ly": round(ours_now, 3),
            "intervals_overlap": agree,
            "sys_shift_ly": round(fnum(s, "sys_dmin_pc") * LY_PER_PC, 3),
        })

    print(f"  reference list: {len(REFERENCE)} stars within 5 ly")
    print(f"  matched in this catalog: {matched}   not present: {len(missing)}")
    if missing:
        print(f"    absent: {', '.join(missing)}")

    print(f"\n  {'star':22s} {'ref d_min':>10s} {'ours':>9s} {'ratio':>6s} "
          f"{'ref epoch':>10s} {'ours':>9s}  ok  grade")
    n_ok = n_epoch_ok = 0
    for r in sorted(rows, key=lambda x: x["ref_dmin_ly"]):
        dt = abs(r["our_epoch_kyr"] - r["ref_epoch_kyr"])
        ep_ok = dt <= max(30.0, 0.05 * abs(r["ref_epoch_kyr"]))
        n_ok += r["intervals_overlap"]
        n_epoch_ok += ep_ok
        print(f"  {r['star'][:22]:22s} {r['ref_dmin_ly']:9.3f} {r['our_dmin_ly']:9.3f} "
              f"{r['d_ratio']:6.2f} {r['ref_epoch_kyr']:10.1f} {r['our_epoch_kyr']:9.1f}"
              f"  {'..' if r['intervals_overlap'] else 'XX'}"
              f"{'' if ep_ok else '  epoch!'}  {r['grade']}")

    if rows:
        ratios = sorted(r["d_ratio"] for r in rows)
        med = ratios[len(ratios) // 2]
        print(f"\n  distance: {n_ok}/{len(rows)} agree within the stated intervals; "
              f"median ours/reference = {med:.3f}")
        print(f"  epoch:    {n_epoch_ok}/{len(rows)} agree to 5% or 30 kyr")
        worst = sorted(rows, key=lambda r: -abs(math.log(r["d_ratio"] or 1)))[:5]
        print("\n  largest disagreements:")
        for r in worst:
            print(f"    {r['star'][:24]:24s} ref {r['ref_dmin_ly']:6.3f} ly   "
                  f"ours {r['our_dmin_ly']:6.3f} [{r['our_lo_ly']:.3f},{r['our_hi_ly']:.3f}] ly"
                  f"   x{r['d_ratio']:.2f}   Galaxy-model +/-{r['sys_shift_ly']:.2f} ly")

    if "--csv" in sys.argv:
        out = pathlib.Path(sys.argv[sys.argv.index("--csv") + 1])
        with open(out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        print(f"\n  -> {out}")


if __name__ == "__main__":
    if "--verify-source" in sys.argv:
        raise SystemExit(verify_source("--refresh" in sys.argv))
    main()
