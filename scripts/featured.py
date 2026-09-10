"""Curated list of stars that must appear in the final catalog.

Some of these fail the automatic cuts - Beta Centauri is 120 pc away and never
comes close; Scholz's star matters only because of a past flyby - but they are
the objects people actually look for, so they are pinned in by name.

Resolved through SIMBAD to Gaia DR3 / HIP identifiers by 015_featured.py.
"""

# name -> short label used in the UI
FEATURED: dict[str, str] = {
    # --- explicitly requested ---
    "NAME Proxima Centauri": "Proxima Centauri",
    "* alf Cen A": "Alpha Centauri A",
    "* alf Cen B": "Alpha Centauri B",
    "* bet Cen": "Beta Centauri",
    "Ross 248": "Ross 248",
    "Ross 128": "Ross 128",
    "Ross 154": "Ross 154",
    "NAME Sirius": "Sirius",
    "NAME Barnard Star": "Barnard's Star",
    "GJ 445": "Gliese 445",
    "LAL 21185": "Lalande 21185",
    "GJ 710": "Gliese 710",
    "EGGR 290": "EGGR 290",
    "LSPM J2146+3813": "LSPM J2146+3813",
    "UPM J0812-3529": "UPM J0812-3529",

    # --- the rest of the current top 20 nearest ---
    "GJ 406": "Wolf 359",
    "GJ 65": "Luyten 726-8",
    "* alf CMa B": "Sirius B",
    "GJ 280": "Procyon",
    "* eps Eri": "Epsilon Eridani",
    "GJ 887": "Lacaille 9352",
    "GJ 411": "Lalande 21185",
    "GJ 866": "EZ Aquarii",
    "* 61 Cyg A": "61 Cygni A",
    "* 61 Cyg B": "61 Cygni B",
    "GJ 725 A": "Struve 2398 A",
    "GJ 725 B": "Struve 2398 B",
    "GJ 15 A": "Groombridge 34 A",
    "GJ 15 B": "Groombridge 34 B",
    "* eps Ind": "Epsilon Indi",
    "GJ 1111": "DX Cancri",
    "* tau Cet": "Tau Ceti",
    "GJ 1061": "GJ 1061",
    "NAME Teegarden's Star": "Teegarden's Star",
    "GJ 273": "Luyten's Star",
    "NAME Kapteyn's Star": "Kapteyn's Star",
    "NAME Luhman 16": "Luhman 16",
    "WISE J085510.83-071442.5": "WISE 0855-0714",
    "TRAPPIST-1": "TRAPPIST-1",
    "NAME van Maanen's Star": "Van Maanen's Star",
    "GJ 674": "GJ 674",
    "GJ 876": "GJ 876",
    "GJ 832": "GJ 832",
    "GJ 682": "GJ 682",

    # --- notable past encounters ---
    "NAME Scholz's Star": "Scholz's Star",
    "HD 7977": "HD 7977",
    "* alf Aql": "Altair",
    "* alf Lyr": "Vega",
    "* alf PsA": "Fomalhaut",
    "* alf Boo": "Arcturus",
    "* alf Aur": "Capella",
    "* bet Gem": "Pollux",
    "* alf Tau": "Aldebaran",
    "GJ 208": "Gliese 208",
    "GJ 451": "Groombridge 1830",
}
