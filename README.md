# The Sun's Changing Neighborhood

Which stars are our neighbors, and which ones will be? This project builds the
dataset and an interactive visualization of how the stars around the Sun move
toward and away from us across **±5 million years**.

The two widely-shared reference charts this started from plot about a dozen
stars over ±80,000 years using straight-line motion (they are not
redistributed here — they are not mine to republish). This extends that by
roughly 60× in time and 1000× in star count, and replaces the straight lines with orbits integrated
through the gravitational field of the Milky Way.

```bash
python3 scripts/serve.py
```

> **Status: exploratory.** The dynamics and the uncertainty propagation are
> consistent with each other and the numbers are reproducible, but the intervals
> are *measurement-sensitivity ranges*, not full confidence intervals: Gaia's
> astrometric covariance is propagated, the parallax zero-point is not, and
> uncertainty in the Galactic potential is measured separately rather than
> folded in — for a tenth of the candidates it exceeds the quoted interval.
> Treat individual predictions as indicative, and check both the interval and
> the model shift before quoting one. See [Known limitations](#known-limitations).
>
> Validated against the refereed Bailer-Jones encounter catalogs on Gaia
> source_id: **1740 paired encounters, median ratio 1.000**, and 53 same-release
> encounters whose epochs agree 53 for 53. The catalog also reproduces the
> published **encounter rate** — [20.7 ± 0.9 ± 2.4 against 19.7 ± 2.2 per Myr
> within 1 pc](#encounter-rate), a 0.3σ agreement reached independently.

Then open http://localhost:8777. (`python3 -m http.server --directory web` also
works, but it lets the browser cache `app.js`, so edits can appear to do
nothing; `serve.py` sends no-store.)

---

## What is in the dataset

**26,567 stars** with complete six-dimensional astrometry (position, distance,
proper motion, radial velocity), each integrated over ±5 Myr relative to the
Sun. **12,218** of them are close enough to be worth displaying and ship with
the interface; the rest are the screening margin and live in the CSV.

| | |
|---|---|
| Stars screened and integrated | 26,567 |
| …shipped to the interface | 12,218 |
| Within 10 pc today / 30 pc today | 289 / 7,108 |
| Matched to a SIMBAD object | 9,571 |
| Candidate co-moving systems | 649, covering 1,321 shipped stars |
| Quality grade A / B / C | 15,970 / 6,323 / 4,274 |
| Time span | ±5 Myr, integrated at 2,000 yr steps |
| Uncertainty | 256 draws per star, each integrated from Gaia's error ellipsoid (59 min to build) |

Selection happens in two stages, deliberately.

**Screening** pulls every star whose *straight-line* closest approach is under
**20 pc** within ±6 Myr, plus everything within 30 pc today. Twenty parsecs is
twice any distance one would call an encounter, and the margin is the point:
Galactic deflection moves a multi-Myr perihelion by of order a parsec, and in
one measured case by 4 pc, so screening at the same 10 pc used for scientific
interest would discard real close passes before they were ever computed.

**Interest** is decided afterwards, from the integrated orbit and the resampled
distribution rather than the straight-line value used to find the star. The web
interface ships those within 30 pc today or approaching within 10 pc; the full
screened catalog stays in `data/star_encounters.csv`.

## Where the data comes from

| Source | Used for |
|---|---|
| **Gaia DR3** | Astrometry for 26,409 stars, radial velocities for 25,503 |
| **Hipparcos-2** (van Leeuwen 2007) | 157 bright stars Gaia cannot measure — Sirius, α/β Centauri and similar saturate its detectors and are absent from DR3 entirely |
| **Published, by hand** | The parallax of α Centauri A and B (Akeson et al. 2021, 750.81 ± 0.38 mas). Two rows, listed in `PUBLISHED_PARALLAX`, because Hipparcos gives the two components of a 23-AU binary parallaxes that disagree by 5.6% |
| **SIMBAD** | Literature radial velocities for 1,064 stars (938 component values, 126 barycentric system values) — most nearby M dwarfs are too faint and red for Gaia's RVS — plus names and aliases, and the **full astrometry for 1 star**: Scholz's star (WISE 0720−0846) has DR3 photometry but no astrometric solution, being a very red L/T binary, so its discovery-paper parallax and proper motion are used instead |

Rather than downloading Gaia's 33-million-star radial-velocity catalog, the
closest-approach solution is written directly into the ADQL query, so ESA's
servers evaluate the filter and return only the 24,006 candidates that pass.

The three catalogs are deduplicated using the official Gaia–Hipparcos
cross-match **and** a positional match, because the official table covers only
~99.5k of ~118k Hipparcos entries and the gaps cluster among exactly the bright
stars needed. Positions are propagated between catalog epochs before matching:
Hipparcos is epoch 1991.25 and Gaia 2016.0, and ε Eridani moves 24 arcsec in
between, so naive matching misses every nearby star.

## The dynamical model

Both the Sun and every star are integrated as test particles through a Milky Way
potential — Miyamoto–Nagai disc, Hernquist bulge and nucleus, NFW halo (gala's
`MilkyWayPotential`, Price-Whelan 2017) — using RK4 at 2,000-year steps.

Solar position and motion: R₀ = 8.122 kpc (GRAVITY 2018), z☉ = 20.8 pc (Bennett
& Bovy 2019), V_c = 229 km/s (Eilers 2019), peculiar velocity (11.1, 12.24,
7.25) km/s (Schönrich, Binney & Dehnen 2010).

The closest approach is **refined between integration steps**, not taken as the
smallest sampled separation. The catalog contains genuinely fast stars — one
hypervelocity star moves at 806 km/s, covering 1.65 pc per 2000-year step — and
for those, picking the smallest sampled distance overstated the perihelion
badly: 1.583 pc against a converged 1.380 pc. Over a single step the relative
motion is effectively straight, which makes the squared separation an exact
parabola in time, so interpolating the three samples that bracket the minimum
recovers 1.379658 pc — matching a 20-year-step integration to six decimals, at
no extra cost.

### Does the potential actually matter?

Yes, and increasingly so with time. Comparing the integrated closest approach
against the straight-line model, for stars whose perihelion falls inside the
window:

| Encounter epoch | N | Median error of the straight-line model |
|---|---|---|
| within 0.1 Myr | 79 | 0.0000035 pc |
| 0.25 – 1 Myr | 390 | 0.0041 pc |
| 1 – 2.5 Myr | 558 | 0.070 pc |
| 2.5 – 4.5 Myr | 466 | **0.530 pc** |

At the ±80 kyr of the reference charts a straight line is perfectly adequate.
At 5 Myr the error exceeds the entire closest-approach distance of Gliese 710.

## Data quality: what a naive pipeline gets wrong

Closest-approach rankings are adversarial. A large inward radial velocity is
exactly what the filter rewards, so every bad measurement in the catalog
surfaces at the top of the results. Six failure modes were found and handled;
the last one was mine.

**1. Gaia's white-dwarf radial velocities are junk.** Gaia's RVS pipeline
cross-correlates spectra against normal-star templates, which fails completely
on a pressure-broadened degenerate atmosphere. Within 50 pc, only 3 white dwarfs
have any DR3 radial velocity, and their median |RV| is 374 km/s at S/N ≈ 6 —
against 17.7 km/s at S/N ≈ 367 for 1,974 comparable main-sequence stars.

This matters directly: **EGGR 290** and **UPM J0812-3529**, both requested for
this project and both circulating as dramatic close-encounter predictions, are
white dwarfs whose entire encounter rests on such a value. UPM J0812-3529 is
cataloged at −374 km/s radially while moving only 3.8 km/s transversely — a
combination that is geometrically implausible rather than merely unusual. Both
are retained but graded C and flagged.

**2. SIMBAD's radial velocity field mixes velocities and redshifts.** Entries
with `rvz_type = 'z'` are redshift-derived and reach 287,000 km/s for stars 18 pc
away. Filtered on `rvz_type = 'v'`, then bounded at 500 km/s, which also removes
genuinely wrong values such as EZ Aquarii at 6,824 km/s (an M dwarf 3.4 pc away)
and HD 163318 at −830 km/s — the latter had otherwise produced a spurious
0.54 pc "encounter".

**3. Close binaries contaminate each other's velocities.** Catalogs list a
component's instantaneous velocity, which in a tight pair includes orbital
motion about the barycenter rather than the system's actual travel. α Centauri A
and B are cataloged 7.3 km/s apart in radial velocity — not because the system
is doing anything unusual, but because each measurement catches a different
phase of an 80-year orbit. Propagated separately they gave closest approaches
0.23 pc apart, which is impossible for one bound system.

Bound pairs are identified by the classic common-proper-motion test — close on
the sky, consistent parallax, shared proper motion — not by 3-D separation.
Grouping in 3-D fails on precisely the pairs that matter: Hipparcos gives α Cen
A and B parallaxes disagreeing by 5.6%, placing them 0.07 pc apart in Cartesian
space although they orbit at 23 AU.

What actually travels through the Galaxy is the barycenter, and SIMBAD often
carries a system-level entry holding exactly that velocity: `* alf Cen` is
listed at −22.3 km/s against components at −15.25 and −22.59. So for the 404
systems whose components disagree by more than 2 km/s, the parent object is
looked up through SIMBAD's `h_link` hierarchy and its barycentric velocity
adopted where one exists — **126 components across the catalog**. Component
values are kept alongside in `rv_component`.

The effect on α Centauri is direct: component A's closest approach moves from
3.62 ly to **3.12 ly** against the reference charts' 2.97, and the A–B
disagreement falls from 0.232 pc to 0.075 pc.

**What remained was the parallax, and that is now fixed too.** α Cen saturates
Gaia and is absent from DR3, so it arrives from Hipparcos-2 — which gives A and
B parallaxes of **754.81 and 796.92 mas, disagreeing by 5.6%**. That put the two
components **0.228 ly apart in present distance**, for a pair that orbits at
23 AU: a true separation of 0.00036 ly, so the catalog was wrong about their
relative distance by a factor of 630, and it fed straight into their closest
approaches. Akeson et al. (2021) measured the system to **750.81 ± 0.38 mas** by
millimeter astrometry, which is both far more precise and — being one
measurement of one system — self-consistent between components. Adopting it:

| | present distance | closest approach |
|---|---|---|
| α Cen A, Hipparcos | 4.3210 ly | 3.122 ly |
| α Cen B, Hipparcos | 4.0927 ly | 2.876 ly |
| **both, published parallax** | **4.3441 ly** | **3.144 / 3.145 ly** |

The published present distance is 4.3441 ly, which this now reproduces exactly,
and the two components agree with each other to 0.001 ly instead of 0.246. The
override lives in `PUBLISHED_PARALLAX` in `02_merge.py` — two entries, sourced,
with the original catalog value preserved in `parallax_catalog` and the
citation in `parallax_source`, because a hand-set parallax is exactly the kind
of thing that should be auditable rather than buried.

**4. A quality grade is not an uncertainty.** Astrometry can be excellent while
the extrapolation from it is not. Every star's perihelion is therefore resampled
within its published errors (256 Monte Carlo draws) and reported as a
distribution. **Every draw is integrated through the same Galactic potential as
the nominal orbit, and the distance and the epoch both come from that one draw
population.**

That last point was got wrong for several revisions and is worth spelling out.
An earlier version solved the draws with straight-line motion, justified by
stage 4's measurement that the potential shifts a perihelion by under 0.005 pc —
which is true for encounters a few tens of kyr away, and was verified on exactly
those stars. The same stage also measured a median shift of **0.53 pc at 2.5–4.5
Myr**, and the closest-approach rankings are dominated by multi-Myr candidates.
The interface then ranked on that straight-line distance, labeled it "closest
approach", and paired it with the epoch from the *integrated* orbit: three
models in one row. Integrating the draws changes real answers —

| Star | straight-line draws | integrated draws |
|---|---:|---:|
| HD 7977 | 0.152 pc | **0.037 pc** |
| TYC 5584-552-1 | 0.899 pc | **1.622 pc** |
| Gaia DR3 6913732624445112832 | 3.356 pc | **0.935 pc** |
| 2MASS J05250205+0137210 | 5.702 pc | **1.622 pc** |
| Proxima, Ross 248, Scholz's | unchanged to 3 dp | |

— and moves them toward the literature: HD 7977's integrated range reaches the
published DR3 value, while the straight-line estimate was off by a factor of
four.
`scripts/mc_integrated_check.py` reproduces this on demand, computing **both**
models from one shared set of draws. (It used to read the shipped column and
call it "straight-line", which silently became a comparison of the integrated
model against itself the moment stage 5b changed. Its medians differ from the
shipped ones in the third decimal because it draws its own sample of 256.)

The width scales sharply with present distance, for stars whose closest approach
is under 2 pc:

| Distance today | N | Median 68% interval width | Median predicted approach |
|---|---|---|---|
| < 25 pc | 42 | 0.025 pc | 1.32 pc |
| 25 – 100 pc | 102 | 0.103 pc | 1.47 pc |
| > 100 pc | 166 | **0.356 pc** | 1.46 pc |

**4a. Gaia's astrometric errors are correlated, and the draws now say so.**
Parallax and both proper-motion components come out of one astrometric solution,
so their errors are not independent: across this catalog the median |r| is
**0.15** and the largest is **0.91**. Drawing them independently samples an
axis-aligned box instead of the real error ellipsoid. Measured on 400 candidates
within 5 pc, both ways, that is a second-order but not negligible effect — the
median interval width changes by 0.5%, but **46 of 400 widen by more than 10%**
and 35 narrow by more than 10%. The draws now use Gaia's published correlations,
factorised per star, falling back to independence for the Hipparcos and SIMBAD
rows, which publish none. `scripts/test_draws.py` asserts the properties that
matter: zero correlation reproduces the old independent draws exactly, a
prescribed correlation comes back out of the sample, and an inconsistent triple
degrades to the diagonal rather than raising.

**4b. For distant candidates the Galaxy model matters more than Gaia does.**
<a id="galaxy-model"></a>
The interval above answers "how much does this move when the measurements are
resampled". It does not answer "how much does it move if the Milky Way is not
quite the mass model we assumed", and that turns out to be the larger question
for exactly the stars that dominate the rankings. `scripts/09_sensitivity.py`
re-integrates the whole catalog under 21 perturbations — disc mass, scale
length and scale height, halo mass and scale radius, bulge mass, R₀, z☉, the
circular speed and the solar peculiar motion — and records the largest shift any
one of them produces:

| Distance today | N approaching within 5 pc | Median Galaxy-model shift | Median parallax zero-point shift | Median 68% interval | Model shift larger |
|---|---|---|---|---|---|
| < 25 pc | 271 | 0.00008 pc | 0.004 pc | 0.042 pc | **0** |
| 25 – 100 pc | 594 | 0.008 pc | 0.018 pc | 0.253 pc | 27 |
| > 100 pc | 819 | 0.082 pc | 0.053 pc | 0.562 pc | **138** |

So the nearby encounters — Proxima, Ross 248, Scholz's star, Gliese 710 — are
genuinely measurement-limited: the model moves them by less than a thousandth of
a parsec. The distant multi-Myr candidates are not. **165 of the 1,697 stars
approaching within 5 pc move further under a single plausible change to the disc
or halo than under the entire resampling of their astrometry**, and the second-
ranked future encounter is one of them. The parallax zero-point, over its
plausible range, does the same for a further 15:

| Star | closest approach | 68% measurement | model shift | ratio |
|---|---|---|---|---|
| Gliese 710 | 0.052 pc | 0.0063 | 0.0013 | 0.21 |
| HD 7977 | 0.038 pc | 0.0441 | 0.0206 | 0.47 |
| LAMOST J045521.42+114441.2 | 0.561 pc | 0.0491 | 0.0136 | 0.28 |
| Gaia DR3 1132123559268892288 | 0.322 pc | 0.0943 | **0.2637** | **2.80** |
| Gaia DR3 6913732624445112832 | 0.910 pc | 0.2441 | **0.5327** | **2.18** |

The dominant variant is almost always the **disc scale height**, which sets how
hard the disc pulls a trajectory back toward the midplane over millions of
years. Those stars carry a flag in the interface saying the Galaxy model matters
more than the measurements, and `sys_dmin_pc` is in the CSV for every star. It
is an envelope over chosen variants, not a posterior, and it should not be added
in quadrature to the measurement interval — the two answer different questions.

**5. Not every trajectory has a closest approach.** Some stars are still closing
at the edge of the ±5 Myr window. Their minimum separation is a *bound*, not an
event, and reporting it as a perihelion invents a fact. Where more than half a
star's draws reach their minimum at the boundary, the interface shows "≤ *d*"
instead of a closest-approach distance and epoch, and says which edge: at +5 Myr
the star is *still approaching, minimum lies later*; at −5 Myr it is *already
receding, minimum lies earlier*. Both were called "still approaching" until a
reviewer pointed out that this is simply false at the past boundary. The
fraction of draws that hit the edge is shown as a percentage, so a star with 51%
is not presented like one with 100%.

Seven stars in the catalog are censored in a majority of draws while their
median epoch lands *inside* the window, because their draws split between the
two boundaries: the astrometry cannot say whether the encounter is before
−5 Myr or after +5 Myr. One of them ships. Those are reported as "outside
±5 Myr · side unresolved" rather than being assigned a side, and a build
assertion checks that every censored star's epoch distribution actually reaches
a boundary.

**6. A finding of mine that did not survive its own fix.** Earlier versions of
this document reported that the perihelion computed from nominal astrometry was
a *biased estimator* — that for distant candidates it fell outside its own 68%
interval, and that rankings must therefore use the Monte Carlo median. The
evidence was real and the conclusion was wrong: the nominal value came from an
integrated orbit and the interval came from straight-line draws, so the two
disagreed because they were different models, not because the estimator was
biased.

With one coherent model, of the 1,697 stars approaching within 5 pc whose
perihelion is inside the window, **9 (0.5%) have their nominal value outside
the 68% interval**, and the median disagreement
between nominal and Monte Carlo median is **0.2%**. The star cited as proof
(Gaia DR3 1132123559268892288) reads nominal 0.319 pc against 0.322
[0.274, 0.368] — comfortably inside.

Rankings still use the Monte Carlo median, which is the right summary of a
distribution. The justification given for it was not.

## Validation

Against the reference charts supplied with the project, distances in light years:

| Star | now (ref / ours) | closest (ref / ours) | epoch kyr (ref / ours) |
|---|---|---|---|
| Gliese 445 | 17.60 / 17.14 | 3.45 / 3.34 | 45.0 / 44.0 |
| Ross 128 | 10.90 / 11.01 | 6.23 / 6.39 | 71.0 / 72.0 |
| Ross 248 | 10.30 / 10.31 | 3.02 / 3.05 | 36.0 / 36.0 |
| Ross 154 | 9.70 / 9.71 | 6.39 / 6.28 | 157.0 / 152.0 |
| Sirius | 8.60 / 8.60 | 7.86 / 8.17 | 64.0 / 45.7 |
| Lalande 21185 | 8.30 / 8.30 | 4.65 / 4.68 | 20.5 / 20.0 |
| Wolf 359 | 7.78 / 7.86 | 7.35 / 7.38 | −13.8 / −14.0 |
| Barnard's Star | 5.98 / 5.96 | 3.74 / 3.77 | 9.8 / 10.0 |
| α Centauri A | 4.36 / 4.32 | 2.97 / 3.12 | 28.4 / 27.8 |
| Proxima Centauri | 4.24 / 4.25 | 2.90 / 3.12 | 27.4 / 26.6 |

Present-day distances agree throughout, and closest approaches to within a few
percent.

### Against the published sub-5-light-year list

The strongest check available. `scripts/compare_published.py` compares this
catalog star by star against the compiled list of every star known to pass
within 5 ly of the Sun (Bailer-Jones et al. 2018/2022; Bobylev & Bajkova),
matching on sky position rather than name — the reference abbreviates 2MASS
designations, and these are by construction high-proper-motion stars, so
positions are propagated to a common epoch before matching.

| | |
|---|---|
| Reference stars | 48 |
| **Matched in this catalog** | **46** |
| **Median ratio, ours / reference** | **1.000** |
| **Epochs agreeing to 5% or 30 kyr** | **46 / 46** |
| Within the stated intervals | 27 / 46 |
| Disagreeing by more than 15% | 13 / 46 |

The epochs are effectively exact and the distances are unbiased. Where the two
disagree, they disagree for a reason this project already measured: **10 of the
13 discrepancies are stars whose Galaxy-model shift exceeds 0.15 ly**, all of
them presently 100–670 ly away with encounters 1.7–4.3 Myr out — precisely the
regime where [the mass model matters more than the astrometry](#galaxy-model).
The nearby, short-timescale encounters agree essentially perfectly: Gliese 445,
Ross 248, GJ 3379, UPM J1121-5549, LSPM J2146+3813 and Lalande 21185 all land
within 0.5%.

**Two corrections to what this file used to claim.** Both were resolved by
going to the refereed source instead of the compilation — see
[the Bailer-Jones cross-match](#bailer-jones) below.

*Gliese 710.* Earlier versions said this sat "about 20% low" against a published
0.0636 pc, then later that the 0.0636 pc figure was one I could no longer
source. It is sourceable: it is `dphmed` in Bailer-Jones (2022). Against *this*
compilation the agreement is still fine — reference 0.167 ± 0.012 ly, ours
0.170 ly, a 1.5% difference — but the disagreement with the refereed value is
real, and the cross-match below shows it is one no potential can produce.

*HD 7977.* Earlier versions called this "the one genuine, unexplained
disagreement", at 0.123 ly against the compilation's 0.478 +0.104/−0.078 ly.
It is explained, and the explanation is arithmetic. The compilation's
0.478 ly is **0.14656 pc**; this project's *straight-line* perihelion for
HD 7977 is **0.14659 pc**. They agree to 0.02%. The compilation entry is a
linear-motion value, and HD 7977 happens to be the most Galactic-deflection-
sensitive star in the whole sub-5-ly list: integrating it moves the perihelion
by **85%**, from 0.147 pc to 0.022 pc. Against the refereed integrated value
(Bailer-Jones 2022, 0.0641 pc with a 5–95% range of 0.019–0.117 pc) this
project's 0.0378 pc [0.0193, 0.0634] sits comfortably inside. HD 7977 was never
an outlier; it was a linear number being compared to an integrated one.

### Against the Bailer-Jones encounter catalogs
<a id="bailer-jones"></a>

The comparison above uses a hand-transcribed compilation of 48 stars.
`scripts/compare_bailer_jones.py` goes to the refereed catalogs that
compilation draws on and joins on **Gaia source_id**, so there is no
position-matching step to get wrong.

| | BJ2022 | BJ2018 |
|---|---|---|
| Reference | Bailer-Jones (2022), ApJL 935, L9 | Bailer-Jones et al. (2018), A&A 616, A37 |
| Gaia release | DR3 — same as ours | DR2 — two releases back |
| Reference rows | 61 (within 1 pc) | 3379 (within 10 pc) |
| Joined to this catalog | 55 | 1834 |
| Comparable (inside our ±5 Myr window) | 53 | 1740 |
| **Median ratio, ours / reference** | **1.014** | **1.000** |
| 16–84% of that ratio | 0.889 – 1.173 | 0.958 – 1.043 |
| Within 10% | 31 / 53 (59%) | 1464 / 1740 (84%) |
| **Epoch agrees to 5% or 30 kyr** | **53 / 53 (100%)** | 1585 / 1740 (91%) |

BJ2018 ids are translated DR2→DR3 through `gaiadr3.dr2_neighbourhood`, keeping
only unambiguous matches (2820 of 3379). A BJ2018 disagreement therefore mixes
method with two releases of astrometry, so that column is for reading the
spread, not for adjudicating individual stars. **1740 paired encounters with an
unbiased median ratio of 1.000 is the broadest validation this project has.**

**Epochs are exact and distances are not.** Every one of the 53 DR3 encounters
agrees on *when*, most to a fraction of a percent. That is the strong result: it
means the geometry, the coordinate transform, the potential and the integration
all agree with an independent group. What differs is *how close*, and only for
the long-baseline encounters.

**Is a disagreement even dynamically reachable?** This project ships both a
straight-line and an integrated perihelion for every star, so `|integrated −
straight line|` is the *entire* budget the Galactic potential has over the
window — the most gravity could possibly do is turning it off and back on. Of
the 53:

| | |
|---|---|
| Reference inside our interval + Galaxy-model envelope | 32 |
| Outside it, but within the gravity budget | 8 |
| **Beyond what any potential could produce** | **13** |

That last row is the interesting one, and **Gliese 710 is in it**. Bailer-Jones
publish 0.0636 pc [0.0595, 0.0678] alongside a parallax of 52.433 mas, a proper
motion of 0.428 mas/yr and an RV of −14.42 km/s — all of which match our inputs
to better than 0.1%, and an epoch matching ours to 0.3%. But those inputs give a
straight-line perihelion of **0.0512 pc**, and for this star gravity is worth
only **2.0%** over 1.29 Myr: switching the potential off entirely moves it from
0.05226 to 0.05129 pc, and no variant tested — `z_sun` 0–30 pc, `W_sun` ±3 km/s,
disc mass ×0.5–1.5, `V_circ` 200–250, `R₀` 7.5–8.5 kpc — moves it by more than
0.6%. Reaching 0.0636 pc requires 24%. We cannot reproduce their number from
the columns they publish next to it, and the same factor (~1.24) appears in
their 2018 table. This is stated as an unreconciled difference, not as an error
on their part: their method may use inputs the table does not expose.

What makes it an unreconciled difference rather than an error *here* is that the
published literature itself splits along the same line, and **this catalog is
not alone on its side of it**:

| study | data | perihelion | epoch |
| --- | --- | --- | --- |
| Berski & Dybczyński 2016 | TGAS | 0.0648 ± 0.0303 pc | 1.35 Myr |
| Bailer-Jones et al. 2018 | Gaia DR2 | 0.0676 ± 0.0157 pc | 1.281 Myr |
| de la Fuente Marcos ×2, 2020 | Gaia EDR3 | **0.051 ± 0.003 pc** | 1.29 ± 0.04 Myr |
| Bailer-Jones et al. 2022 | Gaia DR3 | 0.0636 [0.0595, 0.0678] pc | 1.292 Myr |
| Fernandez-Puig et al. 2026 | Gaia DR3 | 0.0621 ± 0.0023 pc | 1.3446 ± 0.0022 Myr |
| **this catalog** | Gaia DR3 | **0.0523 pc** [0.0491, 0.0554] | 1.294 Myr |

Two clusters, near 0.051–0.052 pc and near 0.062–0.068 pc, and the split does
not track the data release: DR3 appears on both sides. Our 0.0523 pc is within
1σ of the refereed EDR3 result of de la Fuente Marcos & de la Fuente Marcos
([RNAAS 4, 12, 2020](https://iopscience.iop.org/article/10.3847/2515-5172/abd18d)),
which is the value the *List of nearest stars* and *Gliese 710* Wikipedia
articles lead with. The most recent study, Fernandez-Puig et al.
([A&A, 2026](https://www.aanda.org/articles/aa/full_html/2026/06/aa59497-26/aa59497-26.html);
[arXiv:2605.16496](https://arxiv.org/html/2605.16496v1)), lands on the far side
at 0.0621 pc, so the disagreement is live rather than settled — which is why the
narrated tour quotes the spread instead of a single figure.

The integrator itself is not the suspect. Re-integrating Gliese 710 with
SciPy's DOP853 at `rtol=1e-12` reproduces the pipeline's RK4 result to five
decimal places (0.05226 pc, 1294.5 kyr), and for the 22 encounters inside
±1 Myr — where the potential has no time to act — the median ratio against
BJ2022 is **0.999**.

### Completeness within 20 light-years

The same source lists the known systems within 20 ly. Counting only the rows
that carry their own coordinates — one clean entry per system, which avoids the
spanning-cell trap that makes component rows inherit the wrong column — the page
lists **81 systems**, and this catalog contains **63 of them (78%)**. Their
present-day distances agree to the precision either side quotes:

| | |
|---|---|
| Median \|error\| | **0.00%** |
| 90th percentile | **0.00%** |
| Worst | 1.25%, Sirius (8.601 here against 8.709) |

Only two systems disagree by more than 0.3%, both close binaries whose parallax
depends on which component you ask: Sirius at 1.25% and Gliese 65 at 1.14%.
Everything else is exact.

The 18 absences are one story, not eighteen: **no radial velocity**. Twelve are
L, T or Y brown dwarfs — Luhman 16, WISE 0855, DENIS 0255, WISE 1541 and the
rest — far too faint and red for Gaia's RVS. The remaining six are M-dwarf
binaries and one white dwarf (EZ Aquarii, Kruger 60, Ross 614, Wolf 424,
Gliese 440, GJ 1005). This is the "2,324 stars excluded for want of a radial
velocity" limitation with names attached: **the missing neighbors are
systematically the coolest objects, and no ordinary single star within 20 ly is
missing for any other reason.**

**Two reference stars are absent.** 2MASS J0542+3217 encounters at **+5.82 Myr**,
outside the ±5 Myr window this project integrates — a correct exclusion, not a
miss. BD-21 1529 I could not locate in Gaia DR3 at the tabulated position: no
source within 72″ has a parallax consistent with 368 ly, and the brightest is
G = 14.7 against a quoted V = 9.67. Unexplained, and left that way.

**Sirius** remains at 7.86 / 8.17 ly and 64 / 45.7 kyr. This one cannot be
resolved with better data handling: `* alf CMa` is SIMBAD's top-level object,
so there is no barycentric entry to fall back on, and its cataloged −5.5 km/s
is a single epoch of a spectroscopic binary whose measured velocity varies over
a 50-year orbit. Inverting the encounter geometry, the reference chart's 64 kyr
implies a radial velocity near **−9.3 km/s** rather than −5.5. Since the
present-day distance matches exactly, the difference is entirely in the adopted
velocity. Sirius and every other cataloged spectroscopic binary now carries a
flag saying so, rather than being presented as a settled number.

**Scholz's Star** at 0.281 pc / −78.9 kyr sits within the published range
(0.25 – 0.30 pc at 70 – 80 kyr ago) without being fitted to it. Gliese 710 does
not reproduce its published value; see the reconciliation above. An earlier
version of this section claimed both did.

## The encounter rate
<a id="encounter-rate"></a>

Everything else here is a catalog: this star, that perihelion. This is a
*measurement* — how often any star passes within a given distance of the Sun,
which is the quantity the encounter literature reports because it sets the rate
of Oort cloud perturbation. `scripts/10_rate.py`.

**Result: 20.7 ± 0.9 (stat) ± 2.4 (sys) encounters per Myr within 1 pc**, against
a published 19.7 ± 2.2 (Bailer-Jones et al. 2018, from Gaia DR2). **0.3σ.**

### How it is measured

One physical assumption, and the data are made to test it rather than asked to
assume it: over ±5 Myr the solar neighborhood is statistically steady — the Sun
moves about 100 pc and the local density and kinematics do not measurably change
— so the true encounter rate is **constant in time**, and a histogram of
perihelion epochs must be flat.

It is not:

| \|t_ph\| Myr | encounters within 3 pc | rate / Myr | median present distance |
|---|---|---|---|
| 0–1 | 232 | 116 | 27 pc |
| 1–2 | 174 | 87 | 80 pc |
| 2–3 | 103 | 52 | 125 pc |
| 3–4 | 84 | 42 | 166 pc |
| 4–5 | 66 | 33 | 204 pc |

The fourth column is the cause. A star whose perihelion is far from the present
is far from the Sun *today*, and distant stars are the ones missing a radial
velocity and so missing from this catalog. **That turns the deficit from a
nuisance into an instrument**: near t = 0 the encountering stars are nearby and
the sample is nearly complete, so the rate there is very close to the truth.

The rate stops rising as the window shrinks — the plateau is the complete
regime:

| window \|t\| < | 0.125 | 0.1875 | **0.25** | 0.375 | 0.5 | 0.75 | 1.0 Myr |
|---|---|---|---|---|---|---|---|
| rate / Myr | 21.2 | 21.8 | **20.7** | 20.9 | 19.2 | 18.1 | 16.6 |

Flat to within Poisson error out to 0.375 Myr, falling steadily beyond. Inside
that window each encounter is weighted by 1/C(d_now), where C is this
catalog's completeness against the Gaia census — **measured against the census,
not against Gaia's radial-velocity column**, because this pipeline supplements
Gaia RVs from SIMBAD and using Gaia's own fraction would overcorrect the rate by
about 20%:

| present distance | census | Gaia has RV | this catalog | C |
|---|---|---|---|---|
| 0–5 pc | 50 | 68% | 44 | **88%** |
| 10–15 pc | 756 | 73% | 675 | **89%** |
| 20–25 pc | 2697 | 67% | 2025 | 75% |
| 25–30 pc | 4032 | 64% | 2809 | 70% |

Finally the counts are fitted to **N = k·d_ph²**, the scaling any homogeneous
flux of stars must obey. That fit is not decoration — if the counts did not
scale as the cross-sectional area, the sample would not be behaving like a
uniform flux and the rate would mean nothing.

### What moves the answer

| varied | range | rate / Myr |
|---|---|---|
| plateau window, 0.125–0.375 Myr | — | 20.7 – 21.8 |
| perihelion column (16th pct / median / 84th pct / nominal / straight-line) | — | 20.1 – 21.7 |
| **fit range, 1–2 pc up to 1–5 pc** | — | **20.7 – 25.5** |

Measurement uncertainty barely registers: swapping the median perihelion for the
16th or 84th percentile moves the rate by less than one encounter per Myr. The
**fit range dominates**, because residual incompleteness still grows with d_ph
after correction, which is why it is carried as the systematic.

### An independent check

Kinetic theory gives Γ = π d² n ⟨v⟩ with no reference to any perihelion at all.
With ⟨v⟩ = 48.4 km/s measured from the neighborhood sample, that is 12.4, 15.5
or 18.7 per Myr for a local density of 0.08, 0.10 or 0.12 stars/pc³. Consistent,
but it carries its own ~30% density uncertainty and the mean speed is biased low
by the radial-velocity sample favoring brighter, kinematically colder stars — so
it is quoted as a range and used to correct nothing.

### What this does not correct for

The correction is against the **Gaia census**, so stars missing from that census
— the faintest M dwarfs, brown dwarfs, anything lost beside a bright neighbor —
are not corrected for. **20.7 per Myr remains a lower bound.** The published
figure it agrees with has the same limitation, which is part of why they agree.

One number worth keeping in view. Applying this identical estimator to
Bailer-Jones's own DR2 table gives **3.1 ± 0.3 per Myr raw**, needing a factor of
**6.4** to reach their published 19.7. This catalog's DR3 raw plateau is
**16.9 ± 0.7**, needing **1.23**. That is what 33 million radial velocities buys
over 7 million — and it is why a rate can be quoted here while leaning on a
correction of 23% rather than one of 540%.

## Notable results

Closest **future** approaches, ranked by the median of the integrated draws.
Grade A only, and excluding trajectories that never reach a minimum inside the
window. Every candidate is listed, named or not — the unnamed Gaia sources are
the point, not an embarrassment.

| Star | closest approach | 68% range | when | distance today |
|---|---|---|---|---|
| Gliese 710 | **0.052 pc** (0.17 ly) | 0.05 – 0.06 | +1.30 Myr | 19 pc |
| Gaia DR3 1132123559268892288 † | 0.322 pc (1.05 ly) | 0.27 – 0.37 | +2.60 Myr | 207 pc |
| LAMOST J045521.42+114441.2 | 0.561 pc (1.83 ly) | 0.54 – 0.59 | +1.69 Myr | 107 pc |
| LSPM J2146+3813 | 0.569 pc (1.86 ly) | 0.57 – 0.57 | +83 kyr | 7 pc |
| 2MASS J14594591+8347085 † | 0.641 pc (2.09 ly) | 0.53 – 0.82 | +4.49 Myr | 170 pc |
| Gaia DR3 6913732624445112832 † | 0.910 pc (2.97 ly) | 0.79 – 1.04 | +2.95 Myr | 316 pc |
| Ross 248 | 0.934 pc (3.04 ly) | 0.93 – 0.94 | +36 kyr | 3.2 pc |
| HIP 94512 | 0.948 pc (3.09 ly) | 0.91 – 0.98 | +3.02 Myr | 103 pc |
| Proxima Centauri | 0.958 pc (3.12 ly) | 0.95 – 0.96 | +27 kyr | 1.3 pc |

† The Galactic mass model moves this closest approach further than resampling
its astrometry does, so the quoted range understates the real uncertainty — see
[the model-sensitivity section](#galaxy-model).

Closest **past** approaches, same basis:

| Star | closest approach | 68% range | when | distance today |
|---|---|---|---|---|
| 2MASS J06105543−4246035 | 0.210 pc (0.68 ly) | 0.20 – 0.22 | −1.16 Myr | 98 pc |
| 2MASS J22415098−2759470 | 0.582 pc (1.90 ly) | 0.56 – 0.60 | −2.78 Myr | 126 pc |
| 2MASS J06342938−7449471 | 0.624 pc (2.03 ly) | 0.62 – 0.63 | −896 kyr | 65 pc |
| TYC 1662-1962-1 | 0.820 pc (2.67 ly) | 0.81 – 0.83 | −1.54 Myr | 88 pc |
| SAO 172059 | 0.822 pc (2.68 ly) | 0.81 – 0.83 | −1.66 Myr | 113 pc |
| UPM J1121-5549 | 0.859 pc (2.80 ly) | 0.86 – 0.86 | −285 kyr | 22 pc |

Familiar stars, for scale: Proxima Centauri 0.958 pc at +27 kyr, α Centauri A/B
0.959 / 0.887 pc at +28 kyr, Gliese 445 1.024 pc at +44 kyr, Barnard's Star
1.155 pc at +10 kyr, Scholz's star 0.281 pc at −79 kyr, HD 7977 0.038 pc
[0.019 – 0.063] at −2.76 Myr. These tables are regenerated by
`scripts/05_validate.py` into `validation.json`, and a build assertion fails if
they stop matching the shipped column.

**Gliese 710 remains the most consequential encounter in the window** and one of
the best constrained, passing deep inside the Oort cloud. Ross 248 and Gliese
445 each displace Proxima Centauri as the Sun's nearest star, in roughly 33,000
and 42,000 years.

**990 trajectories never reach a minimum inside ±5 Myr** in more than half their
draws. Those are reported as bounds, and directionally: a star pinned to the
+5 Myr edge is "still approaching, minimum lies later", one pinned to −5 Myr is
"already receding, minimum lies earlier". Neither is an encounter.

## Pipeline

| Stage | Does |
|---|---|
| `scripts/01_fetch.py` | Queries Gaia, Hipparcos and SIMBAD into `data/raw/` |
| `scripts/015_featured.py` | Resolves the curated star list to catalog ids |
| `scripts/02_merge.py` | Deduplicates, resolves RVs, grades quality → `merged.csv` |
| `scripts/03_names.py` | Attaches SIMBAD names and aliases → `named.csv` |
| `scripts/04_orbits.py` | Integrates orbits, fits Chebyshev polynomials → `orbits.csv`, `cheb_coeffs.npy` |
| `scripts/055_uncertainty.py` | Per-star Monte Carlo, every draw integrated through the potential |
| `scripts/05_validate.py` | Structural + physical assertions, reference comparison → `validation.json` |
| `scripts/06_export.py` | Selects what to display and packs `web/data/` |
| `scripts/07_csv.py` | Writes `data/star_encounters.csv` and its column reference |
| `scripts/09_sensitivity.py` | Re-integrates under 21 model variants → `sys_dmin_pc`, `sensitivity.json` |
| `scripts/10_rate.py` | Completeness-corrected encounter rate → `encounter_rate.json`. `--scan` prints the systematics |
| `scripts/08_verify.py` | Checks the published artifacts agree → `build_report.json` |
| `scripts/astro.py` | Coordinate transforms, Galactic potential |
| `scripts/gates.py` | The completeness gate both exporters apply before publishing |
| `scripts/test_gates.py` | Tests for that gate — 12 injected defects, 7 valid inputs, plus the real catalog |
| `scripts/compare_published.py` | Star-by-star comparison against the published sub-5-ly encounter list. `--verify-source` re-reads the live page and fails if the hardcoded reference has drifted |
| `scripts/compare_bailer_jones.py` | Cross-match against the refereed Bailer-Jones encounter catalogs on Gaia source_id — 1740 paired encounters, plus a test of whether each disagreement is dynamically reachable |
| `scripts/systems.py` | Bound-system grouping and barycentric velocities |
| `scripts/featured.py` | Curated star list pinned into the catalog |
| `scripts/tap.py` | Minimal TAP client (stdlib only) |
| `scripts/serve.py` | No-cache dev server for `web/` |
| `scripts/fit_study.py` | Diagnostic: Chebyshev degree/precision vs perihelion error. Not part of the build |
| `scripts/mc_integrated_check.py` | Diagnostic: straight-line vs integrated Monte Carlo, both computed from the same draws. Not part of the build |

Run in order:

```bash
python3 scripts/01_fetch.py        # ~10 min, hits Gaia + SIMBAD
python3 scripts/015_featured.py
python3 scripts/02_merge.py
python3 scripts/03_names.py        # first run ~10 min, then cached
python3 scripts/04_orbits.py       # ~4 min
python3 scripts/055_uncertainty.py # ~80 min, integrates every draw
python3 scripts/09_sensitivity.py  # ~5 min, 21 Galaxy-model variants
python3 scripts/10_rate.py         # completeness-corrected encounter rate
python3 scripts/05_validate.py     # exits non-zero if an assertion fails
python3 scripts/06_export.py
python3 scripts/07_csv.py
python3 scripts/08_verify.py       # artifacts must agree with each other
python3 scripts/test_gates.py --catalog   # the gate itself still rejects
```

`pip install -r requirements.txt` first; the pipeline needs only numpy and
scipy.

Two gates, and both fail the build rather than printing plausible-looking
numbers. `05_validate.py` checks the science — required fields, finite values,
ordered percentiles, censor flags consistent with the censored fraction, no
object under two identities, and the published rankings derivable from the
shipped column — and records every assertion, with an overall status, in
`validation.json`. `08_verify.py` checks the *build*: that `orbits.csv`,
`star_encounters.csv`, `stars.json`, `cheb.bin` and `meta.json` describe the
same catalog, and that the embedded validation says pass. The exporters
refuse to run at all if any row is missing its Monte Carlo columns, carries a
blank, `nan` or `inf` where a number belongs, or is not `mc_model=integrated`,
so a premature export cannot publish blanks. `test_gates.py` runs nineteen
cases — twelve injected defects that must be caught and seven valid inputs
(zeros, negatives, scientific notation) that must not be — plus the real
catalog. A gate that quietly stopped rejecting would otherwise look identical
to one that was working, and a gate that started rejecting good data would fail
the build for nothing.

Every stage is re-runnable and SIMBAD name lookups are cached in
`data/processed/simbad_cache.json`; delete it to force a refresh. `data/raw/`
and `data/processed/` are gitignored as regenerable (~44 MB); `web/data/` is
the shipped artifact.

## How the visualization works

Storing sampled trajectories for 12,218 stars would run to tens of megabytes.
Instead each star's heliocentric offset is fitted with a degree-12 Chebyshev
polynomial in time, and the whole catalog compresses to 1.9 MB.

The figure that matters is not the positional residual but the error in the
recovered closest approach, since that is what the project reports. Measured
against the refined integration, using the exact float32 coefficients that ship:

| | error in closest approach |
|---|---|
| median | 6.3×10⁻⁷ pc |
| 99th percentile | 2.4×10⁻⁵ pc |
| worst case | 3.5×10⁻³ pc (724 AU) |
| worst among the 310 stars closing within 2 pc | 1.2×10⁻³ pc (243 AU) |

Raising the degree to 16, 20 or 24 changes none of these, and neither does
storing in float64 — the compression is not the limiting factor anywhere.

Because a polynomial is cheap to evaluate, the coefficients are uploaded as a
floating-point texture and **the trajectory is evaluated in the vertex shader**.
Every star moves on real integrated dynamics at frame rate, with no per-frame
CPU work: the browser sends 13 basis values per frame and the GPU does the rest.

The distance–time chart uses a symmetric-log time axis (linear within ±10 kyr,
logarithmic beyond) and a logarithmic distance axis. The reference charts' linear
axes cannot show a 26,000-year encounter and a 1.3-million-year one on the same
plot; this can. That axis is explained in words beside the chart rather than
left to be inferred from the tick spacing.

Curves and names are budgeted separately. Two dozen labeled trajectories
collapse into a smear around the present epoch, so about 28 curves are drawn for
context but only the closest six are named — plus whatever is selected or
hovered, which also gets its full designation, a brighter stroke and a vertical
bar marking the 16th–84th percentile of its closest approach. That bar is the
only uncertainty drawn rather than written, and it belongs to the star you are
actually looking at.

### From Earth: the sky, not the map

Everything above is a map seen from outside. **From Earth** puts the camera at
the Sun and draws each star at the brightness it would actually have in the
sky, `m = M + 5 log₁₀(d) − 5`, taking `d` from the same Chebyshev evaluation
that already places the star. It costs one extra attribute and a `log` in the
vertex shader; the point size stops being distance-driven and becomes
magnitude-driven, and the view-radius cull is replaced by a magnitude limit,
since on the sky a star is not hidden by being far away, only by being faint.

It is worth having because it answers a question the map cannot. Scrub time and
the brightest object in the sky changes completely:

| epoch | brightest | Gaia G | naked-eye stars in this catalog |
|---|---|---|---|
| −2.765 Myr | **HD 7977** | **−8.2** | 272 |
| −0.298 Myr | Aldebaran | −1.3 | 592 |
| today | Sirius | −1.1 | 657 |
| +1.295 Myr | **Gliese 710** | **−3.8** | 368 |
| +3.020 Myr | HIP 94512 | −3.0 | 234 |

**These are Gaia G magnitudes, and the column is labeled that way for a
reason.** G is what Gaia measures for every star in the catalog, so it is the
only band available here — but it is not the scale the famous numbers are
quoted on. For a Sun-like star the two agree to about a tenth of a magnitude;
for a red one they do not. Gliese 710 is a K7 dwarf: −3.8 in G is roughly −3.2
visually, and the published value is about −2.7, computed at a slightly more
distant perihelion than this catalog finds. Three defensible numbers for one
star, differing by filter and perihelion rather than by any disagreement about
how bright it is. Quoting the G figure against a remembered visual one — Sirius
at −1.46 — would overstate the ratio by a factor of three, so the interface and
the narration both name the band.

HD 7977 is a G0V star, where G and V agree to 0.15 mag, so its roughly −8 holds
on either scale: brighter than Venus's −4.92, bright enough to cast shadows —
with the caveat that its perihelion is the least certain in the close list, so
the figure could be two magnitudes either way. The brightest star is named in
the readout and the name is a button, because naming it and then making you
hunt the sky for it would be useless.

**This is not the night sky, and the interface says so.** The catalog holds
the solar neighborhood and everything screened as an encounter candidate — not
the distant giants that make up most of the constellations. About 660 of its
stars are naked-eye visible today against roughly 9,000 in the real sky, and
stars without Gaia photometry are not drawn at all. It answers "which of the
Sun's *neighbors* can you see, and when", which is a real question, rather
than pretending to be a planetarium.

### The narrated tour

A ten-minute documentary pass through the catalog, on a button in the title
panel. It runs in four movements: the immediate neighbors and how briefly they
hold that title; Gliese 710 and the closest approach anything will make while
our hardware could still be flying; the encounters already behind us; and the
rate that makes all of it ordinary.

**It is written for a reader leaving middle school**, which is a harder
constraint than it sounds, because the subject arrives pre-loaded with jargon.
Every unit is defined at the moment it is first spoken — a parsec as three and
a quarter light years, a red dwarf as a small cool star, a brown dwarf as an
object too small to shine — and the ones that could not be made to earn their
place were cut rather than explained. Proper motion was quoted in
milliarcseconds per year; it is now "watch any other star for a year and it
shifts a little against the ones behind it; Gliese 710 barely shifts at all."

The magnitude passage was the worst offender and needed a structural fix rather
than a reword. It explained the difference between Gaia's green band and the
visual scale and then quoted a figure from each, without ever mentioning the
one fact a listener actually needs: that the scale runs backwards. Both numbers
now come from the band the on-screen readout is already showing — Sirius at
about −1, Gliese 710 at nearly −4 — so the narration and the display agree,
only one scale has to be explained, and the comparison ("more than ten times
brighter") is the 11.6× the two magnitudes really imply.

Plain language costs length: the same content went from 1,105 words to 1,388,
and the button now says ten minutes rather than eight. That is the trade, and
it is the right way round — a listener who has to decode "the visual scale
Sirius's minus one and a half is quoted on" has already stopped listening.

Two further decisions shape the implementation.

**The narration is spoken by the browser, not shipped as audio.** Five minutes
of voice track is several megabytes, would have to be regenerated every time a
number changed, and would be one more artifact free to drift from the data. The
`speechSynthesis` API costs nothing and cannot go stale. The price is that voice
quality varies by platform — so **captions are always drawn** and the tour is
complete with the sound off.

Which voice gets used is not left to the browser. Its default is the OS system
voice, which on macOS is Samantha — a formant voice from 2009 that is genuinely
hard to sit through for five minutes. The narration is written for **UK English
Female** and asks for it by name, falling through an ordered list of the same
register on Edge and macOS. If none of them is installed it takes the best
British voice on the machine, preferring a neural one. **The novelty voices are
excluded outright** — macOS lists Bubbles, Zarvox and Bad News alongside the
real ones and several sort early alphabetically, so any "first English voice"
fallback can land on a joke. In Chrome this gets Google UK English Female; on a
stock Mac with no premium voices installed it gets Daniel.

Each beat is spoken **one sentence at a time** rather than as a paragraph. A
synthesiser handed five sentences at once flattens its intonation and races the
full stops; short utterances with a 260 ms gap between them sound like narration
instead of a screen reader. It also sidesteps Chrome silently truncating
utterances past about fifteen seconds.

**Beats advance when the speech ends, not on a timer.** Speech rate varies with
platform, voice and the user's own accessibility settings; any fixed schedule
would drift out of sync with the words inside a minute. Where there is no voice
to wait for, each beat falls back to an estimate from its own word count, and a
hard ceiling per beat means a synthesiser that never reports back stalls the
narration rather than the tour.

**Interstellar travel times get corrected out loud.** The familiar "Voyager
would take seventy thousand years to reach Proxima Centauri" is distance over
speed with the star held still, and the tour used to repeat it. Curvature is not
what spoils it: integrating a coasting probe through the Galactic potential
departs from a straight line by **0.0001% of the trip** over seventy-five
thousand years, because the Sun and a slow probe fall through the Galaxy
together. What spoils it is that the target moves. **Proxima recedes at 32 km/s
after perihelion, nearly twice Voyager 1's 17 km/s, so a probe launched today
never reaches it** — it falls about 1.1 light years short and loses ground from
then on. Gliese 710 closes at only 14 km/s, slower than the probe, so an
intercept does exist: integrated, about **594,000 years**, some 46% sooner than
the naive figure because the star does most of the traveling. The narration now
says the seventy-thousand-year number is arithmetic on a stationary target and
then says why that fails.

**The screen is cleared for the subject.** The detail card and the chart are
both hidden for the duration: the card is a wall of numbers beside the one
object you are being asked to look at, and it covers a third of the view, while
the chart is a second telling of the same story — and with it visible the dock
is tall enough that the caption lands on whatever the camera just centered. The
narration says the numbers that matter out loud.

**And the rest of the interface is inert while it runs.** Dimming the rail and
the orientation panel to 45% was not enough, because dimmed controls still
work. Pressing Today, dragging the scrubber, picking a camera preset or tapping
a list row all fight the shot the narration is talking over: the tour keeps
animating, the click lands on top of it, and the result reads as broken rather
than interrupted. The keyboard was worse — `J`, `L`, `K`, `E` and the zoom keys
fell straight through the tour's own handler and could start playback or flip
to the sky view mid-sentence.

`#title`, `#rail`, `#nav` and `#dock` now carry the `inert` attribute for the
duration, which removes them from the pointer, focus and accessibility trees in
one attribute and restores them in one line, so the two states cannot drift.
The keyboard handler returns after the tour's own keys rather than falling
through. Nothing is hidden: the epoch readout and the moving scrubber stay at
full contrast because watching time move is part of what the tour is showing,
while every control that is inert is also dimmed, so it looks unavailable
instead of silently ignoring a click. Only the caption stays live, which is
also what keeps **Exit Tour** reachable.

Verified by clicking the dimmed controls with real pointer events rather than
synthetic ones — `element.click()` bypasses hit-testing and fires the handler
anyway, which briefly looked like a bug in `inert` and was a bug in the test.

**The camera turns slowly throughout**, about two degrees a second. A held shot
of two points of light is a photograph: you cannot tell how far apart they are
or which is in front. Drifting around the pair gives the parallax that makes the
geometry read, and it is slow enough that the next beat's fly-to still looks
like a deliberate move rather than a correction. It stops entirely under
`prefers-reduced-motion`.

One thing worth stating because it caused two visible bugs: **a camera aimed at
a fixed point is aimed at the past.** `focusStar` frames a star where it is
*now*, so a shot that flies to Gliese 710 and then scrubs 1.3 million years ends
up pointing at empty space while the star arrives somewhere else. The subtler
version bit the Alpha Centauri beat: it inherited a tight framing from the
previous beat and then ran time forward 28,000 years, over which Alpha Centauri
travels 0.92 pc through a frame 0.46 pc wide — it left the top of the screen
while the narration was still describing it.

A third framing rule came out of the same beat: **the shot is composed for the
strip of canvas the caption is not covering, not for the canvas.** `focusStar`
centered the Sun–star pair in the middle of the viewport, which during the tour
is behind the text. Measured over the succession beat, Ross 248 was behind the
caption in **48 of 84 samples** and settled at y = 669 on a 930 px canvas whose
caption starts at 587.

The correction cannot be a constant, because the strip is a different size on
every layout — the same caption is three lines at 840 px wide and eight on a
phone. `tourFraming()` measures the gap between the bottom of the title bar and
the top of the caption each time a shot is set up, then moves the look-at point
down in z so the pair's center lands in the middle of that strip, and pulls the
camera back in proportion to how much of the height the strip has lost:

| layout | clear strip | lift | pull | Ross 248 settles at |
| --- | --- | --- | --- | --- |
| 840 × 930 | 61 – 633 | 0.127 | 1.35× | y = 507 (was 669) |
| 768 × 1024 | 61 – 747 | 0.105 | 1.30× | y = 586 |
| 375 × 812 | 61 – 451 | 0.185 | 1.47× | y = 384 |

Both ends of the strip matter. Centering on the caption alone over-lifted and
put the Sun one pixel above the title bar in the widest shot on a phone. All ten
star-focusing shots in the tour now keep **both** the subject and the Sun inside
the strip on all three layouts, and outside the tour the framing is unchanged to
five decimal places — `tourFraming()` returns a no-op when the caption is not up.

So the rule every beat now follows: **if a beat moves time, it either tracks a
subject or recenters on the Sun.** Nothing scrubs against a camera fixed on a
world point. Beats that travel to a star re-frame on arrival, and beats about a
*distance* — Alpha Centauri approaching — are framed on the Sun–star pair so the
gap itself is what shrinks on screen: 493 pixels down to 303 over that beat,
against a true closing of 27%.

The same rule has a second half: **if a beat names more than one subject, the
camera hands over with the words** — unless the sentence exists to dismiss them,
as the beat after this one does, where marking each star in turn would argue
against the narration while the shot pulls away. The succession beat is the case that forced
it — the title of nearest star passes from Alpha Centauri to Ross 248 to Gliese
445 and back again inside twenty thousand years, and holding on Ross 248 for the
whole beat says the opposite of what is being narrated. It now runs as five
moves on five sentences, each scrubbed to the epoch that sentence describes:

| epoch | nearest star | distance | held for |
| --- | --- | --- | --- |
| today – 33.0 kyr | Alpha Centauri (Proxima, then B from 31.8 kyr) | 4.246 → 3.172 ly | 33.0 kyr |
| 33.0 – 42.4 kyr | Ross 248 | min **3.045 ly** at 36.5 kyr | 9.4 kyr |
| 42.4 – 50.0 kyr | Gliese 445 | min **3.340 ly** at 44.3 kyr | 7.6 kyr |
| 50.0 – 65.4 kyr | Alpha Centauri B again | 3.949 → 5.126 ly, receding | 15.4 kyr |

Read off this catalog at 200-year resolution over the A- and B-grade stars.
Wikipedia's Ross 248 article gives "about 33,000 years" for the handover and
"about 9,000 years" for the reign, and 3.048 ly at 36,500 years for the minimum;
the *List of nearest stars* table gives Gliese 445 at 3.3400 ly. Alpha Centauri
is never worse than **third** anywhere in this window — it is briefly displaced
to third around 41.8 kyr, when Ross 248 and Gliese 445 hold the first two places
at once, which is why the narration claims the top three rather than second
place. Beyond it the churn continues: the title changes hands about **thirty
times** in the next million years.

### The color of a star

Every point is colored from its effective temperature through a blackbody
approximation, so an M dwarf is orange-red and an A star is blue-white. The
data always carried this — `teff_to_rgb` has been in the export from the start
— but for a long time almost none of it reached the screen. Four separate
things were swallowing it, and they had to be fixed together:

- **The premultiplied alpha overflowed.** The fragment shader emits
  `vec4(col * a, a)` with `a = core + halo`, which peaks at **1.35**. Above 1.0
  every channel of `col * a` clips on write, so the middle of *every* star was
  pure white no matter how carefully its color had been computed. Clamping `a`
  to 1.0 is the single change that made hue visible at all.
- **The hot core was additive.** `v_color + core * 0.55` added the same 0.55 to
  all three channels, and since `teff_to_rgb` pins red at 255 for anything
  cooler than 6600 K, that drove green and blue up to meet it. Ross 248 is
  `rgb(255,142,29)` in the data and was reaching the screen as
  `rgb(255,255,134)`. It now leans toward white with `mix` instead, cubed so
  only the middle ~30% of the radius reads as hot.
- **The selected star was tinted blue** — a 55% mix toward `#8cd9ff`, applied
  to the one star the interface is asking you to look at. The tour narrates
  Gliese 710 as "distinctly orange" over a point that rendered
  `rgb(220,243,250)`. The reticle, the label, a 2.2× size boost and full alpha
  already identify the selection, so the tint is gone.
- **The brightest stars had no temperature.** Gaia saturates on them, so they
  arrive via Hipparcos with neither `teff` nor BP−RP and fell to a flat 5000 K
  default: Sirius, Vega, Arcturus, Altair, Procyon, Capella, Pollux, Fomalhaut
  and both Alpha Centauri components all rendered as the same cream dot — that
  is, essentially every star a viewer would recognize well enough to check.
  `sptype_to_teff` now reads the MK type instead, on the Pecaut & Mamajek
  (2013) dwarf scale with a separate giant scale for luminosity class III.

The spectral-type fallback is checked against measured temperatures:

| star | type | derived | measured | error |
| --- | --- | --- | --- | --- |
| Vega | A0V | 9700 K | 9602 K | +1.0% |
| Sirius | A0mA1Va | 9700 K | 9940 K | −2.4% |
| Procyon | F5IV-V | 6560 K | 6530 K | +0.5% |
| Alpha Centauri A | G2V | 5792 K | 5790 K | +0.0% |
| Alpha Centauri B | K1V | 5137 K | 5260 K | −2.3% |
| Pollux | K0IIIb | 4750 K | 4666 K | +1.8% |
| Arcturus | K1.5III | 4608 K | 4286 K | +7.5% |

Worst case 7.5%, which for Arcturus is the difference between `rgb(255,219,191)`
and `rgb(255,213,179)` — smaller than the error it replaced, since the 5000 K
default put it at `rgb(255,228,205)`. It covers 93 of the 235 stars that were on
the default; the remaining 142 have no parseable spectral type either and keep
5000 K. Real stellar colors are pale, and nothing here saturates them beyond
what the temperature gives.

The point-size floor also went from 2.0 to 2.8. It is in device pixels, so on a
2× display 2.0 was under one CSS pixel — too few samples for a hue to survive.

### One filter that only half worked

The star list, the chart and the sky labels all ask `recomputeVisible()` which
stars pass the filters. The 3-D view did not: the vertex shader re-derived the
same decision from the flag bits it had been given. Two statements of one rule,
with nothing keeping them in step — and they drifted. `recomputeVisible()` also
applies the **Max spread** slider, and the shader never knew about it.

Measured: dragging Max spread to 0.05 pc removed **3,587 stars** from the list,
the chart and the labels, and changed the rendered image by **zero pixels**. The
control appeared to do nothing to the view it most obviously pointed at. The
contrast case proves the measurement rather than the code — *Notable stars only*
removed 10,771 stars and moved 3,592 pixels.

The fix is not to add the missing rule to the shader, which would leave two
statements of it again. Visibility is now decided once on the CPU and uploaded
as a per-star attribute, and the shader reads the answer:

```glsl
bool hidden = a_vis < 0.5;     // was: (u_onlyFeat == 1 && !featured) || …
```

`uploadVis()` is called from inside `recomputeVisible()` rather than beside its
callers, so a filter added later cannot reach the lists without also reaching
the view. The duplicated uniforms, the flag attribute and the renderer-private
grade bit all went with it.

Verified afterwards at the pixel level, on the star's own pixel rather than a
window around it — a neighbor 4 px away made an earlier check read as a false
pass. Grade C: 25 hidden, 733 shown. Max spread: 604 included, 25 filtered. A
selected star that the filters exclude is still drawn, which is what lets Focus
fly to something outside the current view.

### Drawing a close approach that happens far from now

The distance chart plots every curve against a **symlog** time axis, so the
10,000 years either side of the present get a tenth of the width rather than a
five-hundredth. Curves were sampled at 300 points spread evenly across that
axis — evenly in *screen* space, which is the right instinct and the wrong
result for anything that happens a long way from the present.

A megayear out, one of those 300 steps spans several thousand years, while the
approach itself is over in a few hundred. The perihelion is a narrow V and the
polyline stepped clean over it:

| star | dip width, in grid steps | true perihelion | as drawn |
| --- | --- | --- | --- |
| Gliese 710 | 0.23 | 0.170 ly | **1.242 ly** |
| HD 7977 | 0.03 | 0.071 ly | **2.694 ly** |
| Gaia DR3 1132123559268892288 | — | 1.040 ly | **9.180 ly** |

Nine of the 28 curves on a desktop-width chart were wrong by more than 50%, and
Gliese 710 — the encounter the whole project is about — never went below one
light year, while the perihelion dot sat correctly at 0.17 with no curve
anywhere near it.

Adaptive refinement would not have found these. Subdividing where neighboring
samples disagree is the usual fix, and it fails here precisely because the two
samples either side of the spike sit at the *same* height. Nothing in the
sampled data hints that a spike is between them.

The perihelion time is already in the catalog, so the curve samples it
directly: the uniform grid, plus a geometric ladder of offsets from one year to
100,000 years either side of that star's own closest approach. A ladder resolves
the V whatever its width, which matters because the widths span four orders of
magnitude. It costs 42 extra points per curve, and all 28 now reach their true
minimum exactly — worst ratio 1.0000. A full redraw takes 5 ms.

Two deliberate leftovers. The curve is the **nominal** trajectory while the dot
and whisker are the **Monte Carlo median** and its 16th–84th percentile, so for
HD 7977 the curve dips to 0.071 ly below a dot at 0.123 — the curve passing low
through its own error bar is the honest picture, not a mismatch. And the
selected curve is drawn in the interface accent blue rather than the star's
color, as is its 3-D motion trail: blue means *selected* here, and an orange
trail through a field of orange stars would say less than it costs.

Every other curve is stroked with the star's own temperature color, so the
spectral-type work above reaches the chart automatically — Sirius's line is
`rgb(203,219,255)` where it used to be the same cream as everything else.

### Full screen

**⛶ Full screen** beside Methods & limits, or `F`. The whole of `#stage` goes
fullscreen rather than the 3-D canvas alone: the list, the chart, the caption
and the orientation panel are part of the instrument, and fullscreening the
canvas would drop every one of them. `#stage` is already
`position: fixed; inset: 0`, so filling a screen needs no layout change — only
a background, because `<body>` is not painted behind a fullscreen element.

Three details that are easy to get wrong:

- **The button does not own the state.** Esc, the browser's own chrome and the
  window manager all leave fullscreen without going through it, so the label,
  `aria-pressed` and the tooltip are driven by the `fullscreenchange` event and
  the click handler only ever asks. A reload *inside* fullscreen also arrives
  with the right label, because `syncFullscreen()` runs at boot.
- **Escape is shared.** It already stopped the tour, and the browser uses it to
  leave fullscreen. One press now leaves fullscreen and nothing else; a second
  stops the tour. `Cmd`/`Ctrl`-`F` is left to the browser's find.
- **Every measured dimension changes.** The dock height, the detail card's
  ceiling and the chart's label budget are all read from pixels, so the change
  handler re-measures all three. The canvas resizes itself from `clientWidth`
  on the next frame and only needs to be asked for one.

Where the platform cannot do it — iOS Safari fullscreens video and nothing else
— the button is removed rather than left to do nothing when pressed. On narrow
layouts it travels into the filter sheet alongside Methods & limits, for the
same reason and by the same code.

**It also exposed a resize bug that had been latent all along.** Three canvases
decided whether to resize their backing store by testing the width alone:

```js
if (overlay.width !== w * dpr) { overlay.width = w * dpr; overlay.height = h * dpr; }
```

Going full screen from an already-maximized window is the one common resize
that changes *only the height* — the menu bar, tab strip and window chrome
disappear and the width was already the width of the display. So the test said
nothing had changed, the 3-D canvas resized correctly beside them, and the
browser stretched a 1400-pixel-tall overlay buffer over an 1680-pixel-tall box.
Every label and reticle then drifted downward in proportion to its own y:
**85 px at mid-screen**, which is what put the marker below its star.

Reproduced by resizing 1100×700 → 1100×840 and confirming the overlay's backing
height stayed at 1400 while the GL canvas went to 1680; fixed by testing both
dimensions on the overlay, the orientation gizmo and the chart. Both sides now
floor to integers as the GL canvas already did — assigning a fractional size
truncates, so the old comparison could never match again and reallocated the
canvas every frame.

Verified after the fix by reading the framebuffer back and comparing each star's
brightest pixel against where `project()` puts its reticle: 0–1 device pixels
across four viewport shapes, and 0–1 during the tour on isolated subjects. The
one apparent outlier, Alpha Centauri A at 12 px, is the search window finding
Alpha Centauri B beside it.

### Finding your way around

Once the view has been tumbled it is easy to lose track of which way the Galaxy
lies, so orientation is always shown explicitly rather than left to be inferred:

- An **axis gizmo** tracks the camera against the three Galactic axes — GC
  (Galactic center), ROT (direction of Galactic rotation) and NGP (North
  Galactic Pole), with the anticenter and south pole drawn faintly behind.
- **Level** looks along the Galactic plane toward the center with north up;
  **Top-down** looks down from the North Galactic Pole; **Reset** returns to the
  default three-quarter view. All three ease into place rather than cutting, so
  you can see how the new orientation relates to the old one.
- **Zoom** has explicit −/+ buttons and a logarithmic slider spanning 0.6–300 pc,
  with a readout in both parsecs and light years of how wide the view is at the
  Sun's distance. The scroll wheel still works.
- The **Oort cloud** is drawn as three orthogonal great circles at 1.5 light
  years — a sphere rather than another flat ring, because it is a shell you pass
  *inside* rather than a distance you pass. It appears only once the view is 12
  pc across or tighter; at the default 124 ly it would be four pixels of noise.
  This is what makes the close encounters legible: at +1.296 Myr Gliese 710 sits
  visibly within the shell, and HD 7977 is deeper still. The 1.5 ly figure is
  shared with the band the chart shades, from one constant, so the two views
  cannot come to disagree about where the cloud ends. It has no real edge, only
  a consensus about where it thins out.
- The plane rings are **labeled with their radius**, placed wherever the screen
  is clear. Each label is anchored at its ring's widest point on screen: seen
  edge-on a ring collapses to a line and every other point is foreshortened,
  which would otherwise print "20 pc" closer to the Sun than "5 pc". A ring whose
  true extent is off-screen goes unlabeled rather than lying about scale.
- **Nothing marks the Galactic center.** There used to be a ray out to 45 pc
  with a label on it, which is a marker for something 8,000 pc away drawn
  inside a map whose widest ring is 40 — wrong by two orders of magnitude
  wherever it is put, and at close zoom the label landed next to the Oort
  cloud. Direction is the axis gizmo's job; it reads as a compass rather than
  as a position, and it still says GC.
- The rings are drawn **four times, each offset by a fraction of a pixel**.
  `gl.lineWidth()` is capped at 1 on essentially every desktop driver — this
  context reports a maximum of exactly 1 — so a one-pixel, 0.13-alpha blue ring
  was the only thing on offer, and it disappeared on any screen with a little
  glare. Overdrawing shifted copies builds a line about two pixels wide and,
  because blending is additive, brightens it at the same time. The offset is
  applied after projection and scaled by `w`, so it stays a constant width on
  screen instead of growing with distance.

### Finding a star

Stars are small and often overlap: at the default view α Centauri A and B are
**7 pixels apart** and read as a single dot, so selecting one gave almost no
signal — the only feedback was a subtle change in a crowded field while the
detail card sat in the far corner.

Selecting a star now marks it in place with a ring pulse and leaves a steady
reticle and name label on it. Picking one from the list, or pressing **⌖ Focus**
in its card, also flies the view in: the camera looks at the midpoint of the
Sun–star line and swings perpendicular to it, so both ends stay in frame, and
the view radius widens automatically if the star would otherwise be culled. For
α Centauri that takes the view from 124 ly across to 8, at which point A and B
are plainly two stars. Clicking a star directly in the 3-D view marks it but
does not move the camera, since you are already looking at it.

Selection toggles: clicking the same star again clears it, as does Escape, the
card's ×, or clicking empty space. Because selecting can move the camera,
clearing it puts the view back — including the view radius, if focusing had to
widen it to keep the star from being culled. The one exception is when you have
driven the camera yourself since focusing: your view is then the one worth
keeping, so it is left alone. Choosing Level, Top-down or Reset likewise
discards the stored view, since that is an explicit decision about where to
look.

### Scrubbing time

The time cursor is draggable directly — grab the grip at the top of the chart or
the cursor line anywhere down its length, and a readout follows the pointer.
Everywhere else in the plot keeps hover-a-curve and click-to-select, and a drag
that wanders into that region will not select a star on release.

Playback has an explicit direction: **◀** runs into the past, **▶** into the
future, and whichever is running turns into a **❚❚** you can press to pause. A
single button that reversed itself at the ends made "run this forwards"
impossible to express, so reaching either end of the window stops rather than
silently turning around.

Stopping there used to leave a dead control. At +5 Myr, pressing **▶** set
playback going, the next frame clamped it, and the button flicked back —
visibly doing nothing, with no hint why. It now does what every media control
does at the end of a track: **restarts**. Pressing **▶** at +5 Myr jumps to
−5 Myr and runs forward; **◀** at −5 Myr does the mirror image. The button says
so before you press it — its label and tooltip switch to *Replay from −5 Myr* —
because a control that behaves differently at the ends should admit that in
advance rather than surprise you. Those labels are also the buttons' accessible
names, which previously read out as the bare glyph "❚❚" and claimed to be *Play
backwards* while showing a pause icon.
Keyboard: `J` back, `K` pause, `L` forward, and Space resumes in whichever
direction was last used.

The transport slider works in the same symlog units as the chart axis and its
track is inset to exactly the plot area, so the thumb sits pixel-aligned above
the cursor line at every position. A slider linear in years would be useless
here: it would crawl for four megayears and then cross the entire interesting
region — every encounter from Barnard's Star to Ross 248 — inside one pixel.
The readout always shows true elapsed years.

**Playback, however, advances in years, not in axis units.** It used to sweep
the axis at a constant rate, which is right for the cursor and wrong for
everything else. The real-time rate that implies runs from 0.06 to 31 Myr per
axis unit — a factor of **500** across the window — so the stars crawled near
the present and then visibly rocketed, which is not something stars do. Stars
move at a fixed speed through space, so the only way they can look like it is
for time to advance at a fixed rate. Measured on Gliese 710, playback now moves
it 4.42 pc per second of wall clock at every point in the window; the residual
0.2% is real gravitational acceleration rather than an artifact of the axis.

The cost is that the cursor now crosses the compressed middle of the chart
quickly. That is the honest way round: the distortion belongs to the axis, which
exists to make the near-term encounters *reachable by the scrubber*, not to the
physics.

The speeds are labeled with the rate — **10 kyr/s, 50 kyr/s, 250 kyr/s, 1
Myr/s** — because that is the only thing that makes them predictable. At the
default 50 kyr/s a typical 50 km/s star takes about sixteen seconds to cross a
40 pc field, which is a speed you can actually follow; the two fast settings are
for getting somewhere rather than for watching. The full 10 Myr window takes
200 s, 33 s or 10 s accordingly. The tour's own time moves are eased but
likewise linear in years.

## Known limitations
<a id="known-limitations"></a>

**The intervals are measurement-sensitivity ranges, not confidence intervals.**
They answer "how much does this prediction move when the published measurements
are resampled within their stated errors", and no more. Not propagated:

- The Gaia **parallax zero-point correction** is not applied. Lindegren et al.
  (2021) give a magnitude-, color- and position-dependent offset of order
  0.02–0.05 mas; implementing it needs five more Gaia columns and their
  published coefficient table, and a wrong correction is worse than a stated
  omission. What is done instead is to measure it: stage 9 re-integrates the
  catalog at Z = −0.017, −0.030, −0.050 and +0.010 mas and reports the largest
  shift per star as `zp_dmin_pc`, with a flag where it exceeds the measurement
  interval.
- Radial-velocity errors are **assumed** where a catalog publishes none
  (2 km/s), and implausibly small published errors are floored at 0.1 km/s.
  The CSV keeps the untouched published value in `rv_error`, so which rows were
  substituted is reconstructible, but there is no per-row provenance column
  saying where each error came from.
- Uncertainty in the **Galactic potential** and in the **solar parameters** is
  not folded into these intervals — it is measured separately and reported per
  star as `sys_dmin_pc`. See [Does the Galaxy model matter more than
  Gaia?](#galaxy-model) — for a tenth of the candidates within 5 pc, it does.
- 256 draws is adequate for a median and a 68% range, thin for tails.

**There is no asset fingerprinting.** `meta.json`, `stars.json` and `cheb.bin`
are fetched with `cache: 'no-cache'`, so the data can never go stale behind a
current `app.js` — but that policy lives *in* `app.js`, which is itself served
under whatever rules the host applies. A hosted deploy should fingerprint or
explicitly no-cache `index.html` and `app.js`; `serve.py` does the latter, and
nothing in this repository does the former.

**Accessibility is partial.** Controls are labeled and keyboard-operable, rows
are focusable, focus is visible, contrast meets 4.5:1, reduced motion is
respected, and the layout carries breakpoints for phone and short-laptop
windows. But the canvases only carry `role="img"` pointing at the text lists,
which is a floor, not a substitute for an accessible chart, and the phone layout
is a set of breakpoints rather than a deliberately designed small-screen view.

- Stars are integrated as test particles. Encounters between stars are ignored,
  as is the Galactic bar and spiral structure.
- Astrometric errors are treated as independent in the Monte Carlo; Gaia
  publishes a full correlation matrix that is not used here.
- Radial velocities are epoch measurements. For unresolved binaries they include
  orbital motion, which is flagged but not corrected.
- **2,324 stars are excluded outright for want of a radial velocity**, in
  either Gaia RVS or the literature. Position and proper motion are not enough
  to propagate a star; without the line-of-sight component the trajectory is
  undetermined. They are absent from the catalog rather than shown at rest,
  which would be worse, but it does mean "12,218 stars" is not a complete
  census of future neighbors — it is a census of those whose motion is
  actually known.
- Six curated objects are missing for want of any radial velocity: WISE
  0855−0714 (no Gaia astrometry — detected only in the infrared), Luhman 16,
  and Sirius B among them.
- Predictions for stars currently beyond ~100 pc rest on extrapolating small
  parallaxes over millions of years; consult the Monte Carlo interval before
  taking any of them seriously.
- **The star count is only complete near the present.** Selection required a
  star to be within 30 pc today or to approach within 10 pc at some point, so
  the sample thins out as you scrub away from t = 0 — the visualization reports
  779 stars within 20 pc at +1.3 Myr against 2,233 today (counting the whole
  shipped set; the default view enables "Hide flagged measurements" and so
  reports fewer). Some of that is real
  dispersal, but much of it is stars that will arrive in the neighborhood and
  were never selected because they are far away and unremarkable now. Treat
  counts away from the present as a lower bound.
- Closest approaches are clipped to the ±5 Myr window, so a star whose true
  perihelion lies beyond it is reported at the window edge rather than
  extrapolated.

## Serving it

`scripts/serve.py` binds **127.0.0.1** only (pass `--public` if you really want
the LAN to reach it), disables directory listings, and sends the headers a
hosted deployment should also send:

```
Content-Security-Policy: default-src 'none'; script-src 'self';
    style-src 'self' 'unsafe-inline'; img-src 'self' data:; connect-src 'self';
    base-uri 'none'; form-action 'none'; frame-ancestors 'none'
X-Content-Type-Options: nosniff
Referrer-Policy: no-referrer
Cache-Control: no-cache, must-revalidate
```

The page makes **no outbound requests** — no CDN, no fonts, no analytics, no
telemetry — so `default-src 'none'` with `connect-src 'self'` costs nothing and
means an injected script could not exfiltrate anything either. `'unsafe-inline'`
is needed for styles only because the stylesheet lives in a `<style>` block in
`index.html`; scripts get no such exemption. Catalog names are external data
and reach the DOM as text nodes, never as markup.

### Deploying it somewhere you do not control the headers

A hosted copy is a file tree on someone else's server, so the policy above is
**also embedded in `index.html` as a `<meta http-equiv>`**. Without that, the
headers exist only under `serve.py` and the public site has no policy at all —
which is exactly what happened on the first deployment and is why the meta tag
is there. Two directives cannot be set that way and still need server config if
you want them: **`frame-ancestors`** (clickjacking) and
**`X-Content-Type-Options: nosniff`**.

**Keep injected third-party scripts off this path.** Cloudflare Web Analytics,
and anything else that rewrites HTML at the edge, adds a `<script>` the CSP will
refuse — the policy holds, no request goes out, but the browser console fills
with violations and the analytics silently records nothing. The decision here is
to leave the CSP strict and keep the injection turned off rather than
allow-listing a third party, so that "no outbound requests" stays literally
true. If you ever want the analytics, allow-listing it means the sentence above
is no longer accurate and should be amended in the same commit.

After any deploy, the one-line check that the policy actually shipped:

```bash
curl -s https://equn.me/star-neighbors/ | grep -c 'http-equiv="Content-Security-Policy"'
```

## Requirements

The 3-D view needs WebGL2 (float textures and `gl_VertexID`). Without it the
page degrades rather than failing: the 3-D panel and orientation controls hide
themselves and the distance–time chart, search, encounter tables and Monte Carlo
figures — all Canvas 2-D and DOM — keep working, with a notice explaining why.

`web/data/` must be served over HTTP, not opened as a file, and must contain all
three of `meta.json`, `stars.json` and **`cheb.bin` (1.9 MB)**. The binary holds
every trajectory; without it there is nothing to animate. The loader checks that
its length matches the star count and reports precisely what is wrong rather
than failing obscurely.

## Review history

This has been through several rounds of external review and one self-directed
audit against the refereed literature, and most of them changed the science, not
the presentation. Recording that here because the failures are more instructive
than the result:

- **The uncertainty model contradicted the dynamical model.** Draws were solved
  with straight-line motion while the interface displayed the integrated epoch
  beside them and ranked on the straight-line distance. The justification was a
  real measurement — under 0.005 pc of difference — taken from the near-term
  encounters it was verified on and applied to multi-Myr candidates where the
  same stage had measured 0.53 pc. Now every draw is integrated.
- **The candidate screen had no margin.** Screening and scientific interest both
  used 10 pc, so deflection could carry a genuine close pass out of the
  catalog before it was ever computed. The screen is now 20 pc.
- **Boundary cases were reported as events.** Trajectories still closing at the
  edge of the window were quoted as closest approaches. They are now labeled as
  bounds.
- **The future/past split used a different model from the epoch it printed.**
  Rows were assigned to "closest ahead" or "closest behind" by the nominal
  epoch while displaying the Monte Carlo one, so eight stars appeared in the tab
  opposite to the date beside them. Both now come from the same column.
- **Censoring was described in one direction only.** A trajectory pinned to the
  −5 Myr edge was labeled "still approaching", which is true at +5 Myr and
  false at −5 Myr: there the minimum lies *before* the window. The wording is
  now directional, and neutral for the seven stars whose draws straddle both
  boundaries.
- **A duplicated bright star survived the first fix.** ρ Ori arrived once from
  Gaia and once from Hipparcos, 75 pc and 107 pc, because their parallaxes
  disagree by 30.1% and the positional deduplicator requires 30%. Loosening that
  threshold would merge genuinely distinct neighbors, so ambiguous pairs now go
  to SIMBAD's identifier tables and are merged on evidence. A validation
  assertion keeps them merged.
- **The fix was right and the documentation was not.** After the censoring
  wording was made directional in the interface, the CSV column reference and
  the shipped `flag_bits` schema still described the one-directional version,
  and the metadata claimed the sign of `mc_tmin_med` always identifies the edge
  — which is exactly what the straddling cases disprove. Both now describe all
  three outcomes.
- **A gate that did not gate.** The CSV exporter checked that columns existed
  and that `mc_model` was integrated, but would still publish blank or
  non-finite numbers underneath. The rule now lives in one module both
  exporters call, and rejects empty cells, `nan` and `inf` as well;
  `test_gates.py` covers twelve injected defects and seven valid inputs.
- **The Methods panel omitted a source the README had.** `meta.json` credited
  only Gaia and Hipparcos, so the one star whose astrometry comes from neither
  was invisible in the application even after the README named it. That string
  is now counted from the catalog rather than written by hand.
- **Stale data could survive a redeploy.** The three data files were fetched
  with default caching, so a browser holding copies from an earlier build would
  render old science with nothing to show for it. They are now fetched with
  `cache: 'no-cache'`, and `serve.py` now sends `Cache-Control: no-cache`
  rather than `no-store` — storage is allowed, revalidation is required, so an
  unchanged 1.9 MB blob costs a 304 instead of a full refetch. Found because a
  stray `python3 -m http.server` was holding the dev port and serving the bundle
  *with* caching, which is precisely the failure `serve.py` exists to prevent.
- Smaller: barycentric velocities inheriting the top quality grade, implausible
  1e-4 km/s errors, a quality filter conflated with prediction confidence,
  metadata bits set in the data but absent from the schema, exporters that
  turned missing science columns into zeros, a diagnostic script that compared
  the integrated model against itself while calling one side "straight-line",
  a detail card that overlapped the orientation panel on a 720px-tall window,
  and a loader that validated the trajectory blob against its fallback
  polynomial degree before reading the build's actual one.

A later round found a pattern rather than a set of unrelated bugs. Five defects
in a row — the star colors that never reached a pixel, the chart curve that
stepped over its own perihelion, the tour subject sitting behind the caption,
the spectral-type fallback that caught every famous star, and the Max spread
filter that moved nothing — were all the same mistake: **the input was verified
and the output was assumed to follow.** The colors were correct in
`stars.json`. The perihelion was correct in the catalog. The tour ran without
throwing. Each check passed, and each checked the wrong end.

What the last audit changed, accordingly, is that it recomputes ground truth
independently and compares it against what is actually *displayed* — pixels
read back from the framebuffer, text read out of the DOM. That found:

- **Max spread not applied to the 3-D view** (above), the one remaining place
  where two code paths stated the same rule.
- **A factual error in the narration**, introduced while simplifying the
  language a day earlier: HD 7977's uncertainty described as "roughly three
  times nearer or further than our best guess", where the edges are 1.95× and
  1.68× and it is the *span* that is threefold. Now "the far end of that range
  is three times the near end."
- **`NaN` could permanently poison the epoch and the sky field of view.**
  `Math.min/max` propagate it rather than clamping, and nothing downstream
  recovers: every position becomes `NaN`, the view empties and stays empty. No
  control can produce one — every slider coerces — so this is insurance, not a
  repair.

Clean under the same method: the detail card against recomputed values for 159
stars; all 480 list rows across the three tabs, including ordering, the sign of
the epoch and whether each `aria-label` carries the number beside it; 22 of the
23 quantitative claims in the narration; clamping of time, radius and field of
view at their extremes, and out-of-range selection indices; the listbox's
roving tabindex, arrow/Home/End/Enter behavior and canvas text alternatives;
and the no-WebGL fallback, where the only function touching `gl` without a
guard turns out to be reachable only from inside one that has already returned.

Two of those "findings" were the harness being wrong, which is worth recording
too: a grade-C star read as still-drawn because a 9×9 pixel window caught a
neighbor 4 px away, and the Max spread slider read as unlabeled because the
probe looked for `aria-label` and the control uses `<label for>`.

A final self-audit before publication found the rest, none of them raised by a
reviewer at the time:

- **The largest remaining uncertainty was named but never measured.** Every
  earlier draft said the Galactic potential was "not propagated" and left it
  there. Measuring it (stage 9) showed it is negligible nearby and *larger than
  the whole quoted interval* for 165 of the 1,697 candidates within 5 pc,
  including the second-ranked future encounter. The parallax zero-point, the
  other standing omission, was measured the same way — it dominates for a
  further 15. Both are now per-star columns, flags, and a section of this file.
- **Keyboard shortcuts stole keys from focused controls.** Space on a list row
  selected the star *and* started playback; Space on a button pressed it *and*
  started playback; an arrow key on the speed menu changed the speed *and*
  scrubbed time. The global handler only excluded `<input>`.
- **The star list had 160 tab stops.** Every row carried `tabindex="0"`, so
  tabbing past the list took 160 presses. Now a roving tab stop with arrow,
  Home and End keys, which is the listbox pattern.
- **`prefers-reduced-motion` was honored only in CSS.** The camera fly-to and
  the selection pulse are canvas animations and ignored it entirely.
- **An empty list said nothing.** A search with no matches produced a blank
  panel rather than explaining which control had emptied it.
- **The dev server bound every interface and listed directories.** It now binds
  127.0.0.1 unless asked otherwise, refuses listings, and sends a CSP.
- The detail card's two remaining `innerHTML` writes became DOM nodes, so no
  catalog-derived value reaches the parser as markup anywhere.
- The comment claiming this potential gives *v*_c(8.122 kpc) ≈ 229 km/s was
  wrong — it gives 231.5. The consequence was measured (≤ 0.0025 pc) rather than
  assumed.

A fourth review then prompted the last two changes:

- **Alpha Centauri's parallax was wrong by 5.6% between its own components.**
  Found by checking present-day distances against the published list of systems
  within 20 ly — the one comparison the project had never run. This file had
  described the defect and left it, for want of a better number; the reference
  supplied one.
- **Gaia's astrometric correlations are now propagated.** They were the one
  thing Gaia publishes that this pipeline ignored. Median |r| is 0.15 and the
  largest is 0.91; adopting them changes the median interval width by 0.5% but
  widens 46 of 400 close candidates by more than 10%. `test_draws.py` asserts
  the draw code's properties rather than trusting it.
- **A retry that did not retry.** Half the re-fetch died on an
  `http.client.IncompleteRead` — ESA closing a chunked response early — which is
  an `HTTPException` and so slipped past a handler that listed only `URLError`,
  `ConnectionError` and `TimeoutError`. It threw away a completed server-side
  job. Worse, the half-updated raw set would have left the *nearest* stars
  silently drawing without correlations, so validation now asserts that Gaia
  rows carry them, checked separately within 30 pc.
- **A brittle assertion, honestly re-cut.** HD 7977's published 0.0641 pc had
  always sat at the very top of its interval; the correlated draws moved the
  84th percentile from 0.0645 to 0.0634 and the containment test failed. The
  test is now a stated 10% tolerance, and the README says the published value is
  at the edge rather than inside — the fix was to the claim, not just the test.

- **The comparison set was a compilation, not the literature.** Validation ran
  against 48 hand-transcribed rows. It now also runs against the refereed
  catalogs on Gaia source_id: **1740 paired encounters at a median ratio of
  1.000**, and 53 same-data-release encounters whose epochs agree
  100 for 100. Two long-standing claims died on contact with the source.
  HD 7977, called "the one genuine unexplained disagreement", was a *linear*
  reference value being compared with an integrated one — the compilation's
  0.478 ly is 0.14656 pc and this project's straight-line answer is 0.14659 pc,
  and HD 7977 turns out to be the most deflection-sensitive star in the list at
  85%. Against the refereed integrated value it agrees. Gliese 710 went the
  other way: the 0.0636 pc figure earlier declared unsourceable is
  Bailer-Jones (2022), and it cannot be reproduced from the astrometry
  published beside it.

What should be read as open: **13 of the 53 refereed close encounters differ
from this catalog by more than any potential could account for**, Gliese 710
among them, and its own inputs match ours to 0.1%. That is a real unresolved
difference with an independent group, not a rounding disagreement, and it is
[written up in full](#bailer-jones) rather than buried here. The intervals also
still omit the parallax zero-point and the potential's own uncertainty — both
are measured and reported separately, as `zp_dmin_pc` and `sys_dmin_pc`, rather
than folded in.

## Sources

Gaia Collaboration, *Gaia DR3* (2023) · van Leeuwen, *Hipparcos, the New
Reduction* (2007) · Wenger et al., *SIMBAD* (2000) · Price-Whelan, *gala* (2017)
· GRAVITY Collaboration (2018) · Bennett & Bovy (2019) · Eilers et al. (2019) ·
Schönrich, Binney & Dehnen (2010) · Akeson et al. (2021), α Cen AB parallax

Validated against: Bailer-Jones et al., *New stellar encounters discovered in
the second Gaia data release*, A&A 616, A37 (2018) — VizieR
`J/A+A/616/A37/table23` · Bailer-Jones, ApJL 935, L9 (2022) — VizieR
`J/ApJ/935/L9/table12` · Wikipedia, *List of nearest stars*, for the compiled
sub-5-ly list and the 20-ly completeness census. Reference tables are fetched
from VizieR and cached under `data/raw/`; they are not redistributed here.

## Acknowledgments

This work has made use of data from the European Space Agency (ESA) mission
[Gaia](https://www.cosmos.esa.int/gaia), processed by the Gaia Data Processing
and Analysis Consortium
([DPAC](https://www.cosmos.esa.int/web/gaia/dpac/consortium)). Funding for the
DPAC has been provided by national institutions, in particular the institutions
participating in the Gaia Multilateral Agreement.

This research has made use of the SIMBAD database and the VizieR catalog
access tool, CDS, Strasbourg, France (DOI 10.26093/cds/vizier). The original
description of the VizieR service was published in A&AS 143, 23.

## License

The **code** in this repository is MIT (see `LICENSE`).

The **derived catalog** — `data/star_encounters.csv` and everything under
`web/data/` — is released under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). It is a derived
product: Gaia DR3 is published by ESA under CC BY-SA 3.0 IGO, and the
underlying measurements remain ESA's. Cite this repository and the Gaia mission
together, not this repository alone.

Reference tables used only for validation — the Bailer-Jones catalogs from
VizieR and the Wikipedia compilation — are fetched at run time and cached under
`data/raw/`, which is git-ignored. They are **not** redistributed here, so no
part of this repository is subject to their terms.
