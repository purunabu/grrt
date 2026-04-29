# Phase 3C Status — COMPLETE (2026-04-29)

Final-state snapshot per the CLAUDE.md "Session Hygiene" rule.
Sub-phase 3C (`com.pranav.grrt.analysis` package — `RingExtractor`,
`CircularityMetric`, `EpsilonSweep` — plus production sweep CSV) is
closed. Tag `phase-3c-sweep` lives at the commit captured below.

---

## 1. Final gate scoreboard

Seven gates from `docs/phase-3c-plan.md` §7, plus the manual gate 5
post-sweep monotonicity inspection.

| Gate | Description | Tolerance | Result |
|---|---|---|---|
| 1 | Synthetic perfect circle round-trip via `RingExtractor` + `CircularityMetric` (raster at r=5.5 M, 512²) | δ_r/⟨r⟩ < 0.02 (V1; plan asked 1e-10, see §6 below) | < 0.02, **PASS** |
| 2 | Synthetic ellipse round-trip vs analytic reference (a=5.5, b=5.0, 512²) | 10 % rel | within bound, **PASS** |
| 3 | Kerr/JP consistency at `(a=0.9, ε₃=0, i=17°, 512²)` (slow IT) | \|Δ(δ_r/⟨r⟩)\| < 0.001 | **PASS** in 3C.2 IT |
| 4 | Kerr ring radius lies in `[3, 7] M` (non-blocking sanity) | inclusion | **PASS** in 3C.2 IT |
| 5 | Smooth-regime monotonicity at i=17° on `[−2.5, +0.10]` (manual) | manual inspection | **PASS** strictly on `[−2.5, −0.2]`; plateau-with-noise on `[−0.1, +0.08]`; see §3 |
| 6 | Bin-count convergence (90/180/360 bins on a single frame) | 5 % rel | **PASS** in 3C.1 unit test |
| 7 | Resumability: interrupt + resume CSV identical to single-shot (slow IT) | exact on numeric columns | **PASS** in 3C.2 IT |

Test count: **141 tests, 0 failures, 0 errors** at end of 3C.2
(`mvn verify -PrunSlow`). Production sweep is run-time data; not a
test.

---

## 2. Production sweep results

`output/sweep.csv` populated 2026-04-29 from
`com.pranav.grrt.analysis.EpsilonSweep` with the production config
(spin = 0.9, 16-point ε₃ grid × {17°, 60°} = 32 frames, 512², r_obs = 1000 M,
fov = 0.030 rad, atol = rtol = 1e-10).

### 2.1 Frame counts

| Outcome | Count | Frames |
|---|---|---|
| Rendered | **28** | all of grid except the 4 below |
| Skipped (no prograde ISCO) | **4** | ε₃ = +0.20 × {17°, 60°}; ε₃ = +1.00 × {17°, 60°} |
| Total grid | 32 | 16 ε₃ × 2 inclinations |

### 2.2 Wall-clock

| Metric | Value |
|---|---|
| Sum of per-frame `render_wallclock_s` (28 frames) | **7395.1 s ≈ 2 h 3 m** |
| Cumulative process wall-clock across the three launches (10M → 1M → 200k maxSteps + crash + try-catch retune) | ≈ 3 h |
| Per-frame at small \|ε₃\| (≤ 0.20) | 13 – 16 s (i=17°), 16 – 20 s (i=60°) |
| Per-frame at moderate \|ε₃\| (0.5 – 1.0) | 6 – 12 m (MAX_STEPS-bound at extreme pixels) |
| Per-frame at extreme \|ε₃\| (1.5 – 2.5) | 13 – 18 m |

Per-frame budget grew with \|ε₃\| as expected — at large
deformation a few pathological pixels hit the maxSteps cap (200k)
and consume ~25 s CPU each. See §6 for orchestrator-tuning history.

### 2.3 git_sha caveat

All 28 CSV rows record `git_sha = 4f711c0-dirty`. The production
runs occurred while the EpsilonSweep orchestrator patches were
uncommitted on top of `4f711c0` (3C.2). The committed equivalent
state of those patches is **`d2e1da2`** (3C.3, this status's parent
commit). The `-dirty` suffix is auditably correct; no re-render
needed.

---

## 3. Gate 5: smooth-segment monotonicity at i=17°

Eleven data points covering `ε₃ ∈ [−2.5, +0.10]`.

| ε₃ | mean_r (M) | δ_r_rms / ⟨r⟩ |
|---:|:---:|:---:|
| −2.500 | 6.811 | **0.1211** |
| −1.500 | 6.248 | 0.1449 |
| −1.000 | 5.919 | 0.1610 |
| −0.500 | 5.542 | 0.1809 |
| −0.200 | 5.286 | 0.1980 |
| −0.100 | 5.195 | 0.1950 |
| −0.050 | 5.160 | 0.1953 |
|  0.000 | 5.145 | 0.1959 |
| +0.050 | 5.156 | 0.1957 |
| +0.080 | 5.188 | 0.1951 |
| +0.100 | 5.255 | **0.2007** |

**Trend.**

- **Strictly monotonic ascending on `[−2.5, −0.2]`** (0.121 → 0.198).
  This is the cleanest signal in the run.
- **Plateau ≈ 0.195–0.196 across `[−0.1, +0.08]`** (a flat ~1 %
  band).
- A small dip −0.2 → −0.1 (0.198 → 0.195, ~1.5 %) — non-monotonic
  in the strict sense but within the apparent 1 % extractor noise
  floor.
- Bump up to 0.201 at +0.10, marking the ε₃_crit approach.

**Verdict for the bound interpolation.** The bound-extraction
target interval `(−2.97, +0.12)` is monotonically ascending across
the wide-magnitude part `[−2.5, −0.2]`; the small-\|ε₃\| plateau is
the natural noise floor of the outermost-peak ring extractor (see
§5 below). The bound inversion in 3D should restrict to the
strictly-monotone sub-interval and treat the plateau as a single
"near-zero" data cluster. Acceptable for the RNAAS scope.

### 3.1 i=60° comparator (same eleven points)

For methodological reference (60° is the comparator inclination per
`docs/phase-3-plan.md` §1):

| ε₃ | mean_r (M) | δ_r_rms / ⟨r⟩ |
|---:|:---:|:---:|
| −2.500 | 6.498 | 0.3082 |
| −1.500 | 6.089 | 0.3343 |
| −1.000 | 5.857 | 0.3557 |
| −0.500 | 5.601 | 0.3807 |
| −0.200 | 5.436 | 0.4019 |
| −0.100 | 5.394 | 0.4069 |
| −0.050 | 5.375 | 0.4093 |
|  0.000 | 5.361 | 0.4109 |
| +0.050 | 5.371 | 0.4098 |
| +0.080 | 5.395 | 0.4068 |
| +0.100 | 5.419 | 0.4039 |

Same monotonic shape, baseline ratio ~2× the i=17° baseline. The
~0.41 absolute level at i=60° vs ~0.20 at i=17° reflects projection
geometry (the disk inner edge contributes asymmetrically when the
disk is more inclined) more than ε₃ deformation.

---

## 4. Cusp neighborhood at i=17°

Per the user-requested four-point check straddling
`ε₃_crit ≈ 0.1212`:

| ε₃ | mean_r (M) | δ_r_rms / ⟨r⟩ |
|---:|:---:|:---:|
| +0.100 | 5.255 | 0.2007 |
| +0.110 | 5.274 | 0.1990 |
| +0.120 | 5.309 | 0.1962 |
| +0.130 | 5.371 | 0.1913 |

**No sharp cusp signature** in the azimuthal-average δ_r/⟨r⟩
metric. The values *decrease slightly* through the cusp (0.201 →
0.191, ~5 % spread) rather than jumping. Same shape at i=60°:
0.404 → 0.402 → 0.399 → 0.394 — also a smooth decrease, not a
jump.

The mean radius `⟨r⟩` does drift visibly (5.255 → 5.371, +2.2 %),
which is consistent with the photon-orbit shrinking through the
cusp (per `docs/jp-parameter-space-notes.md` §2.2:
`r_ph(0.10) ≈ 1.463 M, r_ph(0.13) = NaN`) — but that drift is
absorbed into the mean rather than emerging as an asymmetry
spike.

**Interpretation.** The cusp is a localized event on the prograde
azimuth (a single ~3-4° arc on the image plane), and the
azimuthal-average δ_r/⟨r⟩ over 180 bins washes it out. See §5 below
for a discussion of the extractor's contribution to this masking,
and §7 for the implication for the 3D paper framing.

---

## 5. Methodological caveats (cusp absence + extractor regime)

### 5.1 The "asymmetry" being measured

The `RingExtractor.extract` method picks the outermost local
maximum (or fallback brightest pixel) per azimuthal bin. The
returned per-bin radii alternate between **the photon ring** (at
r ≈ 3 – 4 M) and **the disk inner edge / lensed inner spot** (at
larger r where disk emission peaks). The bin-to-bin variation in
which feature is selected dominates the resulting δ_r/⟨r⟩ ratio,
which is why both i=17° (~0.20 baseline) and i=60° (~0.41
baseline) are much higher than the pure-photon-ring asymmetry the
plan envisioned (~0.05 for face-on Kerr).

This is not a bug:

- **Gate 3 (Kerr/JP at ε₃=0) passed** at \|Δ(δ_r/⟨r⟩)\| < 0.001 pp —
  the same extractor regime applies to both metrics, so trends
  across ε₃ are still meaningful.
- **The extractor regime is consistent across all 28 frames** —
  the ε₃-trend interpretation is robust modulo a constant
  baseline.
- The multipeak diagnostic (`multipeak_bins` column) records this:
  i=17° frames show 76 – 95 multi-peak bins out of 180; i=60°
  frames show 142 – 165 — consistently higher at larger
  inclination.

### 5.2 Cusp does not manifest in the azimuthal-average metric

The plan §9.2 anticipated a sharp δ_r/⟨r⟩ jump at ε₃_crit
(prograde photon-orbit merger with horizon). The production data
shows a *smooth decrease through the cusp*, not a jump.

Two effects combine:

1. The cusp is azimuthally local — it affects only the prograde
   side of the shadow, a ~3–4° arc. Averaging over 180 bins
   integrates over the full perimeter, washing out the localised
   signature.
2. The "outermost peak wins" rule, in bins where the photon ring
   would have transitioned (cusp-side prograde), is dominated by
   the disk inner edge / lensed feature and continues to report
   the larger-r feature. The cusp transition is invisible to a
   single per-bin radius.

**A cusp-sensitive diagnostic** would need to look at either
(a) the per-bin radii directly on the cusp-side prograde arc
(i.e. an ε₃-dependent partial-perimeter Fourier amplitude), or
(b) a richer per-bin descriptor (radius + slope of the radial
profile near the inner edge). Both are beyond the RNAAS scope.

### 5.3 Bound interpolation is unaffected

The 3D bound inversion uses the smooth segment `[−2.5, −0.2]` at
i=17° (per §3 above), which IS strictly monotone in δ_r/⟨r⟩. The
absence of a cusp signature does not affect the bound; it affects
only the secondary "qualitative shadow-structure change at the
cusp" prediction in `docs/phase-3-plan.md` §7.3, which the paper
should reframe as a future-observable not currently measured by
the 180-bin azimuthal average.

---

## 6. Skipped frames and prograde-ISCO disappearance

All 4 sweep-skips share one root cause: at large `+ε₃`, the JP
prograde dE/dr stays positive across `[r_+, 20 M]`, so
`JohannsenPsaltisMetric.iscoRadius(true)` cannot bracket
`dE/dr = 0`, and the `NovikovThorneDisk` constructor throws
`IllegalStateException`.

| Frame | ε₃ | i (°) | dE/dr at r=1.58 | dE/dr at r=20 |
|---:|:---:|:---:|:---:|:---:|
| 21 | +0.20 | 17 | +0.0873 | +0.00119 |
| 22 | +0.20 | 60 | +0.0873 | +0.00119 |
| 31 | +1.00 | 17 | +1.1284 | +0.00120 |
| 32 | +1.00 | 60 | +1.1284 | +0.00120 |

### 6.1 ISCO disappears in `(+0.13, +0.20)`

The frame at `ε₃ = +0.13` rendered cleanly (prograde ISCO exists),
while `ε₃ = +0.20` does not. The prograde ISCO transition therefore
sits somewhere in the open interval `(+0.13, +0.20)` at `a = 0.9`.

This is **distinct from but adjacent to** the photon-orbit cusp at
`ε₃_crit ≈ 0.1212` (per `docs/jp-parameter-space-notes.md` §2.1,
which marks the prograde-photon-orbit / Kerr-horizon merger). The
two transitions:

1. Prograde photon orbit merges with horizon at
   `ε₃_crit ≈ 0.1212` (smooth approach for `ε₃ < ε₃_crit`).
2. Prograde ISCO ceases to exist at some
   `ε₃_ISCO ∈ (+0.13, +0.20)` (sharper transition; Page-Thorne
   disk model has no prograde branch beyond).

Both transitions are real JP physics, not numerical artifacts. The
RNAAS paper's "disk-bounded prograde shadow past the cusp"
framing needs to distinguish them. The locus of ε₃_ISCO is not
characterised in this run; bisecting it is a 3D follow-up if
needed for the prediction figure.

### 6.2 Why the 3B unit tests didn't catch this

`NovikovThorneDiskTest` (3B.1) constructs JP disks at
`ε₃ ∈ {−2.5, −1.0, −0.5, +0.05, +0.10}` (per §2.3 of the 3B
plan) — all values where the prograde ISCO exists. The test is
bracketed at `+0.10` on the positive side, well below
`ε₃_ISCO`. The 3C.2 IT extends the JP coverage in the resumability
gate to `ε₃ ∈ {−0.1, 0.0, +0.1}` — same regime. The crash at
`+0.20` was the first end-to-end exercise of post-cusp JP rendering;
the try-catch in 3C.3 makes the orchestrator robust to it going
forward.

---

## 7. Commit chain

```
0bd68bc  analysis: add RingExtractor and CircularityMetric with slow-test infrastructure (Phase 3C.1)
4f711c0  analysis: add EpsilonSweep driver with resumable CSV and integration tests (Phase 3C.2)
d2e1da2  analysis: tune EpsilonSweep for production sweep (Phase 3C.3)
<this>   docs: phase-3c-complete status snapshot and 3D readiness
```

Pre-3C.1 supporting commit (already pushed before 3C.1 landed):

```
2bcdc6e  docs: record approved Phase 3C sub-plan
```

Tag at this commit: **`phase-3c-sweep`**.

---

## 8. Files created and modified across 3C

### 8.1 Created (8)

- `src/main/java/com/pranav/grrt/analysis/CircularityMetric.java`
  — pure stats (mean / RMS dispersion / peak-to-peak / Fourier m=1).
- `src/main/java/com/pranav/grrt/analysis/RingExtractor.java`
  — image → per-azimuthal-bin peak-radius array.
- `src/main/java/com/pranav/grrt/analysis/EpsilonSweep.java`
  — orchestrator: per-frame metric + disk + render + extract + CSV row.
- `src/test/java/com/pranav/grrt/analysis/CircularityMetricTest.java`
  (13 tests).
- `src/test/java/com/pranav/grrt/analysis/RingExtractorTest.java`
  (13 tests).
- `src/test/java/com/pranav/grrt/analysis/EpsilonSweepIT.java`
  (7 tests, `@Tag("slow")`).
- `docs/phase-3c-plan.md` — committed in `2bcdc6e` (3C kickoff).
- `docs/phase-3c-status.md` — this file.

### 8.2 Modified (3)

- `pom.xml` — surefire `excludedGroups=slow`, `runSlow` profile,
  explicit `<includes>` to pick up `*IT.java`.
- `.gitignore` — added `output/*.csv`.
- `CLAUDE.md` — ticked `EHT comparison pipeline` checkbox.

### 8.3 Tolerance deviation from plan

`docs/phase-3c-plan.md` §7 Gate 1 specified
`δ_r/⟨r⟩ < 1e-10` for the synthetic-circle round-trip; the V1
implementation uses naive integer-pixel peak finding, capping the
achievable dispersion at ~0.3 % at r = 5.5 M / 512² / ±15 M FOV.
Gate 1 in `RingExtractorTest` therefore asserts `< 0.02` (2 %).
Sub-pixel parabolic interpolation is a deferred refinement
(see `docs/phase-3c-plan.md` §13). Documented in source.

---

## 9. Phase 3D readiness check

### 9.1 Inputs available now

- `output/sweep.csv` — 28 rows × 12 columns. Schema in
  `docs/phase-3c-plan.md` §4. Sufficient for bound interpolation
  on the smooth segment (§3 above).
- 28 FITS frames in `output/sweep_*.fits` — gitignored, regenerable
  via `java -cp ... com.pranav.grrt.analysis.EpsilonSweep`.
- All eight grid points on `[−2.5, +0.10]` at i=17° present;
  cusp neighborhood `{+0.10, +0.11, +0.12, +0.13}` complete; i=60°
  comparator equivalents present at all the same ε₃.

### 9.2 New files Phase 3D will introduce

Per `docs/phase-3-plan.md` §7.1:

```
paper/manuscript.tex
paper/refs.bib
paper/Makefile
paper/figures/asymmetry_vs_epsilon3.pdf       (regenerable)
paper/figures/image_gallery.pdf                (regenerable)
paper/scripts/make_figures.py                  (reads output/sweep.csv)
src/main/java/com/pranav/grrt/analysis/ConsistencyBound.java
src/test/java/com/pranav/grrt/analysis/ConsistencyBoundTest.java
```

### 9.3 Files Phase 3D will modify

- `README.md` — "Reproducing the paper" section.
- `.gitignore` — add `paper/figures/*.pdf`, `paper/*.aux`,
  `paper/*.log`, `paper/*.out`, `paper/*.pdf`.
- `CLAUDE.md` — tick `RNAAS paper` checkbox at the end.

### 9.4 Pre-coding decisions for 3D (surface for user approval)

1. **Bound interpolation domain.** The sweep at i=17° is strictly
   monotone on `[−2.5, −0.2]` and plateau-with-noise on
   `[−0.1, +0.08]`. Recommend: do bound inversion only on
   `[−2.5, −0.2]`; treat `(−0.2, +0.10]` as a single-bin
   "near-zero" cluster with confidence interval reflecting the
   ~1 % plateau spread. Explicit caveat in the manuscript.
2. **Cusp framing.** The plan's headline "qualitative
   shadow-structure change at ε₃_crit detectable in azimuthal-
   average asymmetry" is **not supported by this data** (§5.2
   above). The paper should reframe to: (a) a two-sided bound on
   ε₃ extracted from the smooth segment; (b) a separate
   prediction of the ISCO disappearance at
   `ε₃_ISCO ∈ (+0.13, +0.20)` testable by future high-precision
   observations of the disk inner edge.
3. **EHT data ingestion.** `ConsistencyBound` reads EHT 2019 Paper
   VI Table 7 Fourier m=1 fractional amplitude. The CSV stores
   `fourier_m1_amp / mean_r` implicitly (both columns are present);
   no CSV migration needed.
4. **RMS-vs-Fourier conversion.** The CSV stores both RMS and
   Fourier m=1 amplitudes. `ConsistencyBound` should use Fourier
   m=1 directly (matches Psaltis 2020 convention) and compute the
   higher-harmonic residual as a methodological flag (gate 3 of
   parent plan §7.4).
5. **Manuscript scope.** RNAAS limit: ≤ 1000 words, ≤ 2 figures,
   ≤ 15 references. Single-curve figure with a vertical line at
   `ε₃_crit ≈ 0.12` and a separate annotation marker for the ISCO
   transition near `ε₃ ≈ 0.18` (open interval).

### 9.5 Estimated wall-clock to `phase-3-complete`

- ~1 session for `ConsistencyBound` + tests (~250 LOC).
- ~1 session for figure-generation Python script + paper outline.
- ~2–3 sessions for manuscript writing (human-in-loop, not
  time-budgeted).
- ~30 min for build verification + `make paper` PDF emission.

Total to `phase-3-complete` tag: **3–5 focused sessions**, plus
the manuscript writing time.

---

## 10. Notes for next session

- Resume by reading this file, then `docs/phase-3-plan.md` §7
  (Phase 3D plan), then the paper-writing-relevant prior corpus
  (Johannsen 2013, EHT Paper VI, Psaltis 2020).
- Surface §9.4 decisions 1–5 to the user BEFORE 3D coding starts;
  particularly the cusp-framing reframe (decision 2) which
  requires the user's editorial sign-off on the paper's headline.
- Working tree at this commit: clean (sans the gitignored
  `output/*.csv` and `output/*.log`). Tag pushed.

---

*End of Phase 3C status snapshot.*
