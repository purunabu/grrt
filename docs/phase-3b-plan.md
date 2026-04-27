# Phase 3B Sub-plan — Novikov-Thorne disk emission

Companion to `docs/phase-3-plan.md` §5. Closes the Phase-2-deferred
"disk emission model" checkbox and lands the renderer plumbing
required for the Phase 3C ε₃ sweep.

Parent commit: `f779db1` (`phase-3a-complete`). Sub-phase 3B lands as
three commits leading to tag `phase-3b-nt-disk`. Approved 2026-04-27
following the 3B-readiness review in `docs/phase-3a-status.md` §6.

---

## 1. Scope

Sub-phase 3B delivers four additions and one interface extension:

1. New `com.pranav.grrt.disk` package: `Disk` interface,
   `NovikovThorneDisk`, `DiskEmissionShader`.
2. Dense-output extension to `DormandPrince45` for sub-step
   evaluation at disk equator crossings.
3. Disk-crossing detection in `AdaptiveRayTracer`, gated on a new
   `RayOutcome.HIT_DISK` variant.
4. Position-dependent horizon termination via a new default method
   on the `Metric` interface (`isInsideHorizon`), with a JP-specific
   override that uses `Δ + a² h sin²θ ≤ tol`.
5. New `scripts/page_thorne_reference.py` regenerable Python
   reference for the gate-1 NT flux table, mirroring the
   `jp_photon_orbit_reference.py` pattern from 3A gate 6.

The `Metric` interface gains exactly one method
(`isInsideHorizon`), default-implemented r-based, overridden in
`JohannsenPsaltisMetric` only. Kerr and Schwarzschild inherit the
default — Phase-2 tests stay bit-exact.

---

## 2. Decisions resolved at 3B kickoff

The four decisions left open by `docs/phase-3a-status.md` §§6.3–6.4
are resolved as follows. Each was approved by the user prior to
coding.

### 2.1 DP45 dense-output interpolant: Shampine 5th-order

`DormandPrince45.interpolate(double θ, double[] out)` evaluates a
polynomial in `θ ∈ [0, 1]` over the most recently accepted step,
producing a continuous extension matching the discrete propagation
at the endpoints and globally O(h⁵).

Reference: Hairer, Nørsett & Wanner (1993), *Solving Ordinary
Differential Equations I*, 2nd ed., Eq. II.6.20 (continuous
extension for DP5(4)7M). The 5th-order extension uses the seven
existing `k₁..k₇` plus one additional RHS evaluation `k_dense` at
an intermediate `θ`. The extra `k_dense` is computed lazily on first
`interpolate` call within a step and cached for subsequent calls
within the same step (e.g. multi-iteration bisection during disk
crossing detection).

Cost: ~14% extra `f`-evals per step that bisects. Steps that do not
bisect pay zero. Steps with multiple bisection iterations reuse the
cached `k_dense` for free.

Rejected alternatives:

- **4th-order natural Hermite from `(y, y′)` at endpoints.** Zero
  extra `f`-evals, but gate 5 (1e-8 mid-step accuracy) sits at the
  edge of order-4 capability with `atol=rtol=1e-10`. Bumping into
  the floor would force tolerance relaxation, violating CLAUDE.md
  "never relax a test tolerance to make a test pass". Held as a
  fallback only if 5th-order shows pathological behavior in
  practice.
- **Cubic Hermite.** O(h³) global; needs `h ≲ 0.01 M` to hit 1e-8.
  At long-integration default `h ~ 1 M` it fails by ~5 orders of
  magnitude. Rejected.

If gate 5 fails with the 5th-order extension, fall back to the
4th-order Hermite **before** relaxing tolerance.

### 2.2 Page-Thorne 1974 reference flux: Python script + paste

`scripts/page_thorne_reference.py` numerically integrates
Page-Thorne 1974 eq. (15n) using `scipy.integrate.quad` over a
Sympy-lambdified Kerr equatorial circular-orbit integrand. Output
format: whitespace-delimited with header
`# a r_over_M F_pageThorne` covering 6 radii × 2 spins (a ∈ {0, 0.9})
= 12 reference values.

Output is pasted by hand into `NovikovThorneDiskTest` as a
`double[][]` literal, frozen against the script's output at the time
of paste. Script remains in the repo and is rerunnable; if
`NovikovThorneDisk.surfaceFlux` ever changes form, the developer
regenerates the table and updates the test in lockstep.

Rejected alternatives:

- **Inline computation in the test.** For `a = 0.9`, the integrand
  needs `r_ISCO`, `Ω(r)`, `E(r)`, `L(r)` — exactly the production
  code being tested. Self-referential; defeats the gate's purpose.
- **Tabulate from Page 1974 directly.** Paper prints to 3
  significant figures; gate-1 tolerance of 1e-6 is unachievable
  against rounded printed values.

This mirrors the pattern validated for 3A gate 6
(`jp_photon_orbit_reference.py`).

### 2.3 Disk inner edge: per-frame `Disk` construction, cached `r_ISCO`

`NovikovThorneDisk(Metric metric, double rOuter, double iscoCushion)`
calls `metric.iscoRadius(true)` once in the constructor and stores
the result in a `final double rIsco` field. `Renderer` constructs
one `Disk` instance per frame. The pixel loop reuses the same
`Disk`, so `r_ISCO` is computed exactly once per frame regardless of
pixel count.

For the Phase 3C sweep, `EpsilonSweep` rebuilds
`JohannsenPsaltisMetric(0.9, εᵢ)` per frame and constructs a fresh
`NovikovThorneDisk(metric, 20.0, 0.001)` against it. Each frame uses
its own `r_ISCO(JP(0.9, εᵢ))`, automatically. No cross-frame
contamination, no per-pixel recomputation.

Cost: one ISCO bisection per frame (~1 ms). Free.

A `NovikovThorneDiskTest` smoke test asserts the constructor does
not throw for `ε₃ ∈ {−2.5, −1.0, −0.5, +0.05, +0.10}` at `a = 0.9`,
covering the smooth-regime endpoints + interior of the §6.5 sweep
grid.

### 2.4 JP equatorial-horizon termination: `Metric.isInsideHorizon` (interface extension)

The current `AdaptiveRayTracer` terminates when
`y[1] < horizonRadius() + cushion`. For `JohannsenPsaltisMetric`
with `ε₃ > ε₃_crit ≈ 0.12` at `a = 0.9`, the equatorial horizon
locus `Δ + a² h sin²θ = 0` sits inside the polar-axis Kerr horizon,
so `horizonRadius()` (which returns the polar-axis value as the
safe outer scalar per `JohannsenPsaltisMetric` Javadoc) gives a
**too-conservative** termination radius on the equator: a photon
can pass through what the polar-axis test calls "the horizon"
without the metric being singular on its trajectory. For
`ε₃ < 0`, the situation reverses — the JP equatorial horizon sits
*outside* the Kerr horizon, and the polar-axis fallback would let
rays drift below the singular surface.

**Resolution:** add a default method to the `Metric` interface:

```java
/**
 * @return true iff (x[1], x[2]) is at or inside the metric's event
 *         horizon, within tolerance {@code tol}. Default: r-based,
 *         {@code x[1] - horizonRadius() < tol}, valid for any metric
 *         whose horizon is θ-independent. Position-dependent metrics
 *         (e.g. Johannsen-Psaltis at ε₃ ≠ 0) override.
 */
default boolean isInsideHorizon(double[] x, double tol) {
    return x[1] - horizonRadius() < tol;
}
```

Overridden in `JohannsenPsaltisMetric` to test
`Δ(r) + a² h(r,θ) sin²θ ≤ tol`. `KerrMetric` and `SchwarzschildMetric`
inherit the default — their `horizonRadius()` is θ-independent, so
the test reduces exactly to the prior `r ≤ rCapture` check. Phase-2
DP45 + RK4 + renderer tests stay bit-exact (gate 6).

`AdaptiveRayTracer` constructor caches `horizonCushion` (not
`rCapture`); the per-step termination becomes
`metric.isInsideHorizon(y, horizonCushion)`.

**This is a `Metric` interface change.** Per CLAUDE.md "ask when
touching Metric", explicit user approval was sought and granted at
3B kickoff. The change adds one default method, breaks no existing
implementations, and is the cleaner alternative to an `instanceof`
branch in the renderer (which Phase 1 was explicit about avoiding).

CLAUDE.md will gain a brief note in the same commit that introduces
the method, mirroring the documentation pattern for `horizonRadius`
and `iscoRadius` in 3A.

**Sweep-grid implication:** with position-dependent termination in
place, `EpsilonSweep` can render the full 15-point grid from
`docs/phase-3-plan.md` §6.5 — including the post-cusp transition
points `{+0.13, +0.20, +1.00}` — correctly. The bound extraction in
3D still uses only `ε₃ < ε₃_crit` per `phase-3-plan.md` §7.3; the
post-cusp points populate the figure as a separate curve segment.

---

## 3. Sub-phase decomposition (3 commits, 1 tag)

The split is chosen so the regression boundary falls cleanly:
existing renderer / integrator tests stay bit-exact through 3B.1.

### 3.1 Commit 3B.1 — Disk core + DP45 dense output

Self-contained additions. Renderer untouched.

**Files created:**

| File | Est. lines | Purpose |
|---|---|---|
| `src/main/java/com/pranav/grrt/disk/Disk.java` | ~50 | Interface: `crossedEquator`, `keplerianFourVelocity`, `temperature`, `rIsco`, `rOuter` |
| `src/main/java/com/pranav/grrt/disk/NovikovThorneDisk.java` | ~280 | Page-Thorne 1974 flux via metric-agnostic `Metric.g(x)` and `Metric.iscoRadius` |
| `src/main/java/com/pranav/grrt/disk/DiskEmissionShader.java` | ~140 | `RayShader` impl: `B_bol(T(r_emit)) · g⁴` |
| `src/test/java/com/pranav/grrt/disk/NovikovThorneDiskTest.java` | ~220 | Gates 1, 2, 4 + ε₃ smoke |
| `src/test/java/com/pranav/grrt/disk/DiskEmissionShaderTest.java` | ~110 | Redshift unit + face-on annulus prep |
| `scripts/page_thorne_reference.py` | ~100 | Sympy + scipy NT flux table generator |

**Files modified:**

| File | Modification |
|---|---|
| `src/main/java/com/pranav/grrt/integrator/DormandPrince45.java` | Add `interpolate(double θ, double[] out)`; persist `yPrev`, `hAccepted`, lazy `kDense` cache. `step()` body and signature unchanged. |

**Test budget:** `mvn test` < 15 s. No rendering in 3B.1.

### 3.2 Commit 3B.2 — Renderer integration

**Files modified:**

| File | Modification |
|---|---|
| `src/main/java/com/pranav/grrt/metric/Metric.java` | Add `default boolean isInsideHorizon(double[] x, double tol)` |
| `src/main/java/com/pranav/grrt/metric/JohannsenPsaltisMetric.java` | Override `isInsideHorizon` with `Δ + a² h sin²θ ≤ tol` |
| `src/main/java/com/pranav/grrt/renderer/AdaptiveRayTracer.java` | Add Disk-aware constructor overload; per-step sign test on `θ − π/2`; bisect via `interpolate`; switch termination to `metric.isInsideHorizon` |
| `src/main/java/com/pranav/grrt/renderer/RayOutcome.java` | Add `HIT_DISK` constant |
| `src/main/java/com/pranav/grrt/renderer/BinaryShader.java` and any other `RayOutcome` consumer | Audit and extend exhaustive switches; non-disk consumers map `HIT_DISK` to background per shader semantics |

**Test budget:** `mvn test` < 30 s. Includes one face-on Schwarzschild + NT render at 256² for gate 3.

If 3B.3 is small enough, fold it into this commit.

### 3.3 Commit 3B.3 — Position-dependent termination, all gates, tag

**Files modified:**

| File | Modification |
|---|---|
| `src/test/java/com/pranav/grrt/renderer/RendererTest.java` | Position-dependent termination smoke render at `JP(0.9, +0.20)` to confirm prograde shadow boundary lands on the equatorial JP horizon, not the polar-axis Kerr horizon |
| `CLAUDE.md` | Tick `Disk emission model` checkbox; document `isInsideHorizon` interface method |

After `mvn test` green, tag `phase-3b-nt-disk`.

---

## 4. Validation gates

Six gates per `docs/phase-3-plan.md` §5.5, refined for actual
implementation. All run in `mvn test`; total budget ≤ 30 s.

| # | Test | Tolerance | Lives in |
|---|---|---|---|
| 1 | NT radial flux `F(r)` vs Page-Thorne 1974 eq. (15n) at 6 radii × `a ∈ {0, 0.9}` | 1e-6 | `NovikovThorneDiskTest` |
| 2 | ISCO: `KerrMetric(a=0.9).iscoRadius(true)` → 2.32088 M (Bardeen 1972 Tab. 1) | 1e-6 | `NovikovThorneDiskTest` (sanity floor; uses existing `iscoRadius`) |
| 3 | Face-on Schwarzschild disk at 256² (`a=0`, `i=0°`): annulus is centroid-symmetric, outer edge at 20 M | 0.5 px | `DiskEmissionShaderTest` |
| 4 | Inner-ring redshift, `a=0.9`, `i=85°`: shader `g(r_ISCO)` matches `1/√(−g_tt^eq(r_ISCO))` | 1e-4 | `DiskEmissionShaderTest` |
| 5 | DP45 dense output at `θ = 0.5` mid-step vs DP45 with half-step direct on Kerr photon orbit | 1e-8 | `DormandPrince45Test` (new method) |
| 6 | All Phase-2-and-3A tests pass bit-exactly with `phase-3a-complete` (88 tests) | exact | regression run at end of 3B.2 and 3B.3 |

Gate 5 numerical detail: integrate a Kerr equatorial photon orbit
from `r₀ = 50 M`, `b = 5.5 M` for one step at `atol=rtol=1e-10`;
compare `interpolate(0.5, midState)` against direct DP45 propagation
of two half-steps at the same tolerance. The discrete-step error
floor is ~1e-10 per the existing
`DormandPrince45KerrTest.kerrDeflectingPhotonConservationUnderLongIntegration`
analysis; the interpolant is expected to add at most one order of
margin.

Gate 6 numerical detail: `mvn test` produces a per-class pass/fail
manifest. After 3B.2 and 3B.3, the manifest must show:

- Identical results for all classes that existed at 3A close-out
  (88 tests).
- Net new tests only in `disk/*Test`, `DormandPrince45Test`
  (interpolant), and one in `RendererTest` (3B.3 smoke).

---

## 5. Wall-clock estimates (M3, 10 threads)

| Stage | Wall-clock |
|---|---|
| 3B.1 build + test | ~10 s |
| 3B.2 build + test (adds gate 3 256² face-on render) | ~12 s |
| 3B.3 build + test (adds JP-horizon smoke render) | ~14 s |
| Gate 3: 256² face-on Schwarzschild + NT | ~6 s/frame |
| Gate 4: 256² edge-on Kerr + NT (per-pixel; sampled, not full render) | ~7 s |
| Phase 3C-implied 512² Kerr+NT at `i=17°` per frame | ~24 s |

Coding effort: ~3 focused sessions, ~6–10 hours session time, per
`docs/phase-3a-status.md` §6.5 estimate. This plan does not change
that estimate.

---

## 6. Failure modes

| Failure | Detection | Mitigation |
|---|---|---|
| Missed equator crossing when DP45 step straddles disk | Gate 3 (face-on annulus has no missing wedges) | Cap step when `\|θ − π/2\| < 0.1 rad`; sign-change check on both endpoints; `phase-3-plan.md` §5.7 |
| Plunging-region emission > 0 | Bisect lands `r < r_ISCO` → `RayOutcome.HORIZON` not `HIT_DISK`; no emission | `Disk.crossedEquator` guards against `r < rIsco`; explicit unit test |
| Page-Thorne integrand divergence near ISCO | Gate 1 row at `r = r_ISCO + 0.01 M` | Truncate the integral at exactly `r_ISCO`; document inability to render `r < r_ISCO` |
| DP45 dense output departs from discrete propagation | Gate 5 (1e-8 against half-step direct) | Shampine 5th has matching global order; if gate fails, fall back to 4th Hermite + tighten `atol` to compensate, **never** relax gate 5 |
| `RayOutcome.HIT_DISK` consumer break | Compile error at any non-exhaustive switch site | Audit `grep -r 'RayOutcome\\.'` before modifying enum; add `default ->` arms in every consumer |
| JP equatorial horizon termination wrong for `ε₃ > ε₃_crit` | 3B.3 smoke render at `JP(0.9, +0.20)` 256² — prograde shadow boundary at JP equatorial horizon, not Kerr polar-axis | Position-dependent `isInsideHorizon` per §2.4 |
| `interpolate(θ)` reused across step boundary (stale `kDense`) | Lazy cache invalidation: `kDense` recomputed when caller supplies a new step | Track step generation counter on integrator; `interpolate` checks generation before reuse |
| 3B interface change (`isInsideHorizon`) breaks downstream subclasses | Phase 2 + 3A `mvn test` (gate 6) | Default method body matches existing renderer behavior bit-exactly; KerrMetric/SchwarzschildMetric inherit unchanged |

---

## 7. Commit + tag structure

```
[3B.1] disk: add Novikov-Thorne disk model and DP45 dense output
[3B.2] renderer: plumb disk crossing into AdaptiveRayTracer with HIT_DISK
[3B.3] metric: position-dependent JP horizon termination + 3B gates green
       (followed by `git tag phase-3b-nt-disk`)
```

Author: configured git identity, no attribution trailers per
CLAUDE.md "Agent behavior expectations".

If 3B.3 ends up under ~100 lines, fold it into 3B.2.

---

## 8. Pre-coding checklist (status: all green)

| # | Item | Status |
|---|---|---|
| 1 | DP45 dense-output interpolant: 5th-order Shampine | ✅ §2.1 |
| 2 | Page-Thorne reference: Python script + paste pattern | ✅ §2.2 |
| 3 | Disk inner edge: per-frame `Disk` with cached `r_ISCO` | ✅ §2.3 |
| 4 | JP horizon termination: `Metric.isInsideHorizon` default + JP override | ✅ §2.4 |
| 5 | CLAUDE.md `Disk emission model` reads "Phase 3" | ✅ committed at `45abbd5` |
| 6 | CLAUDE.md `JohannsenPsaltisMetric` checkbox ticked | ✅ committed at `915cd66` |

---

## 9. Pre-3B reading order (for fresh sessions)

1. `docs/phase-3a-status.md` — Phase 3A close-out and 3B readiness
   check.
2. `docs/phase-3b-plan.md` — this file (the authoritative 3B plan).
3. `docs/phase-3-plan.md` §5 — parent plan, superseded by this
   sub-plan where they overlap.
4. `CLAUDE.md` — repo rules, especially "ask when touching Metric"
   and "no allocation in the integrator inner loop".

The `Metric.isInsideHorizon` interface change in §2.4 is the only
load-bearing-abstraction edit this sub-phase makes; it has been
explicitly approved.

---

## 10. Deferred / out of scope for 3B

- Frequency-channel emission (NT shader is bolometric only). Phase
  3C does not need it; the asymmetry metric is wavelength-agnostic
  in the optically-thick limit.
- Disk radiative transfer beyond face-on Planck — no scattering,
  no absorption along the line of sight to infinity. Justified by
  the optically-thick assumption.
- Sub-ISCO plunging emission. NT terminates at `r_ISCO` exactly.
- Adaptive-tolerance per-pixel control. Per-call `atol`/`rtol` stays
  uniform across pixels for 3B; deferred to 3C if needed.

*End of sub-plan.*
