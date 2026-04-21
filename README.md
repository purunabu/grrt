# grrt

A general-relativistic ray tracer in Java 21 that integrates null geodesics
in curved spacetime to image black-hole shadows; Phase 1 ships Schwarzschild,
with the end goal a Research Notes of the AAS paper on photon-ring asymmetry
vs. a Johannsen-Psaltis deviation from Kerr, compared to EHT M87* data.

## Phase 1 result

![Schwarzschild shadow, 256x256, r_obs = 1000 M, fov = 0.02 rad](output/first_shadow.png)

Phase 1 validation: at r_obs = 1000 M, equatorial inclination, fov = 0.02 rad,
256x256, the shadow edge on the middle row is bracketed by pixels 193
(captured) and 194 (escaped). Midpoint **193.5 px** vs analytic prediction
**193.94 px** — Δ = **−0.44 px**. In geometrized units the measured shadow
radius is 5.16 M against the textbook `b_crit = 3√3 M ≈ 5.196 M`, about 0.7 %
low. Two-zone RK4 (`h = 0.01 M` for `r ≤ 10 M`, `h = 0.5 M` beyond), horizon
cushion `0.01 M`, ~15 s on an Apple M3 (8 cores).

## Architecture

Four interfaces, intentionally small — the `Metric` is a drop-in replacement
so Kerr and Johannsen-Psaltis land as single-file extensions without touching
the integrator, camera, or renderer.

- **`Metric`** (`src/main/java/com/pranav/grrt/metric/`) —
  `g(x)`, `gInv(x)`, closed-form Christoffels, and a sparsity-exploiting
  geodesic acceleration. Ships `SchwarzschildMetric` and `MinkowskiMetric`
  (the latter only for integrator validation).
- **`Integrator`** (`integrator/`) — advances the 8-state
  `[t, r, θ, φ, kᵗ, kʳ, kᶿ, kᵠ]` by one step in affine parameter. Phase 1
  ships fixed-step RK4; Dormand-Prince-45 arrives in Phase 2.
- **`Camera`** (`camera/`) — static-observer orthonormal tetrad, pixel
  `(i, j)` → initial null 4-momentum. Future-directed `k^μ` pointing outward
  from the observer; forward integration traces the reversed photon path.
- **`Renderer`** (`renderer/`) — parallel pixel loop
  (`IntStream.parallel()` with `ThreadLocal<RK4>`), terminates rays on the
  horizon cushion or escape radius, runs a `RayShader` on the outcome, and
  writes FITS via `io/FitsWriter`. Phase 1 ships `BinaryShader`
  (1 = escaped, 0 = everything else).

See [CLAUDE.md](CLAUDE.md) for physics conventions (signature, units,
state-vector layout, sign conventions) and the project-wide coding rules.

## Build & run

Prereqs: **Java 21**, **Maven 3.9+**.

```
mvn test        # 30 tests, ~17 s; includes 256x256 shadow validation
mvn package     # jar in target/grrt-0.1.0-SNAPSHOT.jar
mvn clean       # wipe target/
```

Phase 1 does not ship a CLI driver. For a worked
`Camera → RK4 → Renderer` example that renders a 256² image and asserts the
shadow radius against `3√3 M`, see
[`RendererTest.shadowRadiusAt256MatchesCriticalImpactParameter`](src/test/java/com/pranav/grrt/renderer/RendererTest.java).
`Renderer.renderToFits(Path, String)` is the one-call entry point once you
have a `Camera`, `Metric`, and `RenderConfig` wired up —
`RenderConfig.defaults(rObs)` gives the same settings used for the image
above.

## Status

Phase 1 — Schwarzschild:

- [x] `Metric` interface, `SchwarzschildMetric`, `MinkowskiMetric`
- [x] `Integrator` interface + fixed-step `RK4`
- [x] `Camera` (static-observer tetrad, pixel → null 4-momentum)
- [x] `Renderer` (parallel pixel loop, two-zone RK4, binary shader)
- [x] `FitsWriter` (linear WCS, provenance header)
- [x] Validation: shadow radius vs `3√3 M` within 0.5 px at 256²

Coming in Phase 2:

- [ ] `KerrMetric` (Boyer-Lindquist, closed-form Christoffels)
- [ ] Adaptive `DormandPrince45` integrator
- [ ] Disk emission model + radiative-transfer `RayShader`
- [ ] Observer aberration / Doppler if the error budget demands it

Phase 3:

- [ ] `JohannsenPsaltisMetric`
- [ ] EHT M87* comparison pipeline
- [ ] Research Notes of the AAS paper

## See also

- [CLAUDE.md](CLAUDE.md) — physics conventions, coding rules, testing
  requirements, and guidance for agent contributors. Read this before
  editing anything in `src/main/java/com/pranav/grrt/`.
