package com.pranav.grrt.metric;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Phase 3A validation gates 1 through 4 for {@link JohannsenPsaltisMetric}
 * (see {@code docs/phase-3-plan.md} §4.7). Gates 5 (null-norm drift on a
 * tracer orbit) and 6 (closed-form equatorial photon-orbit prediction)
 * are deferred to a separate commit; they require the DP45 integrator
 * and an external reference table.
 *
 * <p>Gate coverage here:
 * <ul>
 *   <li>(1) {@code g · gInv = I} on a (r, θ, ε₃) grid at algebraic precision.</li>
 *   <li>(2) Christoffel symmetry in lower indices.</li>
 *   <li>(3) Kerr reduction at ε₃ = 0 vs {@link KerrMetric} componentwise.</li>
 *   <li>(4) Asymptotic flatness: the JP deformation vanishes at large r.</li>
 * </ul>
 */
class JohannsenPsaltisMetricTest {

    private static final double ALG_TOL      = 1e-12;  // algebraic identities
    private static final double KERR_TOL     = 1e-12;  // ε₃ → 0 reduction
    private static final double ASYMPT_TOL   = 1e-6;   // deformation at r = 1000 M

    private static double[] pos(double r, double theta) {
        return new double[] { 0.0, r, theta, 0.0 };
    }

    // ------------------------------------------------------------------
    // Constructor validation
    // ------------------------------------------------------------------

    @Test
    void constructorRejectsNonPositiveMass() {
        assertThrows(IllegalArgumentException.class,
                () -> new JohannsenPsaltisMetric(0.0, 0.5, 0.5));
        assertThrows(IllegalArgumentException.class,
                () -> new JohannsenPsaltisMetric(-1.0, 0.5, 0.5));
        assertThrows(IllegalArgumentException.class,
                () -> new JohannsenPsaltisMetric(Double.NaN, 0.5, 0.5));
        assertThrows(IllegalArgumentException.class,
                () -> new JohannsenPsaltisMetric(Double.POSITIVE_INFINITY, 0.5, 0.5));
    }

    @Test
    void constructorRejectsNonFiniteSpin() {
        assertThrows(IllegalArgumentException.class,
                () -> new JohannsenPsaltisMetric(1.0, Double.NaN, 0.5));
        assertThrows(IllegalArgumentException.class,
                () -> new JohannsenPsaltisMetric(1.0, Double.POSITIVE_INFINITY, 0.5));
    }

    @Test
    void constructorRejectsExtremalAndSuperextremalSpin() {
        assertThrows(IllegalArgumentException.class,
                () -> new JohannsenPsaltisMetric(1.0, 1.0, 0.5));
        assertThrows(IllegalArgumentException.class,
                () -> new JohannsenPsaltisMetric(1.0, -1.0, 0.5));
        assertThrows(IllegalArgumentException.class,
                () -> new JohannsenPsaltisMetric(1.0, 1.5, 0.5));
    }

    @Test
    void constructorRejectsNonFiniteEpsilon3() {
        assertThrows(IllegalArgumentException.class,
                () -> new JohannsenPsaltisMetric(1.0, 0.5, Double.NaN));
        assertThrows(IllegalArgumentException.class,
                () -> new JohannsenPsaltisMetric(1.0, 0.5, Double.POSITIVE_INFINITY));
        assertThrows(IllegalArgumentException.class,
                () -> new JohannsenPsaltisMetric(1.0, 0.5, Double.NEGATIVE_INFINITY));
    }

    @Test
    void constructorRejectsEpsilon3ThatCausesCoordinatePathology() {
        // For a=0.9, M=1, r+ = 1 + sqrt(0.19) ≈ 1.43589, r+³ ≈ 2.9605.
        // 1 + ε₃ · 1/r+³ ≤ 0 when ε₃ ≤ -r+³ ≈ -2.9605.
        assertThrows(IllegalArgumentException.class,
                () -> new JohannsenPsaltisMetric(1.0, 0.9, -3.0));
        assertThrows(IllegalArgumentException.class,
                () -> new JohannsenPsaltisMetric(1.0, 0.9, -5.0));
    }

    @Test
    void constructorAcceptsApprovedSweepRange() {
        assertDoesNotThrow(() -> new JohannsenPsaltisMetric(1.0, 0.9, -2.0));
        assertDoesNotThrow(() -> new JohannsenPsaltisMetric(1.0, 0.9, -1.0));
        assertDoesNotThrow(() -> new JohannsenPsaltisMetric(1.0, 0.9,  0.0));
        assertDoesNotThrow(() -> new JohannsenPsaltisMetric(1.0, 0.9,  1.0));
        assertDoesNotThrow(() -> new JohannsenPsaltisMetric(1.0, 0.9,  2.0));
    }

    // ------------------------------------------------------------------
    // Accessors and landmarks
    // ------------------------------------------------------------------

    @Test
    void massSpinAndEpsilonAccessors() {
        JohannsenPsaltisMetric jp = new JohannsenPsaltisMetric(2.0, 1.2, 0.5);
        assertEquals(2.0, jp.mass(),     ALG_TOL);
        assertEquals(1.2, jp.spin(),     ALG_TOL);
        assertEquals(0.5, jp.epsilon3(), ALG_TOL);
    }

    @Test
    void unitMassConstructorDefaultsMassToOne() {
        JohannsenPsaltisMetric jp = new JohannsenPsaltisMetric(0.9, 0.5);
        assertEquals(1.0, jp.mass(),     ALG_TOL);
        assertEquals(0.9, jp.spin(),     ALG_TOL);
        assertEquals(0.5, jp.epsilon3(), ALG_TOL);
    }

    @Test
    void horizonRadiusMatchesKerrPolarAxisValue() {
        JohannsenPsaltisMetric jp = new JohannsenPsaltisMetric(1.0, 0.9, 1.0);
        assertEquals(1.0 + Math.sqrt(0.19), jp.horizonRadius(), ALG_TOL);
    }

    // ------------------------------------------------------------------
    // Gate 1: g · gInv = I on a (r, θ, ε₃) grid
    // ------------------------------------------------------------------

    @Test
    void metricInverseIdentityOnParameterGrid() {
        // 80 points: 4 r × 4 θ × 5 ε₃ at a = 0.9 M. Covers equator,
        // off-equator, horizon-adjacent and asymptotic radii.
        double[] rValues   = { 3.0, 8.0, 20.0, 100.0 };
        double[] thetaVals = { Math.PI / 6.0, Math.PI / 4.0, Math.PI / 3.0, Math.PI / 2.0 };
        double[] epsVals   = { -1.0, -0.25, 0.0, 0.25, 1.0 };

        for (double eps : epsVals) {
            JohannsenPsaltisMetric jp = new JohannsenPsaltisMetric(1.0, 0.9, eps);
            for (double r : rValues) {
                for (double th : thetaVals) {
                    assertGtimesGinvIsIdentity(jp, r, th, ALG_TOL);
                }
            }
        }
    }

    @Test
    void metricInverseIdentityAtModerateSpinAndExtremeEps3() {
        // Second grid: a = 0.5, retrograde a = -0.5, and a different M.
        double[][] cases = {
                // M,    a,    ε₃,   r,    θ
                { 1.0,  0.5, +2.0, 10.0, Math.PI / 4.0 },
                { 1.0,  0.5, -2.0, 10.0, Math.PI / 2.0 },
                { 1.0, -0.5,  1.5,  5.0, 1.0 },
                { 2.0,  1.2,  0.3,  8.0, Math.PI / 3.0 },
        };
        for (double[] c : cases) {
            JohannsenPsaltisMetric jp = new JohannsenPsaltisMetric(c[0], c[1], c[2]);
            assertGtimesGinvIsIdentity(jp, c[3], c[4], ALG_TOL);
        }
    }

    private static void assertGtimesGinvIsIdentity(
            Metric m, double r, double theta, double tol) {
        double[] x = pos(r, theta);
        double[][] g  = m.g(x);
        double[][] gi = m.gInv(x);
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                double sum = 0.0;
                for (int k = 0; k < 4; k++) sum += g[i][k] * gi[k][j];
                double expected = (i == j) ? 1.0 : 0.0;
                assertEquals(expected, sum, tol,
                        "g·gInv[%d][%d] at (r=%g, θ=%g)".formatted(i, j, r, theta));
            }
        }
    }

    // ------------------------------------------------------------------
    // Gate 2: Christoffel symmetric in lower indices
    // ------------------------------------------------------------------

    @Test
    void christoffelSymmetricInLowerIndicesAcrossCases() {
        double[][] cases = {
                // a,     r,    θ,           ε₃
                {  0.5,   5.0,  0.8,          0.3 },
                {  0.9,   3.0,  Math.PI / 2, -0.5 },
                { -0.5,   8.0,  1.2,          1.0 },
                {  0.0,  10.0,  0.5,          0.0 },  // Schwarzschild × no deformation
                {  0.9,  15.0,  Math.PI / 4,  2.0 },
                {  0.9,   2.0,  1.0,         -1.5 },  // near-horizon at large |ε₃|
        };
        for (double[] c : cases) {
            JohannsenPsaltisMetric jp = new JohannsenPsaltisMetric(1.0, c[0], c[3]);
            double[][][] G = jp.christoffel(pos(c[1], c[2]));
            for (int alpha = 0; alpha < 4; alpha++) {
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        assertEquals(G[alpha][i][j], G[alpha][j][i], 1e-14,
                                "Γ asymmetric at (a=%g, r=%g, θ=%g, ε₃=%g) [%d][%d][%d]"
                                        .formatted(c[0], c[1], c[2], c[3], alpha, i, j));
                    }
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // Gate 3: Kerr reduction at ε₃ = 0
    // ------------------------------------------------------------------

    @Test
    void zeroEpsilonMetricMatchesKerr() {
        JohannsenPsaltisMetric jp = new JohannsenPsaltisMetric(1.0, 0.9, 0.0);
        KerrMetric k = new KerrMetric(1.0, 0.9);
        double[][] points = {
                { 10.0, Math.PI / 4.0 },
                {  4.0, Math.PI / 3.0 },
                {  2.5, 1.2 },
                { 50.0, Math.PI / 2.0 },
                {100.0, 0.5 },
        };
        for (double[] p : points) {
            double[] x = pos(p[0], p[1]);
            double[][] gJp  = jp.g(x);
            double[][] gK   = k.g(x);
            double[][] giJp = jp.gInv(x);
            double[][] giK  = k.gInv(x);
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    assertEquals(gK[i][j], gJp[i][j], KERR_TOL,
                            "g[%d][%d] at (r=%g, θ=%g)".formatted(i, j, p[0], p[1]));
                    assertEquals(giK[i][j], giJp[i][j], KERR_TOL,
                            "gInv[%d][%d] at (r=%g, θ=%g)".formatted(i, j, p[0], p[1]));
                }
            }
        }
    }

    @Test
    void zeroEpsilonChristoffelMatchesKerr() {
        JohannsenPsaltisMetric jp = new JohannsenPsaltisMetric(1.0, 0.9, 0.0);
        KerrMetric k = new KerrMetric(1.0, 0.9);
        double[][] points = {
                { 10.0, Math.PI / 4.0 },
                {  4.0, Math.PI / 3.0 },
                {  2.5, 1.2 },
                {100.0, 1.0 },
                {  3.0, Math.PI / 2.0 },
        };
        double worst = 0.0;
        for (double[] p : points) {
            double[] x = pos(p[0], p[1]);
            double[][][] gJp = jp.christoffel(x);
            double[][][] gK  = k.christoffel(x);
            for (int alpha = 0; alpha < 4; alpha++) {
                for (int i = 0; i < 4; i++) {
                    for (int j = 0; j < 4; j++) {
                        double d = Math.abs(gJp[alpha][i][j] - gK[alpha][i][j]);
                        if (d > worst) worst = d;
                        assertEquals(gK[alpha][i][j], gJp[alpha][i][j], KERR_TOL,
                                "Γ^%d_%d%d at (r=%g, θ=%g)"
                                        .formatted(alpha, i, j, p[0], p[1]));
                    }
                }
            }
        }
        System.out.println("[3A gate 3] max |Γ_JP(ε₃=0) - Γ_Kerr| over points = " + worst);
    }

    @Test
    void zeroEpsilonGeodesicAccelerationMatchesKerr() {
        JohannsenPsaltisMetric jp = new JohannsenPsaltisMetric(1.0, 0.9, 0.0);
        KerrMetric k = new KerrMetric(1.0, 0.9);
        double[][] trials = {
                { 5.0,  1.1,          0.3,  1.2, -0.4,  0.05, 0.07 },
                { 10.0, Math.PI / 4,  0.0,  1.0,  0.1, -0.2,  0.5  },
                { 3.0,  Math.PI / 2,  1.0,  1.5, -0.2,  0.0,  0.3  },
        };
        double[] aJp = new double[4];
        double[] aK  = new double[4];
        for (double[] t : trials) {
            double[] x  = { 0.0, t[0], t[1], t[2] };
            double[] km = { t[3], t[4], t[5], t[6] };
            jp.geodesicAcceleration(x, km, aJp);
            k.geodesicAcceleration(x,  km, aK);
            for (int mu = 0; mu < 4; mu++) {
                assertEquals(aK[mu], aJp[mu], KERR_TOL,
                        "μ=%d at (r=%g, θ=%g)".formatted(mu, t[0], t[1]));
            }
        }
    }

    @Test
    void kerrReductionResidualIsFloatingPointNoise() {
        // Gate 3 sanity: the Kerr-reduction residual at ε₃ = 0 should scale
        // like machine epsilon, with no systematic structure. If the CSE
        // block had a transcription residual, it would show up as a
        // per-(M, a) pattern that breaks this scaling. One print line per
        // pair for review.
        double[][] pairs = {
                // M,   a
                { 1.0,  0.5 },
                { 2.0,  0.9 },   // a=0.9 at M=2: r+ = 2 + sqrt(4-0.81)=3.786
                { 1.0, -0.9 },   // retrograde at a=0.9
        };
        double[][] points = {
                { 10.0, Math.PI / 4.0 },
                {  4.0, Math.PI / 3.0 },
                {  3.0, Math.PI / 2.0 },
                { 50.0, 0.8 },
        };
        for (double[] pr : pairs) {
            double M = pr[0];
            double a = pr[1];
            JohannsenPsaltisMetric jp = new JohannsenPsaltisMetric(M, a, 0.0);
            KerrMetric k = new KerrMetric(M, a);
            double maxDg = 0.0, maxDgi = 0.0, maxDch = 0.0, maxGmag = 0.0;
            for (double[] p : points) {
                // Scale r by M so test points stay outside r+.
                double[] x = pos(p[0] * M, p[1]);
                double[][]   gJ  = jp.g(x);    double[][]   gK  = k.g(x);
                double[][]   giJ = jp.gInv(x); double[][]   giK = k.gInv(x);
                double[][][] chJ = jp.christoffel(x); double[][][] chK = k.christoffel(x);
                for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) {
                    maxDg   = Math.max(maxDg,   Math.abs(gJ[i][j]  - gK[i][j]));
                    maxDgi  = Math.max(maxDgi,  Math.abs(giJ[i][j] - giK[i][j]));
                    maxGmag = Math.max(maxGmag, Math.abs(gK[i][j]));
                }
                for (int alpha = 0; alpha < 4; alpha++)
                    for (int i = 0; i < 4; i++)
                        for (int j = 0; j < 4; j++)
                            maxDch = Math.max(maxDch,
                                    Math.abs(chJ[alpha][i][j] - chK[alpha][i][j]));
            }
            System.out.printf(
                    "[3A gate 3 sanity] (M=%.1f, a=%+.1f): "
                    + "max|Δg|=%.2e, max|ΔgInv|=%.2e, max|ΔΓ|=%.2e, |g|_max=%.3f%n",
                    M, a, maxDg, maxDgi, maxDch, maxGmag);
            assertTrue(maxDg  < KERR_TOL,
                    "g residual above tolerance for (M=%g, a=%g): %g".formatted(M, a, maxDg));
            assertTrue(maxDgi < KERR_TOL,
                    "gInv residual above tolerance for (M=%g, a=%g): %g".formatted(M, a, maxDgi));
            assertTrue(maxDch < KERR_TOL,
                    "Γ residual above tolerance for (M=%g, a=%g): %g".formatted(M, a, maxDch));
        }
    }

    @Test
    void zeroEpsilonHorizonMatchesKerr() {
        JohannsenPsaltisMetric jp = new JohannsenPsaltisMetric(1.0, 0.9, 0.0);
        KerrMetric k = new KerrMetric(1.0, 0.9);
        assertEquals(k.horizonRadius(), jp.horizonRadius(), ALG_TOL);
    }

    @Test
    void iscoBracketPositionDiagnostic() {
        // Logs the walking-bracket rLo that bisection starts from, for
        // prograde and retrograde ISCO at a=0.9, ε₃=0, along with the
        // Kerr photon-sphere radii (Bardeen 1972 closed form) and the
        // final r_ISCO residual vs Kerr. Bisection convergence is
        // independent of bracket start once the bracket is valid, so
        // these numbers confirm whether any observed r_ISCO residual
        // scales with (rLo − r_photon) — the "bracket-position"
        // hypothesis — or with 1/|d²E/dr²(r_ISCO)| — the "FD-error"
        // hypothesis. No hard assertion: this is a diagnostic print for
        // review.
        JohannsenPsaltisMetric jp = new JohannsenPsaltisMetric(1.0, 0.9, 0.0);
        KerrMetric k = new KerrMetric(1.0, 0.9);

        // Kerr equatorial photon-orbit radii, Bardeen 1972 eq. 2.18.
        double a = 0.9;
        double rPhPro   = 2.0 * (1.0 + Math.cos((2.0 / 3.0) * Math.acos(-a)));
        double rPhRetro = 2.0 * (1.0 + Math.cos((2.0 / 3.0) * Math.acos(+a)));

        double rLoPro   = jp.iscoBracketLowerBound(true);
        double rLoRetro = jp.iscoBracketLowerBound(false);

        double iscoJpPro   = jp.iscoRadius(true);
        double iscoKPro    = k.iscoRadius(true);
        double iscoJpRetro = jp.iscoRadius(false);
        double iscoKRetro  = k.iscoRadius(false);

        System.out.printf(
                "[3A ISCO diag] pro:   rLo=%.6f  r_ph_Kerr=%.6f  gap=%.4f  "
                + "r_isco_JP=%.12f  r_isco_K=%.12f  diff=%+.3e%n",
                rLoPro, rPhPro, rLoPro - rPhPro,
                iscoJpPro, iscoKPro, iscoJpPro - iscoKPro);
        System.out.printf(
                "[3A ISCO diag] retro: rLo=%.6f  r_ph_Kerr=%.6f  gap=%.4f  "
                + "r_isco_JP=%.12f  r_isco_K=%.12f  diff=%+.3e%n",
                rLoRetro, rPhRetro, rLoRetro - rPhRetro,
                iscoJpRetro, iscoKRetro, iscoJpRetro - iscoKRetro);
    }

    @Test
    void zeroEpsilonIscoMatchesKerrBardeenPressTeukolsky() {
        JohannsenPsaltisMetric jp = new JohannsenPsaltisMetric(1.0, 0.9, 0.0);
        KerrMetric k = new KerrMetric(1.0, 0.9);
        // JP uses numeric bisection on dE/dr; Kerr uses the Bardeen closed form.
        // At ε₃ = 0 they solve the same physics; agreement is limited only
        // by the finite-difference step + bisection tolerance (~1e-9).
        assertEquals(k.iscoRadius(true),  jp.iscoRadius(true),  1e-8,
                "prograde ISCO, a=0.9, ε₃=0");
        assertEquals(k.iscoRadius(false), jp.iscoRadius(false), 1e-8,
                "retrograde ISCO, a=0.9, ε₃=0");
    }

    @Test
    void zeroSpinZeroEpsilonMatchesSchwarzschild() {
        // Double reduction: JP(a=0, ε₃=0) → Kerr(a=0) → Schwarzschild.
        JohannsenPsaltisMetric jp = new JohannsenPsaltisMetric(1.0, 0.0, 0.0);
        SchwarzschildMetric s = new SchwarzschildMetric(1.0);

        double[][] points = {
                { 10.0, Math.PI / 4.0 },
                {  4.0, Math.PI / 3.0 },
                {  2.5, 1.2 },
        };
        for (double[] p : points) {
            double[] x = pos(p[0], p[1]);
            double[][] gJp = jp.g(x);
            double[][] gS  = s.g(x);
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    assertEquals(gS[i][j], gJp[i][j], KERR_TOL,
                            "g[%d][%d] at (r=%g, θ=%g)".formatted(i, j, p[0], p[1]));
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // Gate 4: asymptotic flatness — JP deformation vanishes at large r
    // ------------------------------------------------------------------

    @Test
    void deformationVanishesAtR1000M() {
        // |h| ≈ |ε₃| / r³ → ~1e-9 for ε₃ = 1 at r = 1000 M.
        // The JP vs Kerr component difference should be well below 1e-6.
        JohannsenPsaltisMetric jpEps = new JohannsenPsaltisMetric(1.0, 0.9, 1.0);
        KerrMetric             kerr  = new KerrMetric(1.0, 0.9);
        double[] thetas = { Math.PI / 6.0, Math.PI / 3.0, Math.PI / 2.0 };
        double worst = 0.0;
        for (double th : thetas) {
            double[] x = pos(1000.0, th);
            double[][] gJp = jpEps.g(x);
            double[][] gK  = kerr.g(x);
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    double diff = Math.abs(gJp[i][j] - gK[i][j]);
                    if (diff > worst) worst = diff;
                    assertTrue(diff < ASYMPT_TOL,
                            "|g_JP - g_Kerr|[%d][%d] at (r=1000, θ=%g) = %g"
                                    .formatted(i, j, th, diff));
                }
            }
        }
        System.out.println("[3A gate 4] max |g_JP(ε₃=1) - g_Kerr| at r=1000M = " + worst);
    }

    @Test
    void asymptoticMinkowskiInSphericalFormAtVeryLargeR() {
        // At r = 1e6 M, both 2M/r (~2e-6) and h (~1e-18) are negligible.
        // g should recover the Minkowski-in-spherical form
        //   diag(-1, 1, Σ, Σ sin²θ)  with Σ ≈ r².
        JohannsenPsaltisMetric jp = new JohannsenPsaltisMetric(1.0, 0.9, 1.0);
        double r = 1e6;
        double theta = Math.PI / 3.0;
        double[] x = pos(r, theta);
        double[][] g = jp.g(x);
        double s = Math.sin(theta);
        double s2 = s * s;

        assertEquals(-1.0, g[0][0], 1e-5, "g_tt → -1");
        assertEquals(+1.0, g[1][1], 1e-5, "g_rr → 1");
        // g_θθ and g_φφ scale as r²; assert relative error.
        assertEquals(1.0, g[2][2] / (r * r),       1e-10, "g_θθ / r² → 1");
        assertEquals(1.0, g[3][3] / (r * r * s2),  1e-10, "g_φφ / (r² sin²θ) → 1");
        assertEquals(0.0, g[0][3], 1e-5, "g_tφ → 0");
    }
}
