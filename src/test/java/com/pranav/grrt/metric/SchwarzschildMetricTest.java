package com.pranav.grrt.metric;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class SchwarzschildMetricTest {

    private static final double TOL = 1e-12;

    private static double[] pos(double r, double theta) {
        return new double[] { 0.0, r, theta, 0.0 };
    }

    // ------------------------------------------------------------------
    // Constructor + trivial accessors
    // ------------------------------------------------------------------

    @Test
    void constructorRejectsNonPositiveMass() {
        assertThrows(IllegalArgumentException.class, () -> new SchwarzschildMetric(0.0));
        assertThrows(IllegalArgumentException.class, () -> new SchwarzschildMetric(-1.0));
        assertThrows(IllegalArgumentException.class, () -> new SchwarzschildMetric(Double.NaN));
        assertThrows(IllegalArgumentException.class, () -> new SchwarzschildMetric(Double.POSITIVE_INFINITY));
    }

    @Test
    void massAndHorizon() {
        Metric m = new SchwarzschildMetric(2.5);
        assertEquals(2.5, m.mass(), TOL);
        assertEquals(5.0, m.horizonRadius(), TOL);
    }

    // ------------------------------------------------------------------
    // Metric structure
    // ------------------------------------------------------------------

    @Test
    void metricIsDiagonalAndMatchesClosedForm() {
        Metric m = new SchwarzschildMetric(1.0);
        double r = 10.0;
        double theta = Math.PI / 3.0;
        double[][] g = m.g(pos(r, theta));

        double f = 1.0 - 2.0 / r;
        assertEquals(-f, g[0][0], TOL);
        assertEquals(1.0 / f, g[1][1], TOL);
        assertEquals(r * r, g[2][2], TOL);
        assertEquals(r * r * Math.sin(theta) * Math.sin(theta), g[3][3], TOL);

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                if (i != j) assertEquals(0.0, g[i][j], TOL);
            }
        }
    }

    @Test
    void metricInverseCheck() {
        Metric m = new SchwarzschildMetric(1.0);
        double[] x = pos(8.0, 1.2);
        double[][] g = m.g(x);
        double[][] gi = m.gInv(x);

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                double sum = 0.0;
                for (int k = 0; k < 4; k++) sum += g[i][k] * gi[k][j];
                double expected = (i == j) ? 1.0 : 0.0;
                assertEquals(expected, sum, 1e-10, "g * gInv [%d][%d]".formatted(i, j));
            }
        }
    }

    @Test
    void asymptoticFlatness() {
        // At r → ∞: g_tt → -1, g_rr → +1 (θ, φ components scale with r, as expected in spherical coords).
        Metric m = new SchwarzschildMetric(1.0);
        double[][] g = m.g(pos(1e10, Math.PI / 2));
        assertEquals(-1.0, g[0][0], 1e-9);
        assertEquals(1.0, g[1][1], 1e-9);
    }

    // ------------------------------------------------------------------
    // Christoffel symbols
    // ------------------------------------------------------------------

    @Test
    void christoffelSymmetricInLowerIndices() {
        Metric m = new SchwarzschildMetric(1.0);
        double[][][] G = m.christoffel(pos(7.0, 0.7));
        for (int a = 0; a < 4; a++) {
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    assertEquals(G[a][i][j], G[a][j][i], TOL,
                            "Γ^%d_%d%d asymmetric".formatted(a, i, j));
                }
            }
        }
    }

    @Test
    void christoffelClosedFormAtSpecificPoint() {
        double M = 1.0;
        double r = 6.0;
        double theta = Math.PI / 4.0;
        Metric m = new SchwarzschildMetric(M);
        double[][][] G = m.christoffel(pos(r, theta));

        double f = 1.0 - 2.0 * M / r;
        double s = Math.sin(theta);
        double c = Math.cos(theta);

        assertEquals(M / (r * r * f),       G[0][0][1], TOL); // Γ^t_tr
        assertEquals(M * f / (r * r),       G[1][0][0], TOL); // Γ^r_tt
        assertEquals(-M / (r * r * f),      G[1][1][1], TOL); // Γ^r_rr
        assertEquals(-(r - 2.0 * M),        G[1][2][2], TOL); // Γ^r_θθ
        assertEquals(-(r - 2.0 * M) * s * s, G[1][3][3], TOL); // Γ^r_φφ
        assertEquals(1.0 / r,               G[2][1][2], TOL); // Γ^θ_rθ
        assertEquals(-s * c,                G[2][3][3], TOL); // Γ^θ_φφ
        assertEquals(1.0 / r,               G[3][1][3], TOL); // Γ^φ_rφ
        assertEquals(c / s,                 G[3][2][3], TOL); // Γ^φ_θφ
    }

    // ------------------------------------------------------------------
    // Physics: circular photon orbit at r = 3M
    // ------------------------------------------------------------------

    @Test
    void circularPhotonOrbitHasZeroAcceleration() {
        // At r=3M, θ=π/2, with k^r = k^θ = 0 and k^φ/k^t = 1/(3√3 M),
        // the geodesic has a^μ = 0 instantaneously (circular photon orbit,
        // unstable but stationary).
        double M = 1.0;
        Metric m = new SchwarzschildMetric(M);
        double r = 3.0 * M;
        double kt = 1.0;
        double kph = 1.0 / (3.0 * Math.sqrt(3.0) * M);
        double[] x = new double[] { 0.0, r, Math.PI / 2.0, 0.0 };
        double[] k = new double[] { kt, 0.0, 0.0, kph };

        // Null condition sanity
        assertEquals(0.0, m.nullNorm(x, k), 1e-14, "photon not null");

        double[] a = new double[4];
        m.geodesicAcceleration(x, k, a);
        for (int mu = 0; mu < 4; mu++) {
            assertEquals(0.0, a[mu], 1e-14, "a[%d] != 0".formatted(mu));
        }
    }

    // ------------------------------------------------------------------
    // Optimized vs default geodesic acceleration
    // ------------------------------------------------------------------

    @Test
    void optimizedMatchesDefaultAcceleration() {
        // Override acts as ground truth for correctness of the default
        // implementation (which uses christoffel()). Compare on a non-trivial
        // state: off-equator, non-zero radial and angular momentum.
        Metric m = new SchwarzschildMetric(1.0);
        double[] x = new double[] { 0.0, 8.0, 1.1, 0.3 };
        double[] k = new double[] { 1.2, -0.4, 0.05, 0.07 };

        double[] aOpt = new double[4];
        m.geodesicAcceleration(x, k, aOpt);

        double[] aDefault = new double[4];
        Metric wrapper = new Metric() {
            @Override public double mass() { return m.mass(); }
            @Override public double horizonRadius() { return m.horizonRadius(); }
            @Override public double[][] g(double[] x) { return m.g(x); }
            @Override public double[][] gInv(double[] x) { return m.gInv(x); }
            @Override public double[][][] christoffel(double[] x) { return m.christoffel(x); }
            // Intentionally does NOT override geodesicAcceleration → uses default.
        };
        wrapper.geodesicAcceleration(x, k, aDefault);

        for (int mu = 0; mu < 4; mu++) {
            assertEquals(aDefault[mu], aOpt[mu], 1e-12,
                    "mismatch at index %d".formatted(mu));
        }
    }
}
