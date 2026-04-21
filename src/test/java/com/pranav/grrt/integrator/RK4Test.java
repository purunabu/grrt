package com.pranav.grrt.integrator;

import com.pranav.grrt.metric.MinkowskiMetric;
import com.pranav.grrt.metric.Metric;
import com.pranav.grrt.metric.SchwarzschildMetric;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class RK4Test {

    // ---------------------------------------------------------------
    // Minkowski: null geodesics are straight lines in Cartesian coords.
    // Verify by integrating a radial outgoing photon and checking that
    // r advances linearly and θ, φ stay constant.
    // ---------------------------------------------------------------

    @Test
    void radialPhotonInMinkowskiIsStraightLine() {
        Metric m = new MinkowskiMetric();
        RK4 rk = new RK4();

        // Start at r=10, equator, moving radially outward at c=1.
        // Null condition: -(k^t)² + (k^r)² = 0  →  k^t = k^r = 1.
        double[] y = { 0.0, 10.0, Math.PI / 2, 0.0,  1.0, 1.0, 0.0, 0.0 };
        double[] next = new double[8];

        double h = 0.01;
        int n = 1000;
        for (int i = 0; i < n; i++) {
            rk.step(m, y, h, next);
            System.arraycopy(next, 0, y, 0, 8);
        }

        // After λ = 10, expect t = 10, r = 20, θ unchanged, φ unchanged.
        assertEquals(10.0, y[0], 1e-10, "t");
        assertEquals(20.0, y[1], 1e-10, "r");
        assertEquals(Math.PI / 2, y[2], 1e-10, "θ");
        assertEquals(0.0, y[3], 1e-10, "φ");
        // Momentum should be preserved.
        assertEquals(1.0, y[4], 1e-10, "k^t");
        assertEquals(1.0, y[5], 1e-10, "k^r");
    }

    @Test
    void transversePhotonInMinkowskiFollowsStraightLine() {
        // Photon passing with impact parameter b at x=b, moving in +x direction
        // in Cartesian. At closest approach: r=b, θ=π/2, φ=π/2, k_cart = (1,0,0).
        // In spherical: k^r = sin-component toward r; at this instant, purely tangential.
        // Equivalently: start at r=b, φ=π/2, purely in +x: that means k^φ carries the motion.
        //
        // Simpler setup: start at (r=b, θ=π/2, φ=0), and choose initial k such that
        // the trajectory is the line y = b in Cartesian: x = λ, y = b, z = 0.
        // Then r(λ) = √(λ² + b²), φ(λ) = atan2(b, λ).
        Metric m = new MinkowskiMetric();
        RK4 rk = new RK4();

        double b = 5.0;
        // At λ=0: x=0, y=b, z=0  →  r=b, θ=π/2, φ=π/2.
        // Cartesian velocity (1,0,0) in spherical basis at (r,θ=π/2,φ=π/2):
        //   k^r = sinθ cosφ · vx = 0
        //   k^θ = cosθ cosφ /r · vx = 0
        //   k^φ = -sinφ/(r sinθ) · vx = -1/b
        double[] y = { 0.0, b, Math.PI / 2, Math.PI / 2,
                       1.0, 0.0, 0.0, -1.0 / b };
        double[] next = new double[8];

        double h = 0.005;
        int n = 2000;
        double lambda = 0.0;
        for (int i = 0; i < n; i++) {
            rk.step(m, y, h, next);
            System.arraycopy(next, 0, y, 0, 8);
            lambda += h;
        }

        // Expected Cartesian position: (λ, b, 0).
        double expectedX = lambda;
        double expectedR = Math.sqrt(lambda * lambda + b * b);
        double expectedPhi = Math.atan2(b, lambda);

        assertEquals(expectedR, y[1], 1e-6, "r");
        assertEquals(Math.PI / 2, y[2], 1e-8, "θ stays equatorial");
        assertEquals(expectedPhi, y[3], 1e-6, "φ");

        // Cross-check in Cartesian.
        double xCart = y[1] * Math.sin(y[2]) * Math.cos(y[3]);
        double yCart = y[1] * Math.sin(y[2]) * Math.sin(y[3]);
        double zCart = y[1] * Math.cos(y[2]);
        assertEquals(expectedX, xCart, 1e-6, "x_cart");
        assertEquals(b,         yCart, 1e-6, "y_cart");
        assertEquals(0.0,       zCart, 1e-8, "z_cart");
    }

    // ---------------------------------------------------------------
    // Schwarzschild: unstable circular photon orbit at r = 3M.
    // A ray launched exactly tangentially at r=3M should remain (to the
    // integrator's accuracy) at r=3M for at least a few orbits before
    // the instability amplifies numerical noise.
    // ---------------------------------------------------------------

    @Test
    void schwarzschildCircularPhotonOrbitPreservedShortTerm() {
        double M = 1.0;
        Metric m = new SchwarzschildMetric(M);
        RK4 rk = new RK4();

        double r = 3.0 * M;
        double kt = 1.0;
        double kph = 1.0 / (3.0 * Math.sqrt(3.0) * M);
        double[] y = { 0.0, r, Math.PI / 2, 0.0,  kt, 0.0, 0.0, kph };
        double[] next = new double[8];

        // One orbital period in λ: Δφ = 2π at rate kph → Δλ = 2π / kph ≈ 2π·3√3 M ≈ 32.6.
        double period = 2.0 * Math.PI / kph;
        double h = period / 10_000.0;
        int nHalf = 5_000;   // half-orbit

        double maxDev = 0.0;
        for (int i = 0; i < nHalf; i++) {
            rk.step(m, y, h, next);
            System.arraycopy(next, 0, y, 0, 8);
            maxDev = Math.max(maxDev, Math.abs(y[1] - r));
        }

        // Over half an orbit, radial drift should be tiny.
        assertTrue(maxDev < 1e-6, "radial drift too large: " + maxDev);
    }

    // ---------------------------------------------------------------
    // Null-condition drift monitor: g_μν k^μ k^ν should stay near zero.
    // ---------------------------------------------------------------

    @Test
    void nullConditionPreservedInSchwarzschild() {
        Metric m = new SchwarzschildMetric(1.0);
        RK4 rk = new RK4();

        // Radial infalling photon from r=20.
        // -(1-2M/r)(k^t)² + (k^r)²/(1-2M/r) = 0  →  k^r = ±(1-2M/r) k^t.
        // We choose ingoing (k^r < 0), so after setting k^t=1: k^r = -(1-2/20) = -0.9.
        double r0 = 20.0;
        double f = 1.0 - 2.0 / r0;
        double[] y = { 0.0, r0, Math.PI / 2, 0.0,  1.0, -f, 0.0, 0.0 };
        double[] next = new double[8];

        assertEquals(0.0, m.nullNorm(sliceX(y), sliceK(y)), 1e-14, "initial not null");

        double h = 0.01;
        int n = 1500;  // stop before reaching horizon
        double maxNorm = 0.0;
        for (int i = 0; i < n; i++) {
            rk.step(m, y, h, next);
            System.arraycopy(next, 0, y, 0, 8);
            if (y[1] < 2.2) break;   // safety; stop well outside horizon
            maxNorm = Math.max(maxNorm, Math.abs(m.nullNorm(sliceX(y), sliceK(y))));
        }

        assertTrue(maxNorm < 1e-8, "null-norm drift too large: " + maxNorm);
    }

    // ---------------------------------------------------------------
    // Input validation
    // ---------------------------------------------------------------

    @Test
    void rejectsWrongLengthStates() {
        RK4 rk = new RK4();
        Metric m = new MinkowskiMetric();
        assertThrows(IllegalArgumentException.class,
                () -> rk.step(m, new double[7], 0.1, new double[8]));
        assertThrows(IllegalArgumentException.class,
                () -> rk.step(m, new double[8], 0.1, new double[7]));
    }

    // helpers
    private static double[] sliceX(double[] y) { return new double[]{y[0], y[1], y[2], y[3]}; }
    private static double[] sliceK(double[] y) { return new double[]{y[4], y[5], y[6], y[7]}; }
}
