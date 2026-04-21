package com.pranav.grrt.integrator;

import com.pranav.grrt.metric.Metric;

/**
 * Classical 4th-order Runge-Kutta integrator for null geodesics.
 * Fixed step size. Local truncation error O(h⁵); global error O(h⁴).
 *
 * <p>Thread-safety: a single {@code RK4} instance holds scratch buffers and
 * is <b>not</b> thread-safe. The renderer should construct one instance per
 * worker thread (e.g. via {@link ThreadLocal} or by creating a fresh integrator
 * per pixel-batch).
 */
public final class RK4 implements Integrator {

    private final double[] k1 = new double[8];
    private final double[] k2 = new double[8];
    private final double[] k3 = new double[8];
    private final double[] k4 = new double[8];
    private final double[] tmp = new double[8];
    private final double[] accel = new double[4];

    @Override
    public void step(Metric metric, double[] y, double h, double[] out) {
        if (y.length != 8 || out.length != 8) {
            throw new IllegalArgumentException("state vectors must have length 8");
        }

        rhs(metric, y, k1);

        axpy(y, k1, 0.5 * h, tmp);
        rhs(metric, tmp, k2);

        axpy(y, k2, 0.5 * h, tmp);
        rhs(metric, tmp, k3);

        axpy(y, k3, h, tmp);
        rhs(metric, tmp, k4);

        double h6 = h / 6.0;
        for (int i = 0; i < 8; i++) {
            out[i] = y[i] + h6 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
    }

    /**
     * Geodesic RHS: writes f(y) into {@code dy}.
     * f[0..3] = k^μ (the momentum components of y)
     * f[4..7] = -Γ^μ_αβ k^α k^β
     */
    private void rhs(Metric metric, double[] y, double[] dy) {
        // Position derivatives = momentum.
        dy[0] = y[4];
        dy[1] = y[5];
        dy[2] = y[6];
        dy[3] = y[7];

        // Momentum derivatives = geodesic acceleration.
        // geodesicAcceleration reads only x and k; writes accel[0..3].
        metric.geodesicAcceleration(
                /* x = */ sliceX(y),
                /* k = */ sliceK(y),
                accel);
        dy[4] = accel[0];
        dy[5] = accel[1];
        dy[6] = accel[2];
        dy[7] = accel[3];
    }

    /** out = a + s·b, element-wise over length 8. */
    private static void axpy(double[] a, double[] b, double s, double[] out) {
        for (int i = 0; i < 8; i++) out[i] = a[i] + s * b[i];
    }

    /**
     * View of y[0..3] as a length-4 array.
     * Uses a scratch buffer to avoid allocation in the hot loop.
     */
    private final double[] xBuf = new double[4];
    private double[] sliceX(double[] y) {
        xBuf[0] = y[0]; xBuf[1] = y[1]; xBuf[2] = y[2]; xBuf[3] = y[3];
        return xBuf;
    }

    private final double[] kBuf = new double[4];
    private double[] sliceK(double[] y) {
        kBuf[0] = y[4]; kBuf[1] = y[5]; kBuf[2] = y[6]; kBuf[3] = y[7];
        return kBuf;
    }
}
