package com.pranav.grrt.integrator;

import com.pranav.grrt.metric.Metric;

/**
 * ODE integrator for null geodesics in a curved spacetime.
 *
 * <p>State vector layout (length 8):
 * <pre>
 *   y[0..3] = x^μ = (t, r, θ, φ)
 *   y[4..7] = k^μ = dx^μ/dλ
 * </pre>
 *
 * <p>The affine parameter λ is dimensionless in geometrized units; for a
 * photon, integrating toward the past corresponds to negative step sizes
 * (or equivalently, flipped initial k^μ).
 */
public interface Integrator {

    /**
     * Advance the state by one step of size h.
     *
     * @param metric spacetime metric providing geodesic RHS
     * @param y      current state, length 8; not modified
     * @param h      affine-parameter step; may be negative
     * @param out    output state, length 8; overwritten
     */
    void step(Metric metric, double[] y, double h, double[] out);
}
