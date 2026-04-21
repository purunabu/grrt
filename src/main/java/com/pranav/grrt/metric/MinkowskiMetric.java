package com.pranav.grrt.metric;

/**
 * Flat Minkowski spacetime in spherical coordinates (t, r, θ, φ).
 *
 * <p>Line element:
 * <pre>
 *   ds² = -dt² + dr² + r² dθ² + r² sin²θ dφ²
 * </pre>
 *
 * <p>Used as a sanity check for the integrator: null geodesics must be
 * straight lines in Cartesian coordinates, though they exhibit coordinate
 * "acceleration" in spherical coordinates due to the curvilinear chart.
 *
 * <p>Intended for testing only; no physical mass, no horizon.
 */
public final class MinkowskiMetric implements Metric {

    @Override public double mass() { return 0.0; }

    @Override public double horizonRadius() { return 0.0; }

    @Override
    public double[][] g(double[] x) {
        double r = x[1];
        double s = Math.sin(x[2]);
        double[][] g = new double[4][4];
        g[0][0] = -1.0;
        g[1][1] =  1.0;
        g[2][2] =  r * r;
        g[3][3] =  r * r * s * s;
        return g;
    }

    @Override
    public double[][] gInv(double[] x) {
        double r = x[1];
        double s = Math.sin(x[2]);
        double[][] gi = new double[4][4];
        gi[0][0] = -1.0;
        gi[1][1] =  1.0;
        gi[2][2] =  1.0 / (r * r);
        gi[3][3] =  1.0 / (r * r * s * s);
        return gi;
    }

    @Override
    public double[][][] christoffel(double[] x) {
        double r = x[1];
        double theta = x[2];
        double s = Math.sin(theta);
        double c = Math.cos(theta);

        double[][][] G = new double[4][4][4];
        // Γ^r_θθ = -r
        G[1][2][2] = -r;
        // Γ^r_φφ = -r sin²θ
        G[1][3][3] = -r * s * s;
        // Γ^θ_rθ = Γ^θ_θr = 1/r
        G[2][1][2] = 1.0 / r;
        G[2][2][1] = 1.0 / r;
        // Γ^θ_φφ = -sinθ cosθ
        G[2][3][3] = -s * c;
        // Γ^φ_rφ = Γ^φ_φr = 1/r
        G[3][1][3] = 1.0 / r;
        G[3][3][1] = 1.0 / r;
        // Γ^φ_θφ = Γ^φ_φθ = cot θ
        G[3][2][3] = c / s;
        G[3][3][2] = c / s;
        return G;
    }

    @Override
    public void geodesicAcceleration(double[] x, double[] k, double[] out) {
        double r = x[1];
        double theta = x[2];
        double s = Math.sin(theta);
        double c = Math.cos(theta);

        double kr = k[1];
        double kth = k[2];
        double kph = k[3];

        out[0] = 0.0;
        out[1] = r * kth * kth + r * s * s * kph * kph;
        out[2] = -2.0 / r * kr * kth + s * c * kph * kph;
        out[3] = -2.0 / r * kr * kph - 2.0 * c / s * kth * kph;
    }
}
