package com.pranav.grrt.metric;

/**
 * Schwarzschild metric for a non-rotating, uncharged black hole of mass M
 * in Schwarzschild coordinates (t, r, θ, φ).
 *
 * <p>Line element:
 * <pre>
 *   ds² = -(1 - 2M/r) dt² + (1 - 2M/r)⁻¹ dr² + r² dθ² + r² sin²θ dφ²
 * </pre>
 *
 * <p>Landmarks: event horizon at r = 2M, photon sphere at r = 3M, shadow
 * impact parameter at infinity b_crit = 3√3 M ≈ 5.196 M.
 *
 * <p>Valid for r &gt; 2M and θ ∈ (0, π). Behaviour at r = 2M is a coordinate
 * singularity; callers must terminate integration before reaching the horizon.
 * Behaviour at θ = 0, π is a coordinate singularity (cot θ diverges).
 */
public final class SchwarzschildMetric implements Metric {

    private final double M;

    /**
     * @param mass gravitational mass M in geometrized units (G = c = 1);
     *             must be positive and finite
     */
    public SchwarzschildMetric(double mass) {
        if (mass <= 0.0 || !Double.isFinite(mass)) {
            throw new IllegalArgumentException("Mass must be positive and finite: " + mass);
        }
        this.M = mass;
    }

    /** Unit-mass black hole (M = 1 in geometrized units). */
    public SchwarzschildMetric() {
        this(1.0);
    }

    @Override
    public double mass() {
        return M;
    }

    @Override
    public double horizonRadius() {
        return 2.0 * M;
    }

    @Override
    public double[][] g(double[] x) {
        double r = x[1];
        double sinTheta = Math.sin(x[2]);
        double f = 1.0 - 2.0 * M / r;

        double[][] g = new double[4][4];
        g[0][0] = -f;
        g[1][1] = 1.0 / f;
        g[2][2] = r * r;
        g[3][3] = r * r * sinTheta * sinTheta;
        return g;
    }

    @Override
    public double[][] gInv(double[] x) {
        double r = x[1];
        double sinTheta = Math.sin(x[2]);
        double f = 1.0 - 2.0 * M / r;

        double[][] gi = new double[4][4];
        gi[0][0] = -1.0 / f;
        gi[1][1] = f;
        gi[2][2] = 1.0 / (r * r);
        gi[3][3] = 1.0 / (r * r * sinTheta * sinTheta);
        return gi;
    }

    @Override
    public double[][][] christoffel(double[] x) {
        double r = x[1];
        double theta = x[2];
        double f = 1.0 - 2.0 * M / r;
        double sinTheta = Math.sin(theta);
        double cosTheta = Math.cos(theta);

        double rr = r * r;
        double rrf = rr * f;
        double twoM = 2.0 * M;

        double[][][] G = new double[4][4][4];

        // Γ^t_tr = Γ^t_rt = M / (r² f)
        double Gt_tr = M / rrf;
        G[0][0][1] = Gt_tr;
        G[0][1][0] = Gt_tr;

        // Γ^r_tt = M f / r²
        G[1][0][0] = M * f / rr;

        // Γ^r_rr = -M / (r² f)
        G[1][1][1] = -M / rrf;

        // Γ^r_θθ = -(r - 2M)
        G[1][2][2] = -(r - twoM);

        // Γ^r_φφ = -(r - 2M) sin²θ
        G[1][3][3] = -(r - twoM) * sinTheta * sinTheta;

        // Γ^θ_rθ = Γ^θ_θr = 1/r
        double invR = 1.0 / r;
        G[2][1][2] = invR;
        G[2][2][1] = invR;

        // Γ^θ_φφ = -sinθ cosθ
        G[2][3][3] = -sinTheta * cosTheta;

        // Γ^φ_rφ = Γ^φ_φr = 1/r
        G[3][1][3] = invR;
        G[3][3][1] = invR;

        // Γ^φ_θφ = Γ^φ_φθ = cot θ
        double cotTheta = cosTheta / sinTheta;
        G[3][2][3] = cotTheta;
        G[3][3][2] = cotTheta;

        return G;
    }

    /**
     * Optimized geodesic acceleration exploiting the sparsity of Schwarzschild
     * Christoffels. Avoids allocating the 4x4x4 tensor; ~10x faster than the
     * default implementation.
     */
    @Override
    public void geodesicAcceleration(double[] x, double[] k, double[] out) {
        double r = x[1];
        double theta = x[2];
        double f = 1.0 - 2.0 * M / r;
        double sinTheta = Math.sin(theta);
        double cosTheta = Math.cos(theta);

        double kt = k[0];
        double kr = k[1];
        double kth = k[2];
        double kph = k[3];

        double rr = r * r;
        double rrf = rr * f;
        double twoM = 2.0 * M;
        double rMinus2M = r - twoM;

        // a^t = -2 Γ^t_tr k^t k^r
        out[0] = -2.0 * M / rrf * kt * kr;

        // a^r = -Γ^r_tt (k^t)² - Γ^r_rr (k^r)² - Γ^r_θθ (k^θ)² - Γ^r_φφ (k^φ)²
        out[1] = -(M * f / rr) * kt * kt
               + (M / rrf) * kr * kr
               + rMinus2M * kth * kth
               + rMinus2M * sinTheta * sinTheta * kph * kph;

        // a^θ = -2 Γ^θ_rθ k^r k^θ - Γ^θ_φφ (k^φ)²
        out[2] = -2.0 / r * kr * kth + sinTheta * cosTheta * kph * kph;

        // a^φ = -2 Γ^φ_rφ k^r k^φ - 2 Γ^φ_θφ k^θ k^φ
        out[3] = -2.0 / r * kr * kph - 2.0 * cosTheta / sinTheta * kth * kph;
    }
}
