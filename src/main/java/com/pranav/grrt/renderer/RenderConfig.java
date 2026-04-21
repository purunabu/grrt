package com.pranav.grrt.renderer;

/**
 * Immutable configuration for a {@link Renderer} pass.
 *
 * @param horizonCushion  additive cushion ε_h (units of M): rays with
 *                        r &lt; horizonRadius + ε_h are classified CAPTURED.
 *                        Schwarzschild coordinates are singular at r = 2M,
 *                        so we cannot integrate there; ε_h = 0.01 M keeps
 *                        1/(1 − 2M/r) finite while leaving the shadow
 *                        boundary (set by the photon sphere at r = 3M)
 *                        untouched.
 * @param escapeRadius    radial threshold above which a ray is classified
 *                        ESCAPED. Phase 1 uses 2 · r_obs: the observer is
 *                        in the asymptotic regime (M/r_obs ~ 10⁻³ for
 *                        M87*-style imaging), so once a ray passes the
 *                        observer's radius outgoing, further deflection
 *                        is negligible.
 * @param hNear           affine-parameter step size (units of M) for
 *                        r ≤ rTransition. 0.01 M is the step validated by
 *                        RK4Test for null-norm drift &lt; 10⁻⁸ over ~1500
 *                        steps near the BH.
 * @param hFar            affine step for r &gt; rTransition. 0.5 M; the far
 *                        zone is nearly flat, so larger steps are safe.
 * @param rTransition     boundary between near and far zones, in units of M.
 *                        Default 10 M — well outside the photon sphere at
 *                        3 M.
 * @param maxSteps        safety cap per ray. 10⁶ suffices for any physically
 *                        meaningful Phase 1 geodesic; rays that hit this
 *                        cap are near-photon-sphere tangents and are
 *                        shaded as captured.
 * @param shader          maps terminated rays to pixel values
 */
public record RenderConfig(
        double horizonCushion,
        double escapeRadius,
        double hNear,
        double hFar,
        double rTransition,
        int maxSteps,
        RayShader shader
) {

    public RenderConfig {
        if (!(horizonCushion > 0.0))
            throw new IllegalArgumentException("horizonCushion must be > 0: " + horizonCushion);
        if (!(escapeRadius > 0.0))
            throw new IllegalArgumentException("escapeRadius must be > 0: " + escapeRadius);
        if (!(hNear > 0.0 && hFar > 0.0))
            throw new IllegalArgumentException("step sizes must be positive");
        if (!(rTransition > 0.0))
            throw new IllegalArgumentException("rTransition must be > 0: " + rTransition);
        if (maxSteps <= 0)
            throw new IllegalArgumentException("maxSteps must be > 0: " + maxSteps);
        if (shader == null)
            throw new IllegalArgumentException("shader must not be null");
    }

    /**
     * Phase 1 defaults for Schwarzschild at the given observer radius.
     * Cushion 0.01 M, escape radius 2·r_obs, two-zone step size
     * (h_near=0.01 M, h_far=0.5 M, transition 10 M), maxSteps 10⁶,
     * {@link BinaryShader}.
     */
    public static RenderConfig defaults(double rObs) {
        return new RenderConfig(
                0.01, 2.0 * rObs,
                0.01, 0.5, 10.0,
                1_000_000,
                new BinaryShader());
    }

    /** Human-readable summary suitable for the FITS HPOLICY header card. */
    public String stepPolicyString() {
        return String.format(
                "h=%.3g r<=%.1fM, h=%.3g r>%.1fM",
                hNear, rTransition, hFar, rTransition);
    }
}
