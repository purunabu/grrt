package com.pranav.grrt.renderer;

/**
 * Phase 1 shader: 1.0 for escaped rays, 0.0 for everything else (captured,
 * max-steps, or NaN). Produces a flat-background image with the black-hole
 * shadow visible as a dark disk, suitable for the shadow-radius validation
 * milestone.
 *
 * <p>Thread-safe (stateless).
 */
public final class BinaryShader implements RayShader {

    @Override
    public float shade(RayOutcome outcome, double[] endState) {
        return outcome == RayOutcome.ESCAPED ? 1.0f : 0.0f;
    }
}
