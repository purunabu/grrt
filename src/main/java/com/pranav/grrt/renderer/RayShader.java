package com.pranav.grrt.renderer;

/**
 * Maps a terminated null geodesic to a scalar pixel value. Phase 1 uses the
 * binary {@link BinaryShader}; Phase 2 will swap in a radiative-transfer
 * shader that accumulates emissivity along the traced ray.
 *
 * <p>Implementations must be thread-safe (pure functions or immutable
 * state). Renderer calls this from a parallel stream.
 */
@FunctionalInterface
public interface RayShader {

    /**
     * @param outcome  terminal classification of the ray
     * @param endState final 8-state [t, r, θ, φ, k^t, k^r, k^θ, k^φ] at
     *                 termination; read-only for shader purposes, do not mutate
     * @return pixel value
     */
    float shade(RayOutcome outcome, double[] endState);
}
