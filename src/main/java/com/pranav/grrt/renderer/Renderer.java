package com.pranav.grrt.renderer;

import com.pranav.grrt.camera.Camera;
import com.pranav.grrt.integrator.Integrator;
import com.pranav.grrt.io.FitsWriter;
import com.pranav.grrt.metric.Metric;

import java.io.IOException;
import java.nio.file.Path;
import java.util.function.Supplier;
import java.util.stream.IntStream;

/**
 * Pixel loop + ray tracer. For each pixel the Camera provides an initial
 * 8-state, the integrator advances it in affine parameter, and the shader
 * assigns a float value based on the terminal outcome.
 *
 * <h2>Termination</h2>
 *
 * Per {@link RenderConfig}:
 * <ul>
 *   <li>r &lt; horizonRadius + horizonCushion ⇒ CAPTURED</li>
 *   <li>r &gt; escapeRadius ⇒ ESCAPED</li>
 *   <li>any NaN in the state ⇒ NAN (shaded as captured by Phase 1 shader)</li>
 *   <li>step count ≥ maxSteps ⇒ MAX_STEPS (likewise)</li>
 * </ul>
 *
 * <h2>Step size</h2>
 *
 * Two-zone adaptive: {@code h = hNear} when {@code r ≤ rTransition},
 * {@code h = hFar} otherwise. Zone is re-evaluated every step. This is the
 * cheapest useful adaptivity; real adaptive stepping (Dormand-Prince) is
 * deferred to Phase 2.
 *
 * <h2>Parallelism</h2>
 *
 * Pixels are processed in parallel via {@code IntStream.range(0, W*H).parallel()}.
 * Each worker uses its own {@link ThreadLocal} {@link Integrator}
 * (RK4 holds scratch buffers and is not thread-safe) plus thread-local
 * state buffers. {@link Camera} is shared (thread-safe). Because each
 * output pixel is written by exactly one worker and the shader is pure,
 * the result is deterministic regardless of thread scheduling.
 *
 * <h2>Runtime</h2>
 *
 * Estimated wall-clock on Apple M3 (8 cores, 3 GHz), Schwarzschild at
 * r_obs = 1000 M, fov ≈ 20 M/r_obs, two-zone steps:
 * <ul>
 *   <li>256×256  ≈ 16 s</li>
 *   <li>1024×1024 ≈ 4–5 min</li>
 * </ul>
 * Assumes ~5000 RK4 steps/ray average (much higher near the shadow edge
 * due to photon-sphere grazing) and ~500 FLOPs/step for RK4 + Schwarzschild
 * optimized acceleration + state update + termination checks, with Java
 * JIT and GC overhead included.
 */
public final class Renderer {

    private final Camera camera;
    private final Metric metric;
    private final Supplier<Integrator> integratorFactory;
    private final RenderConfig config;

    /**
     * @param camera             provides initial states; must be built against
     *                           {@code metric}
     * @param metric             spacetime metric used by the integrator
     * @param integratorFactory  produces a fresh {@link Integrator} for each
     *                           worker thread (RK4 is not thread-safe). Typical
     *                           value: {@code RK4::new}.
     * @param config             termination, step-size, shader settings
     */
    public Renderer(Camera camera, Metric metric,
                    Supplier<Integrator> integratorFactory,
                    RenderConfig config) {
        if (camera == null || metric == null || integratorFactory == null || config == null) {
            throw new IllegalArgumentException("null argument");
        }
        this.camera = camera;
        this.metric = metric;
        this.integratorFactory = integratorFactory;
        this.config = config;
    }

    /**
     * Render the full image. Returns a newly-allocated {@code float[H][W]} with
     * {@code image[0]} as the top row (j = 0 is top, FITS convention).
     * Pure — no internal state is mutated; two calls with identical configuration
     * produce byte-identical arrays.
     */
    public float[][] render() {
        final int W = camera.width();
        final int H = camera.height();
        final float[][] image = new float[H][W];

        // Per-thread state: one Integrator + two 8-element buffers per worker.
        final ThreadLocal<Integrator> tlInt = ThreadLocal.withInitial(integratorFactory);
        final ThreadLocal<double[]> tlY    = ThreadLocal.withInitial(() -> new double[8]);
        final ThreadLocal<double[]> tlNext = ThreadLocal.withInitial(() -> new double[8]);

        IntStream.range(0, W * H).parallel().forEach(idx -> {
            int i = idx % W;
            int j = idx / W;
            double[] y    = tlY.get();
            double[] next = tlNext.get();
            Integrator rk = tlInt.get();

            camera.initialState(i, j, y);
            RayOutcome outcome = traceRay(rk, y, next);
            image[j][i] = config.shader().shade(outcome, y);
        });

        return image;
    }

    /** Render and write a FITS file with a populated header. */
    public void renderToFits(Path outPath, String scene) throws IOException {
        float[][] image = render();
        FitsWriter.write(outPath, image, fitsMetadata(scene));
    }

    FitsWriter.Metadata fitsMetadata(String scene) {
        return new FitsWriter.Metadata(
                scene,
                metric.mass(),
                camera.rObs(),
                camera.inclination(),
                camera.fovRadians(),
                config.horizonCushion(),
                config.stepPolicyString(),
                config.maxSteps(),
                FitsWriter.gitHash());
    }

    /**
     * Integrate one ray from its initial state until it hits a termination
     * condition. Writes the final state back into {@code y}.
     */
    private RayOutcome traceRay(Integrator rk, double[] y, double[] next) {
        final double rCapture  = metric.horizonRadius() + config.horizonCushion();
        final double rEscape   = config.escapeRadius();
        final double rTrans    = config.rTransition();
        final double hNear     = config.hNear();
        final double hFar      = config.hFar();
        final int maxSteps     = config.maxSteps();

        for (int step = 0; step < maxSteps; step++) {
            double r = y[1];
            double h = (r > rTrans) ? hFar : hNear;
            rk.step(metric, y, h, next);
            // Copy back into y.
            System.arraycopy(next, 0, y, 0, 8);

            double rNow = y[1];
            if (Double.isNaN(rNow) || Double.isNaN(y[4])) {
                return RayOutcome.NAN;
            }
            if (rNow < rCapture) {
                return RayOutcome.CAPTURED;
            }
            if (rNow > rEscape) {
                return RayOutcome.ESCAPED;
            }
        }
        return RayOutcome.MAX_STEPS;
    }
}
