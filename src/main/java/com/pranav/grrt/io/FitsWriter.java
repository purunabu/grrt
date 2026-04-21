package com.pranav.grrt.io;

import nom.tam.fits.BasicHDU;
import nom.tam.fits.Fits;
import nom.tam.fits.FitsException;
import nom.tam.fits.Header;
import nom.tam.util.FitsFile;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Path;
import java.time.Instant;
import java.util.concurrent.TimeUnit;

/**
 * Writes grrt renderer output to a FITS file with a populated header.
 *
 * <h2>Header keys</h2>
 * <ul>
 *   <li>{@code BHMASS, ROBS, INCL, FOV}           - physical scene parameters</li>
 *   <li>{@code NPIXX, NPIXY}                       - image dimensions</li>
 *   <li>{@code CUSHION, HPOLICY, MAXSTEPS}         - renderer numerics</li>
 *   <li>{@code GITHASH, CREATOR, DATE, SCENE}      - provenance</li>
 *   <li>{@code CRPIX1/2, CDELT1/2, CRVAL1/2,
 *       CTYPE1/2, CUNIT1/2}                        - WCS (linear tangent plane)</li>
 * </ul>
 *
 * <p><b>On CDELT2 being negative.</b> FITS 2D images are written with
 * NAXIS1 (columns) fastest and NAXIS2 (rows) next; viewers like DS9 display
 * pixel (1,1) at the bottom-left by default. Our image is row-major with
 * {@code image[0]} as the top row (screen-up convention). Setting
 * {@code CDELT2 < 0} tells the WCS that increasing pixel-j corresponds to
 * decreasing sky-β, so DS9/astropy/CARTA correctly render our top row as
 * the top of the display.
 */
public final class FitsWriter {

    private FitsWriter() {}

    /**
     * Scene + renderer metadata stamped into the FITS header.
     *
     * @param scene           short label, e.g. "schwarzschild-shadow"
     * @param bhMass          BH mass M (geometrized units)
     * @param rObs            observer radius / M
     * @param inclination     observer polar angle / rad
     * @param fovRadians      full horizontal FOV / rad
     * @param horizonCushion  ε_h / M
     * @param stepPolicy      human-readable step-size rule
     * @param maxSteps        integrator safety cap per ray
     * @param gitHash         commit hash, or empty string
     */
    public record Metadata(
            String scene,
            double bhMass,
            double rObs,
            double inclination,
            double fovRadians,
            double horizonCushion,
            String stepPolicy,
            int maxSteps,
            String gitHash
    ) {}

    /**
     * Write {@code image[H][W]} (row 0 = top of image) to {@code outPath}
     * in FITS format, stamping {@code meta} into the header.
     *
     * @throws IOException if the file cannot be written
     */
    public static void write(Path outPath, float[][] image, Metadata meta) throws IOException {
        int H = image.length;
        int W = image[0].length;

        try {
            Fits fits = new Fits();
            BasicHDU<?> hdu = Fits.makeHDU(image);
            fits.addHDU(hdu);

            Header h = hdu.getHeader();
            // Scene parameters
            h.addValue("BHMASS",   meta.bhMass(),         "BH mass M in geometrized units");
            h.addValue("ROBS",     meta.rObs(),           "observer r-coordinate / M");
            h.addValue("INCL",     meta.inclination(),    "inclination theta_obs / rad");
            h.addValue("FOV",      meta.fovRadians(),     "full horizontal FOV / rad");
            h.addValue("NPIXX",    W,                     "image width");
            h.addValue("NPIXY",    H,                     "image height");

            // Renderer numerics
            h.addValue("CUSHION",  meta.horizonCushion(), "horizon cushion / M");
            h.addValue("HPOLICY",  truncate(meta.stepPolicy(), 68), "step-size policy");
            h.addValue("MAXSTEPS", meta.maxSteps(),       "RK4 step cap per ray");

            // Provenance
            h.addValue("SCENE",    truncate(meta.scene(), 68), "scene label");
            h.addValue("CREATOR",  "grrt",                "GR ray tracer");
            h.addValue("DATE",     Instant.now().toString(), "creation time (UTC)");
            h.addValue("GITHASH",  truncate(meta.gitHash(), 68),
                       "commit hash (empty if not a git repo)");

            // WCS: linear tangent-plane coordinates in radians.
            // Reference pixel = geometric center of the image (1-indexed).
            double cdelt = meta.fovRadians() / W;   // rad/pixel (square pixels)
            h.addValue("CTYPE1",   "LINEAR",              "alpha (screen-right)");
            h.addValue("CTYPE2",   "LINEAR",              "beta (screen-up)");
            h.addValue("CUNIT1",   "rad",                 "units of CDELT1");
            h.addValue("CUNIT2",   "rad",                 "units of CDELT2");
            h.addValue("CRPIX1",   0.5 * W + 0.5,         "reference pixel, 1-indexed");
            h.addValue("CRPIX2",   0.5 * H + 0.5,         "reference pixel, 1-indexed");
            h.addValue("CRVAL1",   0.0,                   "alpha at reference pixel");
            h.addValue("CRVAL2",   0.0,                   "beta at reference pixel");
            h.addValue("CDELT1",   cdelt,                 "rad/pixel along axis 1");
            h.addValue("CDELT2",   -cdelt,                "rad/pixel axis 2 (neg: j=0 is top)");

            try (FitsFile ff = new FitsFile(outPath.toString(), "rw")) {
                fits.write(ff);
            }
        } catch (FitsException e) {
            throw new IOException("failed to write FITS to " + outPath, e);
        }
    }

    /**
     * Best-effort git commit hash via {@code git rev-parse HEAD}. Returns the
     * empty string on any failure (not a repo, git not installed, timeout).
     * Never throws.
     */
    public static String gitHash() {
        try {
            Process p = new ProcessBuilder("git", "rev-parse", "HEAD")
                    .redirectErrorStream(true)
                    .start();
            String line;
            try (BufferedReader r = new BufferedReader(
                    new InputStreamReader(p.getInputStream(), StandardCharsets.UTF_8))) {
                line = r.readLine();
            }
            if (!p.waitFor(2, TimeUnit.SECONDS)) {
                p.destroyForcibly();
                return "";
            }
            if (p.exitValue() != 0 || line == null) {
                return "";
            }
            return line.trim();
        } catch (IOException | InterruptedException e) {
            if (e instanceof InterruptedException) {
                Thread.currentThread().interrupt();
            }
            return "";
        }
    }

    private static String truncate(String s, int max) {
        if (s == null) return "";
        return s.length() <= max ? s : s.substring(0, max);
    }

}
