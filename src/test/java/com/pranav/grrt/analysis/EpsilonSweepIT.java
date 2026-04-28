package com.pranav.grrt.analysis;

import com.pranav.grrt.metric.KerrMetric;
import com.pranav.grrt.metric.Metric;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

/**
 * Phase 3C integration tests for {@link EpsilonSweep}, slow-tagged
 * (excluded from default {@code mvn test}; run via
 * {@code mvn verify -PrunSlow}).
 *
 * <p>Validation gates from {@code docs/phase-3c-plan.md} §7:
 * <ul>
 *   <li><b>Gate 3</b> Kerr/JP consistency at
 *       {@code (a = 0.9, ε₃ = 0, i = 17°, 512²)} — direct Kerr render
 *       vs JP(0.9, 0) frame produced by {@link EpsilonSweep},
 *       {@code |Δ(δ_r/⟨r⟩)| &lt; 0.001} (0.1 pp).</li>
 *   <li><b>Gate 4</b> External shadow-diameter cross-check
 *       (non-blocking sanity): Kerr {@code a = 0.9}, {@code i = 17°}
 *       ring radius lies in {@code [3, 7] M}.</li>
 *   <li><b>Gate 7</b> Resumability: interrupt + resume yields a CSV
 *       whose numeric columns match a single-shot run exactly.</li>
 * </ul>
 *
 * <p>Plus a few operational tests: header-mismatch fail-fast,
 * standalone resumability with skip semantics. All slow-tagged.
 */
@Tag("slow")
class EpsilonSweepIT {

    private static final int IT_RES_TINY = 64;
    private static final int IT_RES_SMALL = 256;
    private static final int IT_RES_FULL = 512;

    private static final double SPIN = 0.9;
    private static final double R_OUTER = 20.0;
    private static final double R_OBS = 1000.0;
    private static final double FOV_RAD = 30.0 / 1000.0;
    private static final double ATOL = 1.0e-10;
    private static final double RTOL = 1.0e-10;
    private static final double H_INITIAL = 1.0;
    private static final double HORIZON_CUSHION = 0.01;
    private static final int MAX_STEPS = 10_000_000;
    private static final int N_BINS = 180;
    private static final double PX_THRESHOLD = 2.0;

    // ------------------------------------------------------------------
    // Gate 3: Kerr/JP consistency at ε₃ = 0
    // ------------------------------------------------------------------

    @Test
    void gate3_kerrAndJpEpsilonZeroProduceConsistentAsymmetry(@TempDir Path tmp) throws IOException {
        EpsilonSweep.Config cfg = configFor(IT_RES_FULL, new double[] { 0.0 }, new double[] { 17.0 }, tmp);
        EpsilonSweep sweep = new EpsilonSweep();

        EpsilonSweep.SweepResult sr = sweep.runSweep(cfg);
        assertEquals(1, sr.rendered());
        assertEquals(0, sr.skipped());
        double jpAsym = readAsymmetryFromCsv(cfg.csvPath(), 0.0, 17.0, IT_RES_FULL);

        Metric kerr = new KerrMetric(SPIN);
        Path kerrFits = tmp.resolve("kerr_a0p9_i17p0_res512.fits");
        EpsilonSweep.FrameResult kerrFrame = sweep.renderAndAnalyze(kerr, 17.0, cfg, kerrFits);
        double kerrAsym = kerrFrame.stats().deltaRrms() / kerrFrame.stats().meanR();

        double diff = Math.abs(jpAsym - kerrAsym);
        assertTrue(diff < 0.001,
                "Kerr/JP consistency at ε₃=0: |Δ(δ_r/⟨r⟩)| < 0.001 expected, got "
                        + diff + " (JP=" + jpAsym + ", Kerr=" + kerrAsym + ")");
    }

    // ------------------------------------------------------------------
    // Gate 4: external shadow-diameter cross-check (non-blocking sanity)
    // ------------------------------------------------------------------

    @Test
    void gate4_kerrShadowRingRadiusInPhysicalRange(@TempDir Path tmp) throws IOException {
        EpsilonSweep.Config cfg = configFor(IT_RES_SMALL, new double[] { 0.0 }, new double[] { 17.0 }, tmp);
        EpsilonSweep sweep = new EpsilonSweep();
        Metric kerr = new KerrMetric(SPIN);
        Path fitsPath = tmp.resolve("kerr_a0p9_i17p0_res256.fits");

        EpsilonSweep.FrameResult fr = sweep.renderAndAnalyze(kerr, 17.0, cfg, fitsPath);

        double meanR = fr.stats().meanR();
        assertTrue(meanR >= 3.0 && meanR <= 7.0,
                "Kerr ring mean radius expected in [3, 7] M, got " + meanR);
        assertTrue(fr.stats().validBins() >= 0.5 * N_BINS,
                "expected at least half of bins to be valid, got "
                        + fr.stats().validBins() + " of " + N_BINS);
    }

    // ------------------------------------------------------------------
    // Gate 7: resumability — interrupt + resume vs single-shot
    // ------------------------------------------------------------------

    @Test
    void gate7_resumeMatchesSingleShotOnNumericColumns(@TempDir Path tmp) throws IOException {
        double[] fullGrid = { -0.1, 0.0, +0.1 };
        double[] firstFrame = { 0.0 };
        double[] incl = { 17.0 };
        EpsilonSweep sweep = new EpsilonSweep();

        Path tmpA = tmp.resolve("single-shot");
        Path tmpB = tmp.resolve("resumed");
        Files.createDirectories(tmpA);
        Files.createDirectories(tmpB);

        EpsilonSweep.Config cfgA = configFor(IT_RES_TINY, fullGrid, incl, tmpA);
        EpsilonSweep.SweepResult srA = sweep.runSweep(cfgA);
        assertEquals(3, srA.rendered());
        assertEquals(0, srA.skipped());

        EpsilonSweep.Config cfgB1 = configFor(IT_RES_TINY, firstFrame, incl, tmpB);
        EpsilonSweep.SweepResult srB1 = sweep.runSweep(cfgB1);
        assertEquals(1, srB1.rendered());
        assertEquals(0, srB1.skipped());

        EpsilonSweep.Config cfgB2 = configFor(IT_RES_TINY, fullGrid, incl, tmpB);
        EpsilonSweep.SweepResult srB2 = sweep.runSweep(cfgB2);
        assertEquals(2, srB2.rendered(), "resume should render 2 new frames");
        assertEquals(1, srB2.skipped(), "resume should skip the previously-rendered frame");

        List<String> rowsA = dataRowsSortedByEpsilon(cfgA.csvPath());
        List<String> rowsB = dataRowsSortedByEpsilon(cfgB2.csvPath());
        assertEquals(3, rowsA.size());
        assertEquals(3, rowsB.size());
        for (int i = 0; i < 3; i++) {
            String[] a = rowsA.get(i).split(",", -1);
            String[] b = rowsB.get(i).split(",", -1);
            // Columns 0..8 are numeric (epsilon_3 .. multipeak_bins); 9 is wallclock,
            // 10 is git_sha, 11 is timestamp — skip 9-11 (legitimately differ).
            for (int c = 0; c <= 8; c++) {
                if (c == 2 || c == 8) {
                    assertEquals(Long.parseLong(a[c]), Long.parseLong(b[c]),
                            "column " + c + " mismatch in row " + i);
                } else {
                    double va = Double.parseDouble(a[c]);
                    double vb = Double.parseDouble(b[c]);
                    assertEquals(va, vb, Math.max(Math.abs(va), 1.0) * 1.0e-12,
                            "numeric column " + c + " mismatch in row " + i
                                    + ": " + va + " vs " + vb);
                }
            }
        }
    }

    // ------------------------------------------------------------------
    // CSV header validation + helper smoke
    // ------------------------------------------------------------------

    @Test
    void csvHeaderMismatchOnResumeFailsFast(@TempDir Path tmp) throws IOException {
        Path csvPath = tmp.resolve("sweep.csv");
        Files.writeString(csvPath, "wrong,header,here\n", StandardCharsets.UTF_8);

        EpsilonSweep.Config cfg = configFor(IT_RES_TINY, new double[] { 0.0 }, new double[] { 17.0 }, tmp);
        EpsilonSweep.Config bogus = new EpsilonSweep.Config(
                cfg.spin(), cfg.epsilonGrid(), cfg.inclinationDegrees(), cfg.resolution(),
                cfg.rOuter(), csvPath, cfg.fitsDir(),
                cfg.rObs(), cfg.fovRadians(), cfg.atol(), cfg.rtol(), cfg.hInitial(),
                cfg.horizonCushion(), cfg.maxSteps(), cfg.numBins(), cfg.pxThreshold());

        EpsilonSweep sweep = new EpsilonSweep();
        IOException thrown = assertThrows(IOException.class, () -> sweep.runSweep(bogus));
        assertTrue(thrown.getMessage().contains("CSV header mismatch"),
                "expected header-mismatch IOException, got: " + thrown.getMessage());
    }

    @Test
    void emptyCsvAtPathIsTreatedAsFresh(@TempDir Path tmp) throws IOException {
        Path csvPath = tmp.resolve("sweep.csv");
        Files.writeString(csvPath, "", StandardCharsets.UTF_8);

        Set<String> keys = EpsilonSweep.readCompletedKeys(csvPath);
        assertTrue(keys.isEmpty());
    }

    @Test
    void formatKeyIsLocaleIndependentAndExactMatch() {
        assertEquals("+0.0000|17.0|512", EpsilonSweep.formatKey(0.0, 17.0, 512));
        assertEquals("-2.5000|17.0|512", EpsilonSweep.formatKey(-2.5, 17.0, 512));
        assertEquals("+0.1200|60.0|512", EpsilonSweep.formatKey(+0.12, 60.0, 512));
        // Round-trip: parse the key, reformat, expect identity.
        double eps = -0.05;
        double inc = 17.0;
        int res = 64;
        String once = EpsilonSweep.formatKey(eps, inc, res);
        String[] parts = once.split("\\|");
        double parsedEps = Double.parseDouble(parts[0]);
        double parsedInc = Double.parseDouble(parts[1]);
        int parsedRes = Integer.parseInt(parts[2]);
        assertEquals(once, EpsilonSweep.formatKey(parsedEps, parsedInc, parsedRes));
    }

    // ------------------------------------------------------------------
    // Standalone resumability with explicit skip-count assertions
    // ------------------------------------------------------------------

    @Test
    void resumeSkipsAlreadyCompletedKeysWithoutRerendering(@TempDir Path tmp) throws IOException {
        EpsilonSweep sweep = new EpsilonSweep();
        EpsilonSweep.Config cfg = configFor(IT_RES_TINY, new double[] { 0.0 }, new double[] { 17.0 }, tmp);

        EpsilonSweep.SweepResult first = sweep.runSweep(cfg);
        assertEquals(1, first.rendered());
        assertEquals(0, first.skipped());
        long firstFitsMtime = Files.getLastModifiedTime(
                cfg.fitsDir().resolve("sweep_a0.9_eps+0.0000_i17.0_res64.fits")).toMillis();

        try {
            Thread.sleep(50);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            fail("interrupted");
        }

        EpsilonSweep.SweepResult second = sweep.runSweep(cfg);
        assertEquals(0, second.rendered());
        assertEquals(1, second.skipped());

        long secondFitsMtime = Files.getLastModifiedTime(
                cfg.fitsDir().resolve("sweep_a0.9_eps+0.0000_i17.0_res64.fits")).toMillis();
        assertEquals(firstFitsMtime, secondFitsMtime,
                "skipped frame's FITS file was unexpectedly re-written");

        List<String> lines = Files.readAllLines(cfg.csvPath(), StandardCharsets.UTF_8);
        assertEquals(2, lines.size());
        assertEquals(EpsilonSweep.CSV_HEADER, lines.get(0));
    }

    // ------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------

    private static EpsilonSweep.Config configFor(int resolution, double[] epsGrid,
                                                 double[] inclDeg, Path workDir) {
        Path csvPath = workDir.resolve("sweep.csv");
        Path fitsDir = workDir;
        return new EpsilonSweep.Config(
                SPIN, epsGrid, inclDeg, resolution, R_OUTER, csvPath, fitsDir,
                R_OBS, FOV_RAD, ATOL, RTOL, H_INITIAL, HORIZON_CUSHION,
                MAX_STEPS, N_BINS, PX_THRESHOLD);
    }

    private static double readAsymmetryFromCsv(Path csvPath, double eps, double inclDeg, int resolution)
            throws IOException {
        String key = EpsilonSweep.formatKey(eps, inclDeg, resolution);
        List<String> lines = Files.readAllLines(csvPath, StandardCharsets.UTF_8);
        for (int i = 1; i < lines.size(); i++) {
            String[] parts = lines.get(i).split(",", -1);
            if (parts.length < 5) {
                continue;
            }
            double rowEps = Double.parseDouble(parts[0]);
            double rowInc = Double.parseDouble(parts[1]);
            int rowRes = Integer.parseInt(parts[2]);
            if (EpsilonSweep.formatKey(rowEps, rowInc, rowRes).equals(key)) {
                double meanR = Double.parseDouble(parts[3]);
                double rms = Double.parseDouble(parts[4]);
                assertNotEquals(0.0, meanR, "mean_r should be non-zero");
                return rms / meanR;
            }
        }
        fail("CSV does not contain a row for key " + key);
        return Double.NaN;
    }

    private static List<String> dataRowsSortedByEpsilon(Path csvPath) throws IOException {
        List<String> lines = Files.readAllLines(csvPath, StandardCharsets.UTF_8);
        assertNotNull(lines);
        if (lines.isEmpty()) {
            return List.of();
        }
        return lines.subList(1, lines.size()).stream()
                .filter(s -> !s.isBlank())
                .sorted((a, b) -> {
                    String aEps = a.substring(0, a.indexOf(','));
                    String bEps = b.substring(0, b.indexOf(','));
                    return Double.compare(Double.parseDouble(aEps), Double.parseDouble(bEps));
                })
                .toList();
    }
}
