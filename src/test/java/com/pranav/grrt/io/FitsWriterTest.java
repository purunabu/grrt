package com.pranav.grrt.io;

import nom.tam.fits.BasicHDU;
import nom.tam.fits.Fits;
import nom.tam.fits.Header;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.jupiter.api.Assertions.*;

class FitsWriterTest {

    // ---------------------------------------------------------------
    // d. FITS round-trip: write an image + header, read it back,
    //    verify bit-identical pixel array and all expected header keys.
    // ---------------------------------------------------------------

    @Test
    void roundTripPreservesImageAndHeader(@TempDir Path tmp) throws Exception {
        int W = 17, H = 13;  // deliberately prime/non-power-of-two
        float[][] image = new float[H][W];
        for (int j = 0; j < H; j++) {
            for (int i = 0; i < W; i++) {
                image[j][i] = (float) (j * 1000 + i);
            }
        }

        FitsWriter.Metadata meta = new FitsWriter.Metadata(
                "test-scene",
                1.0, 1000.0, Math.PI / 2, 0.02,
                0.01, "h=0.01 r<10M else h=0.5", 1_000_000,
                "0123456789abcdef0123456789abcdef01234567");

        Path out = tmp.resolve("roundtrip.fits");
        FitsWriter.write(out, image, meta);
        assertTrue(Files.exists(out), "output file not written");
        assertTrue(Files.size(out) > 0, "output file empty");

        try (Fits fits = new Fits(out.toFile())) {
            BasicHDU<?>[] hdus = fits.read();
            assertEquals(1, hdus.length, "expected 1 HDU");

            // Pixel round-trip.
            Object kernel = hdus[0].getData().getKernel();
            assertInstanceOf(float[][].class, kernel);
            float[][] readBack = (float[][]) kernel;
            assertEquals(H, readBack.length);
            assertEquals(W, readBack[0].length);
            for (int j = 0; j < H; j++) {
                assertArrayEquals(image[j], readBack[j], 0.0f, "row " + j);
            }

            // Header spot-checks.
            Header h = hdus[0].getHeader();
            assertEquals(1.0,        h.getDoubleValue("BHMASS"),  0.0);
            assertEquals(1000.0,     h.getDoubleValue("ROBS"),    0.0);
            assertEquals(Math.PI/2,  h.getDoubleValue("INCL"),    1e-15);
            assertEquals(0.02,       h.getDoubleValue("FOV"),     0.0);
            assertEquals(W,          h.getIntValue("NPIXX"));
            assertEquals(H,          h.getIntValue("NPIXY"));
            assertEquals(0.01,       h.getDoubleValue("CUSHION"), 0.0);
            assertEquals(1_000_000,  h.getIntValue("MAXSTEPS"));
            assertEquals("test-scene", h.getStringValue("SCENE"));
            assertEquals("grrt",     h.getStringValue("CREATOR"));
            assertEquals("0123456789abcdef0123456789abcdef01234567",
                         h.getStringValue("GITHASH"));

            // WCS sanity: CDELT2 must be negative so FITS viewers flip
            // top-to-bottom to match our j=0-is-top convention.
            double cdelt1 = h.getDoubleValue("CDELT1");
            double cdelt2 = h.getDoubleValue("CDELT2");
            assertEquals(0.02 / W, cdelt1, 1e-15, "CDELT1");
            assertEquals(-0.02 / W, cdelt2, 1e-15, "CDELT2 (must be negative)");
            assertEquals(0.5 * W + 0.5, h.getDoubleValue("CRPIX1"), 1e-15);
            assertEquals(0.5 * H + 0.5, h.getDoubleValue("CRPIX2"), 1e-15);
        }
    }

    // ---------------------------------------------------------------
    // gitHash() best-effort: in this repo it should return a 40-char
    // SHA1 hex. If run outside a git repo it should return "". Either
    // way it must not throw.
    // ---------------------------------------------------------------

    @Test
    void gitHashIsBestEffortAndNeverThrows() {
        String hash = FitsWriter.gitHash();
        assertNotNull(hash);
        if (!hash.isEmpty()) {
            assertTrue(hash.matches("[0-9a-f]{40}"),
                    "expected 40-char lowercase SHA1, got: " + hash);
        }
    }
}
