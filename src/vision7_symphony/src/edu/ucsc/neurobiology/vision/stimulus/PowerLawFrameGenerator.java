package edu.ucsc.neurobiology.vision.stimulus;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.math.*;


/*
 * This class defines how to create each frame of a pink movie.
 *
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class PowerLawFrameGenerator
    extends FrameGenerator {


    private float[] sigma; //standard deviation of final stimulus.

    private float normalizer; //used to normalize the frames to the correct standard deviation
    private float[] currentColors; //current frame
    private int colorType;
    private MovieType noiseType; //RGB vs GRAY and BINARY_MOVIE vs GAUSSIAN_MOVIE
    public static final int RGB = 0, GRAY = 2;
//    FFT fft;

    int nFrames;

    private int kXSize, kYSize, omegaSize;
    private double spatialPower, temporalPower;
    private String tempFileName;
    FileInputStream[] fis;
    BufferedInputStream[] fileReader;

    private int memoryAvailable = 1 << 30;

    RandomMersenne mersenne;

    int pixelsAtATime;
    int blocks;

    public PowerLawFrameGenerator(int width, int height, int spatialPeriodCutoff,
                                  int temporalPeriodCutoff, double spatialPower,
                                  double temporalPower, float[] sigma, int startSeed,
                                  int colorType, MovieType noiseType, int nFrames,
                                  String tempMovieFileName) throws IOException {

        super(width, height, RandomNumberGenerator.MAC_RANDOM, 0);
        this.sigma = sigma;
        this.colorType = colorType;
        this.noiseType = noiseType;
        this.spatialPower = spatialPower;
        this.temporalPower = temporalPower;
        this.nFrames = nFrames;
        this.tempFileName = tempMovieFileName;

        mersenne = new RandomMersenne(startSeed);
        double powerRe[][][], powerIm[][][];

        kXSize = 1;
        while (kXSize < width) {
            kXSize <<= 1;
        }

        kYSize = 1;
        while (kYSize < height) {
            kYSize <<= 1;
        }

        omegaSize = 1;
        while (omegaSize < nFrames) {
            omegaSize <<= 1;
        }

        //used to put x and y on the same length scale.
        double aspectRatio = (double) kXSize / kYSize;

        double maxSpatialNormalization = Math.pow(2.0 * kXSize / spatialPeriodCutoff *
                                                  kYSize * aspectRatio /
                                                  spatialPeriodCutoff,
                                                  .5 * spatialPower);
        double maxTemporalNormalization = Math.pow( (double) omegaSize /
            temporalPeriodCutoff,
            temporalPower);
        currentColors = new float[width * height * 3];
        pixelsAtATime = Math.min(memoryAvailable / 16 / omegaSize /
                                 (colorType == RGB ? 3 : 1),
                                 kXSize * kYSize);
        while (kXSize * kYSize % pixelsAtATime != 0) {
            if (pixelsAtATime != 0) {
                pixelsAtATime--;
            }
        }
        blocks = kXSize * kYSize / pixelsAtATime;

        if (pixelsAtATime == 0) {
            throw new IllegalStateException(
                "Insufficient memory specified to generate a movie of this length");
        }

        powerRe = new double[ (colorType == RGB ? 3 : 1)][pixelsAtATime][omegaSize];
        powerIm = new double[ (colorType == RGB ? 3 : 1)][pixelsAtATime][omegaSize];
        byte[] byteArray = new byte[2 * 8 * (colorType == RGB ? 3 : 1) * pixelsAtATime];

        int refPixel = 0;
        double totalVariance = 0.0;

        for (int block = 0; block < blocks; block++) {
            FileOutputStream fos = new FileOutputStream(tempFileName +
                Integer.toString(block));
            BufferedOutputStream fileWriter = new BufferedOutputStream(fos);

            for (int pix = 0; pix < pixelsAtATime; pix++) {
                int currentPixel = refPixel + pix;

                int x = currentPixel % kXSize;
                int y = currentPixel / kXSize;
                for (int c = 0; c < (colorType == RGB ? 3 : 1); c++) {

                    //All waves going backwards in time are zero.
                    for (int t = 0; t < omegaSize; t++) {
                        if (t > omegaSize / 2 - 1) {

                            double phase = mersenne.nextPhase();
                            double temporalNormalization = Math.pow(t <
                                omegaSize / 2 + 1 ? t :
                                (omegaSize - t),
                                temporalPower);
                            double spatialNormalization = Math.pow( (x <
                                kXSize / 2 + 1 ?
                                x * x :
                                (kXSize - x) * (kXSize - x)) +
                                (y < kYSize / 2 + 1 ?
                                 y * y * aspectRatio * aspectRatio :
                                 (kYSize - y) * (kYSize - y) * aspectRatio *
                                 aspectRatio),
                                .5 * spatialPower);

                            double normalization = Math.min(temporalNormalization,
                                maxTemporalNormalization) *
                                Math.min(spatialNormalization,
                                         maxSpatialNormalization);

                            powerRe[c][pix][t] = normalization * Math.cos(phase);
                            powerIm[c][pix][t] = normalization * Math.sin(phase);
                            if (c == 0) {
                                totalVariance += .5 * normalization *
                                    normalization; // sigma from one sine wave
                            }
                        } else {
                            powerRe[c][pix][t] = 0;
                            powerIm[c][pix][t] = 0;
                        }

                    }

                    FFT.fft(powerRe[c][pix], powerIm[c][pix], 1);
                }
            }
            System.out.println("refPixel: " + refPixel);
            System.out.println("pixelsAtATime: " + pixelsAtATime);

            Date start = new Date();
            for (int t = 0; t < nFrames; t++) {

                //This method of writing bytes is faster than writing the doubles directly
                //because the bytes can be written as an array.
                int location = 0;
                for (int pix = 0; pix < pixelsAtATime; pix++) {
                    for (int c = 0; c < (colorType == RGB ? 3 : 1); c++) {
                        long re = Double.doubleToLongBits(powerRe[c][pix][t]);
                        long im = Double.doubleToLongBits(powerIm[c][pix][t]);

                        byteArray[location++] = (byte) (re >> 56);
                        byteArray[location++] = (byte) (re >> 48);
                        byteArray[location++] = (byte) (re >> 40);
                        byteArray[location++] = (byte) (re >> 32);
                        byteArray[location++] = (byte) (re >> 24);
                        byteArray[location++] = (byte) (re >> 16);
                        byteArray[location++] = (byte) (re >> 8);
                        byteArray[location++] = (byte) re;
                        byteArray[location++] = (byte) (im >> 56);
                        byteArray[location++] = (byte) (im >> 48);
                        byteArray[location++] = (byte) (im >> 40);
                        byteArray[location++] = (byte) (im >> 32);
                        byteArray[location++] = (byte) (im >> 24);
                        byteArray[location++] = (byte) (im >> 16);
                        byteArray[location++] = (byte) (im >> 8);
                        byteArray[location++] = (byte) im;

                    }
                }
                fileWriter.write(byteArray);

            }
            fileWriter.flush();
            fileWriter.close();
            fos.flush();
            fos.close();

            Date end = new Date();
            System.out.println( (end.getTime() - start.getTime()) / 1000.0);

            refPixel += pixelsAtATime;

        }

        fis = new FileInputStream[blocks];
        fileReader = new BufferedInputStream[blocks];

        for (int i = 0; i < blocks; i++) {
            fis[i] = new FileInputStream(tempFileName + Integer.toString(i));
            fileReader[i] = new BufferedInputStream(fis[i]);
        }

        normalizer = (float) Math.sqrt(totalVariance) / sigma[0];
//normalizer=1f;
    }


    public Object createFrameBuffer() {
        return new int[width * height * 3];
    }


    public void nextFrame(Object frame) throws IOException {
        if (frame instanceof float[]) {
            nextFrame( (float[]) frame);
        } else if (frame instanceof int[]) {
            throw new IllegalStateException(
                "Frames are not created as int[] with power law movies.");
        }
    }


    private final void nextFrame(float[] color) throws IOException {
        currentColors = createFrame(currentColors);
        switch (colorType) {
            case RGB:

                //put colors in correct form for display.
                for (int i = 0; i < currentColors.length; i++) {
//                        color[i] = 2*(currentColors[i]/normalizer)*(currentColors[i]/normalizer )+.5f;
                    color[i] = (currentColors[i] / normalizer) + .5f;
                    if (noiseType == MovieType.GAUSSIAN_MOVIE) {
                        if (color[i] < 1.0 / 257.0) {
                            color[i] = (float) (1.0 / 257.0);
                        } else if (color[i] > 256.0 / 257.0) {
                            color[i] = (float) (256.0 / 257.0);
                        }
                    } else {
                        if (color[i] < .5) {
                            color[i] = (float) .5 - sigma[i % 3];
                        } else {
                            color[i] = (float) .5 + sigma[i % 3];
                        }
                    }

                }
                break;

            case GRAY:
                for (int i = 0; i < currentColors.length; i += 3) {
                    color[i] = (currentColors[i] / normalizer) + .5f;
//                           color[i] = 2*(currentColors[i]/normalizer)*(currentColors[i]/normalizer )+.5f;
                    if (noiseType == MovieType.GAUSSIAN_MOVIE) {
                        if (color[i] < 1.0 / 257.0) {
                            color[i] = (float) (1.0 / 257.0);
                        } else if (color[i] > 256.0 / 257.0) {
                            color[i] = (float) (256.0 / 257.0);

                        }
                    } else {

                        if (color[i] < .5) {
                            color[i] = (float) .5 - sigma[i % 3];
                        } else {
                            color[i] = (float) .5 + sigma[i % 3];

                        }
                    }

                    color[i + 1] = color[i];
                    color[i + 2] = color[i];

                }

                break;
        }

    }


    private float[] createFrame(float[] colors) throws IOException {
        Arrays.fill(colors, 0);

        double aveNorm = 0.0;
        byte[] bytes = new byte[2 * 8];
        double powerSpectrum[][][][] = new double[2][ (colorType == RGB ? 3 : 1)][kXSize][
                                       kYSize];

        for (int p = 0; p < kXSize * kYSize; p++) {
            for (int c = 0; c < (colorType == RGB ? 3 : 1); c++) {
                int x = p % kXSize;
                int y = p / kXSize;

                fileReader[p / pixelsAtATime].read(bytes);

                //reading as bytes is faster because BufferedInputSream allows block,
                // buffered reading.
                //block reading of doubles is not directly supported.

                // ? : is necessary because casting a byte as a long moves the sign
                // bit.
                long b0 = (bytes[0] < 0) ? bytes[0] + 256 : bytes[0];
                long b1 = (bytes[1] < 0) ? bytes[1] + 256 : bytes[1];
                long b2 = (bytes[2] < 0) ? bytes[2] + 256 : bytes[2];
                long b3 = (bytes[3] < 0) ? bytes[3] + 256 : bytes[3];
                long b4 = (bytes[4] < 0) ? bytes[4] + 256 : bytes[4];
                long b5 = (bytes[5] < 0) ? bytes[5] + 256 : bytes[5];
                long b6 = (bytes[6] < 0) ? bytes[6] + 256 : bytes[6];
                long b7 = (bytes[7] < 0) ? bytes[7] + 256 : bytes[7];

                long re = b7 +
                          (b6 << 8) +
                          (b5 << 16) +
                          (b4 << 24) +
                          (b3 << 32) +
                          (b2 << 40) +
                          (b1 << 48) +
                          (b0 << 56);

                b0 = (bytes[8] < 0) ? bytes[8] + 256 : bytes[8];
                b1 = (bytes[9] < 0) ? bytes[9] + 256 :
                     bytes[9];
                b2 = (bytes[10] < 0) ? bytes[10] + 256 :
                     bytes[10];
                b3 = (bytes[11] < 0) ?
                     bytes[11] + 256 : bytes[11];
                b4 = (bytes[12] < 0) ?
                     bytes[12] + 256 : bytes[12];
                b5 = (bytes[13] < 0) ? bytes[13] + 256 :
                     bytes[13];
                b6 = (bytes[14] < 0) ? bytes[14] + 256 :
                     bytes[14];
                b7 = (bytes[15] < 0) ? bytes[15] + 256 :
                     bytes[15];

                long im = b7 +
                          (b6 << 8) +
                          (b5 << 16) +
                          (b4 << 24) +
                          (b3 << 32) +
                          (b2 << 40) +
                          (b1 << 48) +
                          (b0 << 56);

                powerSpectrum[0][c][x][y] = Double.longBitsToDouble(re);
                powerSpectrum[1][c][x][y] = Double.longBitsToDouble(im);

            }
        }

        for (int c = 0; c < (colorType == RGB ? 3 : 1); c++) {
            FFT.fft2(powerSpectrum[0][c], powerSpectrum[1][c], 1);
//             algebra.FFTXY(powerSpectrum[c], 1);
        }

        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                for (int c = 0; c < 3; c++) {

                    colors[3 * (y * width + x) +
                        c] = (float) powerSpectrum[0][ (colorType == RGB ? c : 0)][x][y];
                }
            }
        }

        return colors;
    }


    public void deleteTempFile() throws IOException {
        for (int i = 0; i < blocks; i++) {
            fileReader[i].close();
            fis[i].close();
            File file = new File(tempFileName + Integer.toString(i));
            file.delete();
        }
    }

}
