package edu.ucsc.neurobiology.vision.stimulus;

import java.io.*;
import static java.lang.Math.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.util.*;


/*
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class RevGratBarFrameGenerator
    extends FrameGenerator {

    final int nReversals;
    final double tPeriod;
    final double sPeriod;
    final double contrast;
    final int nFrames;
    final double sPhase;
    final double barLocation;

    int frameIndex = -1;

    public String toString() {
        return
            "REVERSING GRATING BAR: " +
            "location " + barLocation + ", " +
            "sPeriod " + sPeriod + ", " +
            "sPhase " + sPhase + ", " +
            "tPeriod " + tPeriod + ", " +
            "contrast " + contrast + ", " +
            "nFrames " + nFrames + ", "
            ;
    }


    public RevGratBarFrameGenerator(int barLocation, double tPeriod, double sPeriod,
                                    double sPhase,
                                    int nReversalsPerCondition, double contrast) {

        super(640, 1, RandomNumberGenerator.JAVA_RANDOM, 0);

        this.tPeriod = tPeriod;
        this.sPeriod = sPeriod;
        this.nReversals = nReversalsPerCondition;
        this.contrast = contrast;
        this.barLocation = barLocation;
        this.sPhase = sPhase;

        nFrames = (int) (nReversalsPerCondition * tPeriod);
    }


    public int getNumberOfFrames() {
        return nFrames;
    }


    public void nextFrame(Object _frame) {
        if (_frame instanceof int[]) {
            throw new IllegalStateException(
                "Frames are not created as int[] with raw movies.");
        }
        float[] frame = (float[]) _frame;

        frameIndex++;

        double x1 = barLocation - sPeriod / 2;
        double x2 = barLocation + sPeriod / 2;
        double tValue = contrast * sin(2 * PI * frameIndex / tPeriod);

        for (int x = 0; x < width; x++) {
            double v;

            if (x >= x1 && x < x2) {
                v = 0.5 + tValue * sin(2 * PI * (x - x1 + sPhase) / sPeriod);
            } else {
                v = 0.5f;
            }

            frame[3 * x + 0] = frame[3 * x + 1] = frame[3 * x + 2] = (float) v;
        }
    }


    public Object createFrameBuffer() {
        return new int[width * height * 3];
    }


    public static void main(String args[]) throws IOException {
        String movieFileName = "f:\\movies\\barCRG.rawMovie";

        int width = 640;
        int height = 1;
        final double tPeriod = 15;
        final double stepSize = 20;
        final double contrast = 0.48;
        final int hPadding = 100;
        final int nReversalsPerCondition = 11 * 8;
        int[] sPeriod = {24, 48, 96};

        ArrayList<RevGratBarFrameGenerator> list = new ArrayList();
        for (int barLocation = hPadding; barLocation <= 640 - hPadding;
                               barLocation += stepSize) {
            for (int i = 0; i < sPeriod.length; i++) {
                for (int j = 0; j < 5; j++) {
                    list.add(new RevGratBarFrameGenerator(
                        barLocation, tPeriod, sPeriod[i], j * sPeriod[i] / 8.0,
                        nReversalsPerCondition, contrast));
                }
            }
        }

        GrayFrameGenerator gray = new GrayFrameGenerator(15);

        FrameGenerator[] array = new FrameGenerator[2 * list.size()];
        Random r = new Random(11111);
        for (int i = 0; i < array.length / 2; i++) {
            int n = (int) (r.nextDouble() * list.size());
            array[2 * i + 0] = list.remove(n);
            array[2 * i + 1] = gray;
        }

        PrintWriter pw = new PrintWriter("f:\\movies\\barCRG.txt");
        for (int i = 0; i < array.length; i++) {
            pw.println(array[i]);
        }
        pw.close();

        JoiningFrameGenerator g = new JoiningFrameGenerator(array);
        int nFrames = g.getNumberOfFrames();
        System.err.println("Movie Length: " + nFrames / 120. + "s.");

        RawMovieFile file = new RawMovieFile(
            movieFileName, width, height, nFrames, "Reversing Gratings Bar", 0
            , "tPeriod" + AsciiHeaderFile.fieldSeparator + tPeriod
            , "stepSize" + AsciiHeaderFile.fieldSeparator + stepSize
            , "contrast" + AsciiHeaderFile.fieldSeparator + contrast
            , "hPadding" + AsciiHeaderFile.fieldSeparator + hPadding
            ,
            "nReversalsPerCondition" + AsciiHeaderFile.fieldSeparator +
            nReversalsPerCondition
            , "sPeriod" + AsciiHeaderFile.fieldSeparator + StringUtil.toString(sPeriod)
                            );
        
        float[] colorFrame = new float[width * height * 3];

        for (int i = 0; i < nFrames; i++) {
            g.nextFrame(colorFrame);
            file.writeFrame( (float[]) colorFrame);
        }
        
        file.close();


    }

}
