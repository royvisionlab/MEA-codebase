package edu.ucsc.neurobiology.vision.testing;

import java.io.*;

import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.stimulus.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author nobody, anyone can change
 */
public class TestStimuli {
    //This main function is used to check a movie to make sure that it has the proper
    //statistics
    public static void testMovieStatistics() {
        try {

            RawMovie rawMovie = new RawMovie("G://movies2005-03_8//99999RGBlong.rawMovie",
                                             "G://data//2005-05-02-0//data002//data002.globals");
            BufferedMovie movie = new BufferedMovie(rawMovie);
            int staOffset = -1;
            STACalculator c = new STACalculator(1, staOffset, 1, .05, movie);

            double samplesPerFrame = 20 * rawMovie.getRefreshTime();

            int percent = -1, oldPercent = -1;
            int nPoints = 100000; //100000;
            int time[] = new int[nPoints];
            double mean[] = new double[nPoints];
            double variance[] = new double[nPoints];
            for (int i = 0; i < time.length; i++) {
                time[i] = (int) Math.round(i * samplesPerFrame * 1);
            }

            for (int i = 0; i < time.length; i++) {
                c.addSpike(time[i]);

                percent = (int) (100.0 * i / time.length);
                if (percent != oldPercent) {
                    System.out.println(percent + " %");
                    oldPercent = percent;
                }

                c.finish();

                STA sta = c.getSTA();
//                mean[i] = sta.getFrame(0).getPixel(0, 0, 1);
                sta = c.getSTV();
                if (i != 0) {
                    float[] buffer = sta.getFrame(0).getBuffer();
                    variance[i] = MathUtil.mean(buffer); //sta.getFrame(0).getBuffer();
                } else {
                    variance[i] = .0256;
                }
            }

//            PlotUtil.showArray("", mean);
            PlotUtil.showArray("", variance);

        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }


    public static void testFrameFitting() throws Exception {
        STAFile sf = new STAFile("f:\\data\\2005-04-14-0\\data002\\data002.sta");
        STA sta = sf.getSTA(5);
        double stixelWidth = sta.getStixelWidth();
        double stixelHeight = sta.getStixelHeight();
        int w = sta.getWidth();
        int h = sta.getHeight();

        STAFrame frame = new STAFrame(w, h, stixelWidth, stixelHeight);
        for (int k = 0; k <= 29; k++) {
            ImageFrame f = sta.getFrame(k);
            for (int i = 0; i < w; i++) {
                for (int j = 0; j < h; j++) {
                    for (int c = 0; c < 3; c++) {
                        frame.addToPixel(i, j, c, f.getPixel(i, j, c));
                        frame.setPixelError(i, j, c, f.getPixelError(i, j, c));
                    }
                }
            }
        }

//        Gaussian2DFunction g = ImproveSTAFit.fitFrame(frame, 1);
//        System.err.println(g);
//        ParametricEllipse e = new ParametricEllipse(
//            g.getX0() * pixelSize, g.getY0() * pixelSize, g.getSigmaX() * pixelSize,
//            g.getSigmaY() * pixelSize, -g.getTheta());
//
//        STA sta1 = new STA(new ImageFrame[] {frame}, 1);
//        PlotPanel p = STAPlotMaker.makeSTAPanel(sta1, false, 3, false, true, false, e);
//        PlotUtil.showData("", p);
    }


    public static void main(String[] args) throws Exception {
        testFrameFitting();
    }
}
