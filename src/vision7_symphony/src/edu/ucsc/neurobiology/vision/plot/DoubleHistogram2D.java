package edu.ucsc.neurobiology.vision.plot;

import java.io.*;

import java.awt.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class DoubleHistogram2D
    implements Histogram2DData {
    protected String description;

    protected double xMin, xMax, yMin, yMax;
    protected double xBinInterval, yBinInterval;
    //  protected double xUnderflowBin, xOverflowBin, yUnderflowBin, yOverflowBin;
    final int nxBins, nyBins;
    protected double[][] bins;


    public DoubleHistogram2D(
        String description, double xMin, double xMax, double yMin, double yMax,
        double xBinInterval, double yBinInterval) {

        this.description = description;
        this.xMin = xMin;
        this.xMax = xMax;
        this.yMin = yMin;
        this.yMax = yMax;
        this.xBinInterval = xBinInterval;
        this.yBinInterval = yBinInterval;

        nxBins = (int) Math.ceil( (xMax - xMin) / xBinInterval);
        nyBins = (int) Math.ceil( (yMax - yMin) / yBinInterval);
        bins = new double[nxBins][nyBins];
    }


    public DoubleHistogram2D(
        String description, double xMin, double xMax, double yMin, double yMax,
        double[][] bins) {

        this.description = description;
        this.xMin = xMin;
        this.xMax = xMax;
        this.yMin = yMin;
        this.yMax = yMax;
        this.nxBins = bins.length;
        this.nyBins = bins[0].length;
        this.xBinInterval = (xMax - xMin) / nxBins;
        this.yBinInterval = (yMax - yMin) / nyBins;

        this.bins = bins;
    }


    public double getMinX() {
        return xMin;
    }


    public double getMaxX() {
        return xMax;
    }


    public double getMinY() {
        return yMin;
    }


    public double getMaxY() {
        return yMax;
    }


    public double getBinIntervalX() {
        return xBinInterval;
    }


    public double getBinIntervalY() {
        return yBinInterval;
    }


    public double getBin(int x, int y) {
        return bins[x][y];
    }


    public Point getMaxValueBin() {
        Point p = new Point(0, 0);
        double maxValue = bins[0][0];

        for (int i = 0; i < nxBins; i++) {
            for (int j = 0; j < nyBins; j++) {
                if (bins[i][j] > maxValue) {
                    maxValue = bins[i][j];
                    p.x = i;
                    p.y = j;
                }
            }
        }

        return p;
    }


    public Point getMaxAbsValueBin() {
        Point p = new Point(0, 0);
        double maxValue = Double.NEGATIVE_INFINITY;

        for (int i = 0; i < nxBins; i++) {
            for (int j = 0; j < nyBins; j++) {
                if (Math.abs(bins[i][j]) > maxValue) {
                    maxValue = Math.abs(bins[i][j]);
                    p.x = i;
                    p.y = j;
                }
            }
        }

        return p;
    }


    public int getBinCountX() {
        return nxBins;
    }


    public int getBinCountY() {
        return nyBins;
    }


//    public int getBinCount() {
//        return getBinCountX() * getBinCountY();
//    }


    /*
        public double getBinInterval() {
            return binInterval;
        }
     */

    /*
        public void scale(double factor) {
            for (int i = 0; i < bins.length; i++) {
                bins[i] *= factor;
            }
            underflowBin *= factor;
            overflowBin *= factor;
        }
     */

    /*
        public double getBinSum() {
            double binSum = 0;
            for (int i = 0; i < bins.length; i++) {
                binSum += bins[i];
            }
            return binSum;
        }
     */

    public void fill(double x, double y, double value) {
        int xBin = (int) Math.floor( (x - xMin) / xBinInterval);
        int yBin = (int) Math.floor( (y - yMin) / yBinInterval);
//        System.out.println(xBin);
//        System.out.println(yBin);
        try {
            bins[xBin][yBin] += value;
        } catch (ArrayIndexOutOfBoundsException e) {
//            System.out.println("fill - problem");
        }

        /*
                 if (x < min) {
            underflowBin += value;
                 } else  if (x >= max) {
            overflowBin += value;
                 } else {
            int bin = (int)Math.floor((x - min)/binInterval);
            bins[bin] += value;
                 }
         */
    }


    public Point getBinForCoordinates(double x, double y, Point b) {
        b.x = (int) Math.floor( (x - xMin) / xBinInterval);
        b.y = (int) Math.floor( (y - yMin) / yBinInterval);

        return b;
    }


    public void fillBin(int xBin, int yBin, double value) {
        bins[xBin][yBin] += value;
    }


    public void setBin(int xBin, int yBin, double value) {
        bins[xBin][yBin] = value;
    }


    public double getBinValueForCoords(double x, double y) {
        int xBin = (int) Math.floor( (x - xMin) / xBinInterval);
        int yBin = (int) Math.floor( (y - yMin) / yBinInterval);
        return bins[xBin][yBin];
    }


    public void save(String fileName) throws IOException {
        PrintWriter writer = new PrintWriter(new FileWriter(fileName));

        for (int j = 0; j < bins[0].length; j++) {
            for (int i = 0; i < bins.length; i++) {
                writer.print(bins[i][j] + "\t");
            }
            writer.println();
        }

        writer.close();
    }


    /*
        public void fillBin(int bin, double value) {
            bins[bin] += value;
        }
     */

    /*
        public double getValueAt(double x) {
            if (x < min) {
                return underflowBin;
            } else  if (x >= max) {
                return overflowBin;
            } else {
                int bin = (int)Math.floor((x - min)/binInterval);
                return bins[bin];
            }
        }
     */

    /*
        public void clear() {
            Arrays.fill(bins, 0);
        }
     */
    /*
        public void setBin(int bin, double value) {
            bins[bin] = value;
        }
     */

    public String getDescription() {
        return description;
    }


    public void setDescription(String description) {
        this.description = description;
    }


    /*
        public double getMin() {
            return min;
        }
        public double getMax() {
            return max;
        }
     */
    /*
        public double getBin(int i) {
            if ((i < 0)||(i >= bins.length))
                throw new ArrayIndexOutOfBoundsException("Index: " + i);
            return bins[i];
        }
     */

}
