package edu.ucsc.neurobiology.vision.plot;

import java.util.*;

import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.math.fitting.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class IntegerHistogram
    implements HistogramData, Cloneable {

    protected String description;
    protected double min;
    protected double max;
    protected double binInterval;
    protected double underflowBin;
    protected double overflowBin;
    protected int nNaN = 0;
    protected int[] bins;


    public IntegerHistogram() {
    }


    public IntegerHistogram(
        String description, double min, double max, double binInterval) {

        this.description = description;
        this.min = min;
        this.max = max;
        this.binInterval = binInterval;

        int nBins = (int) Math.ceil( (max - min) / binInterval);
        this.bins = new int[nBins];
    }


    public IntegerHistogram(
        String description, int[] bins, double min, double max) {

        this.description = description;
        this.min = min;
        this.max = max;
        this.binInterval = (max - min) / bins.length;

        this.bins = bins;
    }


    public void redefine(double min, double max) {
        this.min = min;
        this.max = max;
        this.binInterval = (max - min) / bins.length;
    }


    public Gaussian1DFunction fitToGaussian(int nLimit) throws FitFailedException {
        Fitter fitter = new Fitter();

        double A = Double.NEGATIVE_INFINITY;
        double xMax = Double.NaN;
        DoubleList x = new DoubleList();
        DoubleList y = new DoubleList();
        DoubleList err = new DoubleList();
        for (int i = 0; i < bins.length; i++) {
            if (bins[i] > nLimit) {
                final double x0 = min + i * binInterval;
                x.add(x0);
                y.add(bins[i]);
                err.add(Math.sqrt(bins[i]));
                if (bins[i] > A) {
                    A = bins[i];
                    xMax = x0;
                }
            }
        }
        Gaussian1DFunction f = new Gaussian1DFunction(A, xMax, 10);
        fitter.fit(f, new double[][] {x.toArray()}
                   , y.toArray(), err.toArray(), x.size());

        return f;
    }


    public void setBins(double[] newBins) {
        System.arraycopy(newBins, 0, bins, 0, bins.length);
    }


    public double getBinInterval() {
        return binInterval;
    }


    public void scale(double factor) {
        if (Double.isNaN(factor)) {
            throw new NumberFormatException("factor is NaN");
        }
        for (int i = 0; i < bins.length; i++) {
            bins[i] *= factor;
        }
        underflowBin *= factor;
        overflowBin *= factor;
    }


    public void scaleBin(int i, double factor) {
        if (Double.isNaN(factor)) {
            throw new NumberFormatException("factor is NaN");
        }

        bins[i] *= factor;
    }


    public double getBinSum() {
        double binSum = 0;

        for (int i = 0; i < bins.length; i++) {
            binSum += bins[i];
        }

        return binSum;
    }


    public double getMaxValue() {
        double max = Double.NEGATIVE_INFINITY;

        for (int i = 0; i < bins.length; i++) {
            if (bins[i] > max) {
                max = bins[i];
            }
        }

        return max;
    }


    public int getMaxValueBin() {
        double max = Double.NEGATIVE_INFINITY;
        int maxBin = -1;

        for (int i = 0; i < bins.length; i++) {
            if (bins[i] > max) {
                max = bins[i];
                maxBin = i;
            }
        }

        return maxBin;
    }


    public double[] getXValues() {
        double[] x = new double[bins.length];
        for (int i = 0; i < bins.length; i++) {
            x[i] = min + (i + 0.5) * binInterval;
        }
        return x;
    }


    public void fill(double x, double value) {
        if (Double.isNaN(x)) {
            nNaN++;
        } else if (Double.isNaN(value)) {
            throw new NumberFormatException("value is NaN");
        } else if (x < min) {
            underflowBin += value;
        } else if (x >= max) {
            overflowBin += value;
        } else {
            bins[ (int) ( (x - min) / binInterval)] += value;
        }
    }


    public void fill(DoubleHistogram h) {
        int n = h.getBinCount();
        double x0 = h.getMin();
        double dx = h.getBinInterval();

        for (int i = 0; i < n; i++) {
            fill(x0 + (i + 0.5) * dx, h.getBin(i));
        }
    }


    public void fillBin(int bin, double value) {
        if (Double.isNaN(value)) {
            throw new NumberFormatException("value is NaN");
        }
        bins[bin] += value;
    }


    public double getValueAt(double x) {
        if (x < min) {
            return underflowBin;
        } else if (x >= max) {
            return overflowBin;
        } else {
            int bin = (int) Math.floor( (x - min) / binInterval);
            return bins[bin];
        }
    }


    public void clear() {
        Arrays.fill(bins, 0);
    }


    public void setBin(int bin, int value) {
        bins[bin] = value;
    }


    public String getDescription() {
        return description;
    }


    public void setDescription(String description) {
        this.description = description;
    }


    public int getBinCount() {
        return bins.length;
    }


    public double getMin() {
        return min;
    }


    public double getMax() {
        return max;
    }


    public double getBin(int i) {
        if ( (i < 0) || (i >= bins.length)) {
            throw new ArrayIndexOutOfBoundsException("Index: " + i);
        }
        return bins[i];
    }


    public double count(double x1, double x2) {
        int n1 = (int) Math.round( (x1 - min) / binInterval);
        int n2 = (int) Math.round( (x2 - min) / binInterval);
        double n = 0;

        for (int i = n1; i < n2; i++) {
            n += bins[i];
        }

        return n;
    }


    public Object clone() {
        try {
            return super.clone();
        } catch (CloneNotSupportedException e) {
            return null;
        }
    }


    public double[] toArray() {
        double[] newBins = new double[bins.length];
        System.arraycopy(bins, 0, newBins, 0, bins.length);
        return newBins;
    }


    /*
     public Tag read(int TagID, TaggedInputStream input, int len) throws IOException {
            DoubleHistogram h = (DoubleHistogram) clone();
            h.description = input.readString();
            h.min = input.readDouble();
            h.max = input.readDouble();
            h.binInterval = input.readDouble();
            h.underflowBin = input.readDouble();
            h.overflowBin = input.readDouble();
            int binCount = input.readInt();
            h.bins = new double[binCount];
            for (int i = 0; i < binCount; i++) {
                h.bins[i] = input.readDouble();
            }
            return h;
        }
        public void write(int TagID, TaggedOutputStream output) throws IOException {
            output.writeString(description);
            output.writeDouble(min);
            output.writeDouble(max);
            output.writeDouble(binInterval);
            output.writeDouble(underflowBin);
            output.writeDouble(overflowBin);
            output.writeInt(bins.length);
            for (int i = 0; i < bins.length; i++) {
                output.writeDouble(bins[i]);
            }
        }
     */

}
