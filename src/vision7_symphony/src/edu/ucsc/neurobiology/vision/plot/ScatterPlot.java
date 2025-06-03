package edu.ucsc.neurobiology.vision.plot;

import java.util.*;

import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ScatterPlot
    implements ScatterPlotData {

    private static int scatteplotCount = 1;

    private DoubleList xList, yList, xErrList, yErrList;
    private String description;
    private boolean hasYErrors, hasXErrors;


    public ScatterPlot(String description) {
        this.description = description;
        this.xList = new DoubleList();
        this.yList = new DoubleList();
        this.xErrList = new DoubleList();
        this.yErrList = new DoubleList();

        hasXErrors = true;
        hasYErrors = true;
    }


    public ScatterPlot() {
        this("ScatterPlot " + scatteplotCount++);
    }


    public ScatterPlot(double[] y) {
        this();

        this.xList = new DoubleList();
        this.yList = new DoubleList();

        for (int i = 0; i < y.length; i++) {
            xList.add(i);
            yList.add(y[i]);
        }

        hasXErrors = false;
        hasYErrors = false;
    }


    public ScatterPlot(double[] x, double[] y, double[] yError) {
        this(x, y, yError, null);
    }


    public ScatterPlot(double[] x, double[] y, double[] yError, double[] xError) {
        this();

        this.xList = new DoubleList();
        this.yList = new DoubleList();

        if (y == null) {
            throw new IllegalArgumentException("y cannot be null");
        } else {
            for (int i = 0; i < y.length; i++) {
                yList.add(y[i]);
            }
        }

        if (x == null) {
            for (int i = 0; i < y.length; i++) {
                xList.add(i);
            }
        } else {
            for (int i = 0; i < y.length; i++) {
                xList.add(x[i]);
            }
        }

        if (xError == null) {
            hasXErrors = false;
        } else {
            hasXErrors = true;
            this.xErrList = new DoubleList();
            for (int i = 0; i < y.length; i++) {
                xErrList.add(xError[i]);
            }
        }

        if (yError == null) {
            hasYErrors = false;
        } else {
            hasYErrors = true;
            this.yErrList = new DoubleList();
            for (int i = 0; i < y.length; i++) {
                yErrList.add(yError[i]);
            }
        }
    }


    public ScatterPlot(float[] x, float[] y, float[] yError, float[] xError) {
        this();

        this.xList = new DoubleList();
        this.yList = new DoubleList();

        if (y == null) {
            throw new IllegalArgumentException("y cannot be null");
        } else {
            for (int i = 0; i < y.length; i++) {
                yList.add(y[i]);
            }
        }

        if (x == null) {
            for (int i = 0; i < y.length; i++) {
                xList.add(i);
            }
        } else {
            for (int i = 0; i < y.length; i++) {
                xList.add(x[i]);
            }
        }

        if (xError == null) {
            hasXErrors = false;
        } else {
            hasXErrors = true;
            this.xErrList = new DoubleList();
            for (int i = 0; i < y.length; i++) {
                xErrList.add(xError[i]);
            }
        }

        if (yError == null) {
            hasYErrors = false;
        } else {
            hasYErrors = true;
            this.yErrList = new DoubleList();
            for (int i = 0; i < y.length; i++) {
                yErrList.add(yError[i]);
            }
        }

    }


    public ScatterPlot(short[] x, short[] y, short[] yError, short[] xError) {
        this();

        this.xList = new DoubleList();
        this.yList = new DoubleList();

        if (y == null) {
            throw new IllegalArgumentException("y cannot be null");
        } else {
            for (int i = 0; i < y.length; i++) {
                yList.add(y[i]);
            }
        }

        if (x == null) {
            for (int i = 0; i < y.length; i++) {
                xList.add(i);
            }
        } else {
            for (int i = 0; i < y.length; i++) {
                xList.add(x[i]);
            }
        }

        if (xError == null) {
            hasXErrors = false;
        } else {
            hasXErrors = true;
            this.xErrList = new DoubleList();
            for (int i = 0; i < y.length; i++) {
                xErrList.add(xError[i]);
            }
        }

        if (yError == null) {
            hasYErrors = false;
        } else {
            hasYErrors = true;
            this.yErrList = new DoubleList();
            for (int i = 0; i < y.length; i++) {
                yErrList.add(yError[i]);
            }
        }

    }


    public ScatterPlot(int[] x, int[] y, int[] yError, int[] xError) {
        this();

        this.xList = new DoubleList();
        this.yList = new DoubleList();

        if (y == null) {
            throw new IllegalArgumentException("y cannot be null");
        } else {
            for (int i = 0; i < y.length; i++) {
                yList.add(y[i]);
            }
        }

        if (x == null) {
            for (int i = 0; i < y.length; i++) {
                xList.add(i);
            }
        } else {
            for (int i = 0; i < y.length; i++) {
                xList.add(x[i]);
            }
        }

        if (xError == null) {
            hasXErrors = false;
        } else {
            hasXErrors = true;
            this.xErrList = new DoubleList();
            for (int i = 0; i < y.length; i++) {
                xErrList.add(xError[i]);
            }
        }

        if (yError == null) {
            hasYErrors = false;
        } else {
            hasYErrors = true;
            this.yErrList = new DoubleList();
            for (int i = 0; i < y.length; i++) {
                yErrList.add(yError[i]);
            }
        }

    }


    public boolean hasXErrors() {
        return hasXErrors;
    }


    public boolean hasYErrors() {
        return hasYErrors;
    }


    public String getDescription() {
        return description;
    }


    public void setDescription(String description) {
        this.description = description;
    }


    public ArrayList<String> getLegend() {
        ArrayList<String> legend = new ArrayList<String>();

        legend.add("nPoints: " + xList.size());

        return legend;
    }


    /**
     * The error is assumed to be zero.
     */
    public void add(double x, double y) {
        xList.add(x);
        yList.add(y);
        xErrList.add(0);
        yErrList.add(0);
    }


    public void add(double x, double y, double yErr) {
        xList.add(x);
        yList.add(y);
        xErrList.add(0);
        yErrList.add(yErr);
    }


    public void add(double x, double y, double yErr, double xErr) {
        xList.add(x);
        yList.add(y);
        xErrList.add(xErr);
        yErrList.add(yErr);
    }


    public void clear() {
        xList.clear();
        yList.clear();
        xErrList.clear();
        yErrList.clear();
    }


    public void setDrawingControl(DrawingControl drawingControl) {
    }


    public int getPointCount() {
        return xList.size();
    }
    
    public int getPointSize() {
        return yList.size();
    }


    public void getDataPoint(int i, double[] point) {
        point[0] = xList.get(i);
        point[1] = yList.get(i);
        point[2] = (hasYErrors) ? yErrList.get(i) : 0;
        point[3] = (hasXErrors) ? xErrList.get(i) : 0;
    }


    public void maxOne() {
        double max = yList.get(0);
        for (int i = 0; i < yList.size(); i++) {
            if (yList.get(i) > max) {
                max = yList.get(i);
            }
        }

        scale(1 / max);
    }


    public void scale(double factor) {
        for (int i = 0; i < yList.size(); i++) {
            yList.set(i, yList.get(i) * factor);
        }
    }


    public void print() {
        for (int i = 0; i < xList.size(); i++) {
            System.err.println(xList.get(i) + "\t" + yList.get(i));
        }
    }
}
