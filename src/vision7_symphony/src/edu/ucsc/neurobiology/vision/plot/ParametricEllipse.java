package edu.ucsc.neurobiology.vision.plot;

import java.util.*;

import java.awt.geom.*;

import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ParametricEllipse implements ParametricFunctionData, SelectionRegion {
    
    private double[] point = new double[2];
    private double x0, y0, a, b, theta;
    private String description;
    private double cosT, sinT;
    private double scaleX = 1, scaleY = 1;
    private String label = null;
    
    /**
     * Creates a new ParametricEllipse class with given parameters.
     *
     * @param x0 X location of the center
     * @param y0 Y location of the center
     * @param a the X semiaxis
     * @param b the Y semiaxis
     * @param theta the rotation angle in radians
     */
    public ParametricEllipse(double x0, double y0, double a, double b, double theta) {
        setParameters(x0, y0, a, b, theta, 1, 1);
    }
    
    /**
     * Creates a new ParametricEllipse class with given parameters.
     *   
     * Must use built in scaling functions for parametric ellipse.
     * Otherwise, asymmetry between scaleX and scaleY will break the function.
     * Major and minor axes do not point in the direction of x and y.
     * 
     * Most commonly, Gaussian fit is in units of stixels, while we 
     * wish to display in microns.
     * 
     * @param x0 X location of the center
     * @param y0 Y location of the center
     * @param a the a semiaxis
     * @param b the b semiaxis
     * @param theta the rotation angle in radians
     * @param scaleX  X scale factor
     * @param scaleY  Y scale factor
     */    
    public ParametricEllipse(double x0, double y0, double a, double b, double theta, 
            double scaleX, double scaleY) {
        setParameters(x0, y0, a, b, theta, scaleX, scaleY);
    }
    
    
    public void setParameters(double x0, double y0, double a, double b, double theta,
            double scaleX, double scaleY) {
        this.x0 = x0;
        this.y0 = y0;
        this.a = a;
        this.b = b;
        this.theta = theta;
        this.scaleX = scaleX;
        this.scaleY = scaleY;

        cosT = Math.cos(theta);
        sinT = Math.sin(theta);
        this.theta = theta;
    }


    public ArrayList<String> getLegend() {
        ArrayList<String> legend = new ArrayList<String>();
        legend.add("x0 = " + StringUtil.format(x0, 5));
        legend.add("y0 = " + StringUtil.format(y0, 5));
        legend.add("a = " + StringUtil.format(a, 5));
        legend.add("b = " + StringUtil.format(b, 5));
        legend.add("theta = " + StringUtil.format(theta, 5));
        return legend;
    }


    public ScatterPlot getCardinalPoints() {
        ScatterPlot sp = new ScatterPlot();
//        for (int i = 0; i < 2; i++) {
//            double[] p = getPointFor(i * Math.PI / 2);
//            sp.add(p[0], p[1]);
//        }
        double dx = Math.sqrt(a * a * cosT * cosT + b * b * sinT * sinT);
        double dy = Math.sqrt(a * a * sinT * sinT + b * b * cosT * cosT);
        sp.add((x0 + dx)*scaleX, y0*scaleY);
        sp.add((x0 - dx)*scaleX, y0*scaleY);
        sp.add(x0*scaleX, (y0 + dy)*scaleY);
        sp.add(x0*scaleX, (y0 - dy)*scaleY);
        return sp;
    }


    public double[] getPointFor(double phi) {
        double x = a * Math.cos(phi);
        double y = b * Math.sin(phi);
        // this is the correct rotation
        point[0] = (x0 + x * cosT - y * sinT)*scaleX;
        point[1] = (y0 + y * cosT + x * sinT)*scaleY;
        return point;
    }


    public double getMinParamValue() {
        return 0;
    }


    public double getMaxParamValue() {
        return 2 * Math.PI;
    }


    public String getDescription() {
        return description;
    }


    public double getMinimumDistance(double x, double y) {
        return 0;
    }


    public boolean contains(double x, double y) {
        x/=scaleX;
        y/=scaleY;
        double sin = Math.sin(theta);
        double cos = Math.cos(theta);
        double d = Math.pow( ( (x - x0) * cos+ (y - y0) * sin) / a, 2) +
                   Math.pow( ( (y - y0) * cos - (x - x0) * sin) / b, 2);
        return (d < 1) ? true : false;
    }


    public Rectangle2D getBounds2D(double userScale) {
        double cos = Math.cos(theta);
        double sin = Math.sin(theta);
        double dx = Math.sqrt(a * a * cos * cos + b * b * sin * sin);
        double dy = Math.sqrt(a * a * sin * sin + b * b * cos * cos);
        return new Rectangle2D.Double((x0-dx*userScale)*scaleX, (y0-dy*userScale)*scaleY, 2*dx*scaleX*userScale, 2*dy*scaleY*userScale);
    }
    
    public Rectangle2D getBounds2D() {
        return getBounds2D(1.0);
    }
    
    /**
     * Expands one dimension to meet targetAspectRatio.  If targetAspectRatio < 0 then do nothing.
     * @param userScale
     * @param targetAspectRatio
     * @return
     */
    public Rectangle2D getBounds2D(double userScale, double targetAspectRatio) {
        Rectangle2D bounds2D = getBounds2D(userScale);
        if (targetAspectRatio < 0) return bounds2D;
        
        // Expand dimensions as necessary
        double w = bounds2D.getWidth();
        double h = bounds2D.getHeight();
        double fixW = Math.max(w, h*targetAspectRatio);
        double fixH = Math.max(h, w/targetAspectRatio);
        
        // Calculate new top left
        double fixX = bounds2D.getCenterX() - fixW/2;
        double fixY = bounds2D.getCenterY() - fixH/2;
        
        // Calculate new top left corner
        return new Rectangle2D.Double(fixX, fixY, fixW, fixH);
    }
    

    public String toString() {
        return
            "x0 " + x0 + "\n" +
            "y0 " + y0 + "\n" +
            "a " + a + "\n" +
            "b " + b + "\n" +
            "theta " + theta + "\n";
    }
    
    public double getX0() {
        return x0;
    }
    
    public double getY0() {
        return y0;
    }
    
    public double getA() {
        return a;
    }
    
    public double getB() {
        return b;
    }
    
    public double getTheta() {
        return theta;
    }
    
    public void setLabel(String l) {
        label = l;
    }
    
    public String getLabel() {
        return label;
    }
    
    public double[] getLabelPosition() {
        return new double[]{x0*scaleX, y0*scaleY};
    }
}