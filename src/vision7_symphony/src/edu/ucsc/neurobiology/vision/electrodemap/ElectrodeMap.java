package edu.ucsc.neurobiology.vision.electrodemap;

import java.util.*;

import java.awt.geom.*;

import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * Class which defines an electrode map.  Electrodes are numbered
 * from zero to (getNumberOfElectrodes()-1).  The electrode map
 * provides information about the position of each electrode and
 * their adjacency. The positions
 * of all electrodes must be characterized by two integers.
 * Subclasses of ElectrodeMap may be able to cut the map in a number of submaps
 * as appropriate.
 *
 * @see getElectrodeMap(int part, int nParts)
 * @author Charles A. Loomis, University of California, Santa Cruz
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ElectrodeMap implements PlotData {

    /** The number of electrodes in this map. */
    protected int nElectrodes;

    /** The x-position of each electrode. */
    protected float[] xPosition;

    /** The y-position of each electrode. */
    protected float[] yPosition;

    protected boolean[] isDisconnected;

    protected int[] parentElectrodeNumbers = null;

    protected double pitch = 1;

    public double micronsPerPixelX, micronsPerPixelY, centerX, centerY, angle;
    public boolean flipX, flipY;
    protected double xAdjacencyScale, yAdjacencyScale;



    public ElectrodeMap(float[] xCoord, float[] yCoord, int[] parentElectrodeNumbers) {
        this(xCoord, yCoord);
        this.parentElectrodeNumbers = parentElectrodeNumbers;
    }


    /**
     * This constructor contains all the information necessary for the
     * electrode map.  Generally this is called by some utility routine.
     */
    public ElectrodeMap(float[] xCoord, float[] yCoord) {
        if (xCoord.length != yCoord.length) {
            throw new IllegalArgumentException("xCoord.length != yCoord.length");
        }

        // Check that the sizes of the arrays are OK.
        if (xCoord.length == 0) {
            throw new IllegalArgumentException("Coordinate arrays are too short (0).");
        }

        this.nElectrodes = xCoord.length;

        // Cache the actual positions (given in microns).
        xPosition = new float[nElectrodes];
        yPosition = new float[nElectrodes];
        for (int i = 0; i < nElectrodes; i++) {
            xPosition[i] = (float) (xCoord[i] /* * xPitch*/);
            yPosition[i] = (float) (yCoord[i] /* * yPitch*/);
        }
        if (parentElectrodeNumbers == null) {
            parentElectrodeNumbers = new int[nElectrodes];
            for (int i = 0; i < nElectrodes; i++) {
                parentElectrodeNumbers[i] = i;
            }
        }

        this.isDisconnected = new boolean[nElectrodes];
    }


    /**
     * Return the number of electrodes in this map.
     */
    public int getNumberOfElectrodes() {
        return nElectrodes;
    }


    public void setDisconnected(int[] id) {
        for (int i = 0; i < id.length; i++) {
            isDisconnected[id[i]] = true;
        }
    }


    /**
     * Return the position of the given electrode.  This will reuse the given
     * point; if this is null, then a new Point2D will be created and
     * returned.
     */
    public Point2D getPosition(int electrode, Point2D point) {
        if (point == null) {
            point = new Point2D.Double();
        }
        point.setLocation(xPosition[electrode], yPosition[electrode]);
        return point;
    }


    /**
     * Return whether or not the two electrodes are adjacent. */
    //	public abstract boolean areAdjacent(int electrode1, int electrode2);


    public boolean[] getDisconnectedElectrodesList() {
        return isDisconnected;
    }


    public float getXPosition(int electrode) {
        return xPosition[electrode];
    }


    public float getYPosition(int electrode) {
        return yPosition[electrode];
    }


    public int countDisconnectedElectrodes() {
        return MathUtil.countValues(true, getDisconnectedElectrodesList());
    }


    public boolean isDisconnected(int electrode) {
        if (electrode == 0 || getDisconnectedElectrodesList()[electrode]) {
            return true;
        } else {
            return false;
        }
    }


    /**
     * Return whether or not the two electrodes are near. */
    //	 public abstract boolean areNear(int electrode1, int electrode2);


    public int[] getAdjacentsTo(int electrode) {
        return getAdjacentsTo(electrode, 1);
    }



    public double getPitch() {
        return pitch;
    }


    public int[] getAdjacentsTo(int electrode, int radius) {
        double a;
        double b;
        int n = getNumberOfElectrodes();
        if (n == 68) {
            a = xAdjacencyScale * getPitch();
            b = yAdjacencyScale * getPitch();
        } else {
            a = xAdjacencyScale * radius * getPitch();
            b = yAdjacencyScale * radius * getPitch();
        }

        Point2D.Double p0 = new Point2D.Double();
        getPosition(electrode, p0);
        
        IntegerList list = new IntegerList();
        list.add(electrode);
        Point2D.Double p = new Point2D.Double();
        for (int i = 1; i < n; i++) {
            getPosition(i, p);
            if (i != electrode && !isDisconnected(i) &&
                    Math.pow( (p.x - p0.x) / a, 2) + Math.pow( (p.y - p0.y) / b, 2) <= 1.01) {
                list.add(i);
            }
        }

        return list.toArray();
    }


    /**
     * Returns the bounds of the array in the format: xMin, xMax, yMin, yMax
     */
    public double[] getBounds() {
        double xMin = getXPosition(1);
        double xMax = getXPosition(1);
        double yMin = getYPosition(1);
        double yMax = getYPosition(1);

        int n = getNumberOfElectrodes();
        for (int i = 1; i < n; i++) {
            double x = getXPosition(i);
            double y = getYPosition(i);
            if (x < xMin) {
                xMin = x;
            }
            if (x > xMax) {
                xMax = x;
            }
            if (y < yMin) {
                yMin = y;
            }
            if (y > yMax) {
                yMax = y;
            }
        }

        return new double[] {xMin, xMax, yMin, yMax};
    }


    public String getDescription() {
        return "Electrode Map";
    }


    public double getDistance(int i, int j) {
        Point2D.Double pi = new Point2D.Double();
        Point2D.Double pj = new Point2D.Double();
        getPosition(i, pi);
        getPosition(j, pj);

        return Math.sqrt( (pi.x - pj.x) * (pi.x - pj.x) + (pi.y - pj.y) * (pi.y - pj.y));
    }


    public String toString() {
        String s = "Electrode Map [" + "nElectrodes=" + nElectrodes;

        for (int i = 0; i < nElectrodes; i++) {
            s += "\n" + xPosition[i] + ", " + yPosition[i];
        }

        return s;
    }


    public int[] getParentElectrodeNumbers() {
        return parentElectrodeNumbers;
    }

    /**
     * Measures and corrects for differences between electrode array and display computer positioning and scaling
     * Flip, then rotation, then displacement.
     * 
     * @param xCorners double[]  X-coordinates from corners of electrode map, in pixels.
     * @param yCorners double[]  Y-coordinates from cornes of electrode map, in pixels.
     * @param flipX boolean  Are the optics flipped in the X direction?
     * @param flipX boolean  Are the optics flipped in the Y direction?
     * 
     */
    public void SetAlignment(int[] xCorners, int[] yCorners, boolean flipX, boolean flipY) throws Exception {

        if(xCorners == null) {
            SetAlignment(5.8, 5.8, Double.NaN, Double.NaN, flipX, flipY, 0.0);
            return;
        }


        //clip hexagon down to a rectangle.  Assumes that bottom and top of hexagon are horizontal.
        if(this instanceof Hexagonal61ElectrodeMap) {
            double mean = MathUtil.mean(yCorners);
            double max = MathUtil.max(yCorners);
            double halfHeight = Math.abs(max-mean);
            int[] xCornersNew = new int[4];
            int[] yCornersNew = new int[4];
            int currentCorner = 0;

            for(int i=0; i<xCorners.length; i++) {
                if(Math.abs(yCorners[i] - mean) > halfHeight/2) {
                    xCornersNew[currentCorner] = xCorners[i];
                    yCornersNew[currentCorner] = yCorners[i];
                    currentCorner++;
                }
            }
            xCorners = xCornersNew;
            yCorners = yCornersNew;
        }
        //clip hexagon down to a rectangle.  Assumes that bottom and top of hexagon are a single point.		
        if(this instanceof Hexagonal519ElectrodeMap) {

            double max = MathUtil.max(yCorners);
            double min = MathUtil.min(yCorners);

            int[] xCornersNew = new int[4];
            int[] yCornersNew = new int[4];
            int currentCorner = 0;

            for(int i=0; i<xCorners.length; i++) {
                if((yCorners[i] != max) && (yCorners[i] != min)) {
                    xCornersNew[currentCorner] = xCorners[i];
                    yCornersNew[currentCorner] = yCorners[i];
                    currentCorner++;
                }
            }
            xCorners = xCornersNew;
            yCorners = yCornersNew;
        }
        

        //Used to determine if a segment is oriented in the x direction
        double xTest = 0;
        for(int i=0; i<xCorners.length; i++) {
            if(xCorners[(i+1)%xCorners.length] - xCorners[i] > xTest) {
                xTest = xCorners[(i+1)%xCorners.length] - xCorners[i];
            }
        }
        xTest/=2;



        double angle = 0.0; //in radians,  positive is CCW
        double pixelWidth = 0.0, pixelHeight = 0.0;  // in microns/pixel
        double centerX = 0.0, centerY = 0.0; // in um


        //invert y axis.  screen origin is in upper left.  physical origin is in lower left.
        for(int i=0; i<xCorners.length; i++) {
            yCorners[i] = 480 - yCorners[i];
        }

        double polygonArea = 0.0;
        for(int i=0; i<xCorners.length; i++) {
            polygonArea+= (xCorners[i]*yCorners[(i+1)%yCorners.length]
                                                -xCorners[(i+1)%xCorners.length]*yCorners[i]);
        }
        polygonArea/=2;
        //if vertices are ordered CW, switch to CCW.  Angle calculation assumes this.
        if(polygonArea < 0) {
            int[] tempX = new int[xCorners.length];
            int[] tempY = new int[yCorners.length];
            for(int i=0; i<xCorners.length; i++) {
                tempX[i] = xCorners[xCorners.length-1-i];
                tempY[i] = yCorners[yCorners.length-1-i];
            }
            for(int i=0; i<xCorners.length; i++) {
                xCorners[i] = tempX[i];
                yCorners[i] = tempY[i];
            }
            polygonArea = -polygonArea;
        }




        for (int i = 0; i < xCorners.length; i++) {

            double x = (xCorners[(i+1)%xCorners.length]-xCorners[i]);
            double y = (yCorners[(i+1)%yCorners.length]-yCorners[i]);
            double hyp = Math.sqrt(
                    Math.pow(x, 2)+ 
                    Math.pow(y,2));
            //if x oriented segment
            if(x > xTest || x < -xTest) {
                angle += ((x > 0) ? 1 : -1) * Math.asin(y/hyp);				
            } else {
                angle += ((y > 0) ? -1 : 1) * Math.asin(x/hyp);
            }


            centerX+=xCorners[i];
            centerY+=yCorners[i];

        }
        angle/=xCorners.length;
        centerX/=xCorners.length;
        centerY/=yCorners.length;
        int xMeasurements=0, yMeasurements=0;

        for (int i=0; i<xCorners.length; i++) {
            double x = Math.abs((xCorners[(i+1)%xCorners.length]-xCorners[i]));
            double y = Math.abs((yCorners[(i+1)%yCorners.length]-yCorners[i]));
//			double hyp = Math.sqrt(
//					Math.pow(x, 2)+ 
//					Math.pow(y,2));
            //if x oriented segment
            if(x > xTest) {

                if(this instanceof Rectangular512ElectrodeMap) {
                    pixelWidth += 31.5*pitch*Math.cos(angle)/x;
                }
                else if(this instanceof Hexagonal61ElectrodeMap ) {
                    pixelWidth += 4 * pitch * Math.cos(angle)/x;
                } else if(this instanceof Hexagonal519ElectrodeMap ) {
                    pixelWidth += 24 * pitch * Math.cos(angle)/x;
                }
                else {
                    throw new IllegalStateException("Set Alignment failed.  Requested electrode map not implemented.");
                }
                xMeasurements++;

            } else {
                if(this instanceof Rectangular512ElectrodeMap) {
                    pixelHeight += 15*pitch*Math.cos(angle)/y;
                }
                else if(this instanceof Hexagonal61ElectrodeMap ) {
                    pixelHeight += 8*pitch*Math.cos(angle)/y;
                } 	else if(this instanceof Hexagonal519ElectrodeMap ) {
                    pixelHeight += 14*pitch*Math.cos(angle)/y;
                } else {
                    throw new IllegalStateException("Set Alignment failed.  Requested electrode map not implemented.");
                }
                yMeasurements++;

            }
        }


        pixelWidth/= (xMeasurements);
        pixelHeight/= (yMeasurements);

        centerX *= pixelWidth;
        centerY *= pixelHeight;


        SetAlignment(pixelWidth, pixelHeight, centerX, centerY, flipX, flipY, angle);

    }





    public void SetAlignment(double micronsPerPixelX,double micronsPerPixelY, double centerX,
            double centerY, boolean flipX, boolean flipY, double angle) {
        this.micronsPerPixelX = micronsPerPixelX;
        this.micronsPerPixelY = micronsPerPixelY;
        this.centerX = centerX;
        this.centerY = centerY;
        this.flipX = flipX;
        this.flipY = flipY;
        this.angle = angle;

        flipAndRotate(flipX, flipY, angle, centerX, centerY);
    }

    //flip occurs before rotation
    /**
     * @param flipX
     * @param flipY
     * @param angle
     * @param centerX
     * @param centerY
     * 
     * NOTE: If you make any change here, be sure to make matching change to Polygon2D#flipAndRotate, or 519
     * streaming will break.
     * 
     * FIXME: Should actually use THE SAME code to do both, so make an abstracted call or encapsulate the x and y points
     * in both objects into a Coordinates object that has this method...
     * 
     */
    public void flipAndRotate(boolean flipX, boolean flipY, double angle, double centerX, double centerY) {
        for (int i = 0; i < nElectrodes; i++) {
            xPosition[i] = flipX ? - xPosition[i] : xPosition[i];
            yPosition[i] = flipY ? - yPosition[i] : yPosition[i];
        }

        double newX, newY;
        for (int i = 0; i < nElectrodes; i++) {
            //modifies variables in super.  Local xCoord, yCoord are final "pure" definition.
            newX = xPosition[i] * Math.cos(angle) - yPosition[i] * Math.sin(angle) + centerX;
            newY = xPosition[i] * Math.sin(angle) + yPosition[i] * Math.cos(angle) + centerY;
            xPosition[i] = (float) newX;
            yPosition[i] = (float) newY;
        }		 
    }
    
    /**
     * Cuts the map into "nParts" submaps and returns the "part" submap.
     * the method getRegionShape() is used to do the cutting. See specific classes
     * for allowed cuts.
     *
     * @see getRegionShape(int part, int nParts)
     * @param part int  from 1 to nParts
     * @param nParts int the number of parts this map is to be split in
     * @return ElectrodeMap
     */
    public ElectrodeMap getElectrodeMap(int part, int nParts) {
        if (nParts == 0) {
            throw new Error("nParts cannot be zero");
        }
        if (part <= 0 || part > nParts) {
            throw new Error("wrong 'part' number " + part);
        }

        Polygon2D region = getRegionShape(part, nParts);

        IntegerList electrodeNumberList = new IntegerList();
        float[] x = new float[nElectrodes];
        float[] y = new float[nElectrodes];
        boolean[] disc = new boolean[nElectrodes];
        x[0] = xPosition[0];
        y[0] = yPosition[0];

        // create the actual electrode map
        int n = 1;
        for (int i = 1; i < nElectrodes; i++) {
            if (region.contains(xPosition[i], yPosition[i])) {
                x[n] = xPosition[i];
                y[n] = yPosition[i];
                disc[n] = isDisconnected[i];
                electrodeNumberList.add(i);
                n++;
            }
        }

        int[] electrodeNumbers = electrodeNumberList.toArray();
        Arrays.sort(electrodeNumbers);

        float[] _x = new float[n];
        float[] _y = new float[n];
        boolean[] _disc = new boolean[n];
        for (int i = 0; i < n; i++) {
            _x[i] = x[i];
            _y[i] = y[i];
            _disc[i] = disc[i];
        }

        ElectrodeMap map = new ElectrodeMap(_x, _y, electrodeNumbers);
        map.isDisconnected = _disc;

        return map;
    }


    public Polygon2D getRegionShape(int part, int nParts) {
        /*
         * Vincent Deo, 02/01/2016
         * Making a hack here to be able to use single electrode map with ID 9999
         */
        if (nElectrodes != 2) {
            // Historical version
            throw new Error("Regions are not supported for this electrode array");
        }
        if (!(part == 1 && nParts == 1)) {
            throw new Error("!(part == 1 && nParts == 1)");
        }
        Polygon2D p = new Polygon2D();
        p.addPoint(1, 0);
        p.addPoint(0, 1);
        p.addPoint(-6, -5);
        p.addPoint(-5, -6);
        return p;
    }

}
