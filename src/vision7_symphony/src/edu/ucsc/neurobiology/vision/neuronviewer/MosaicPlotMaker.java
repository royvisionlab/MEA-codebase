package edu.ucsc.neurobiology.vision.neuronviewer;

import static java.lang.Math.round;



import java.awt.*;
import java.io.File;
import java.util.*;

import javax.swing.*;
import javax.swing.tree.*;

import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;
import edu.ucsc.neurobiology.vision.io.chunk.*;


/**
 * @author Matthew Grivich, The Salk Institute
 */
public class MosaicPlotMaker extends CompoundPlotMaker {

    public boolean showGaussians = true;
    public boolean showIDs = false;
    public boolean showEIGaussians = true;
    public boolean showDisplacements = true;
    public boolean showContours = false;
    public boolean showLegend = true;
    public boolean showAxes = false;
    
    private boolean movieExists, calParamsExists;
    private double micronsPerPixelX, micronsPerPixelY;
    private double centerX, centerY; //microns
    private boolean flipX, flipY; 
    private double angle;//radians		
    
    private FunctionStyle gaussianFitStyle =
        new FunctionStyle("STA Gaussian Fit", Color.black, 1 / 2f);

    private  FunctionStyle gaussianFitStyleEI =
        new FunctionStyle("EI Gaussian Fit", Color.red, 1 / 2f);

    private FunctionStyle simpleContour = new FunctionStyle("Simple Contour", Color.black, .5f);
    private FunctionStyle fullContour = new FunctionStyle("Contour", Color.red, .5f);

    public MosaicPlotMaker() {
        super("Mosaic Plot", CLASS_PLOT);
    }


    public Component makePlot(IntegerList list, int plotType, TreePath classPath) {
        if (list.size() == 0) return new JLabel("No neurons selected");
        final int[] neurons = list.toArray();

        try {
            movieExists = viewer.globalsFile.runTimeMovieParamsExists();
            calParamsExists = viewer.globalsFile.imageCalibrationParamsExists();
        } catch (Exception e) {
            e.printStackTrace();
            return new JLabel("Could not make mosaic plot.");
        }

        try {
            GlobalsFile.ImageCalibrationParams calParams = null;
            if (calParamsExists) {
                calParams = viewer.globalsFile.getImageCalibrationParams();
                micronsPerPixelX = calParams.micronsPerPixelX;
                micronsPerPixelY = calParams.micronsPerPixelY;
                centerX = calParams.centerX;
                centerY = calParams.centerY;
                flipX = calParams.flipX;
                flipY = calParams.flipY;
                angle = calParams.angle;
            } else {
                micronsPerPixelX = 5.8;
                micronsPerPixelY = 5.8;
                centerX = Double.NaN;
                centerY = Double.NaN;
                flipX = false;
                flipY = false;
                angle = 0;

            }
            GlobalsFile.RunTimeMovieParams runParams = null;
            if (movieExists) runParams = viewer.globalsFile.getRunTimeMovieParams();

            HashMap<Integer, Double> xMap = new HashMap<Integer, Double>();			
            HashMap<Integer, Double> yMap = new HashMap<Integer, Double>();
            HashMap<Integer, Double> sigmaXMap = new HashMap<Integer, Double>();
            HashMap<Integer, Double> sigmaYMap = new HashMap<Integer, Double>();
            HashMap<Integer, Double> thetaMap = new HashMap<Integer, Double>();

            HashMap<Integer, Double> xEIMap = new HashMap<Integer, Double>();
            HashMap<Integer, Double> yEIMap = new HashMap<Integer, Double>();
            HashMap<Integer, Double> sigmaXEIMap = new HashMap<Integer, Double>();
            HashMap<Integer, Double> sigmaYEIMap = new HashMap<Integer, Double>();
            HashMap<Integer, Double> thetaEIMap = new HashMap<Integer, Double>();

            HashMap<Integer, double[]> contourXMap = new HashMap<Integer, double[]>();
            HashMap<Integer, double[]> contourYMap = new HashMap<Integer, double[]>();
            HashMap<Integer, double[]> simpleContourXMap = new HashMap<Integer, double[]>();
            HashMap<Integer, double[]> simpleContourYMap = new HashMap<Integer, double[]>();

            int totalSpikes = 0;
            for (int i = 0; i < neurons.length; i++) {
                totalSpikes += neuronFile.getSpikeCount(neurons[i]);

                if ((showGaussians || showDisplacements) && movieExists) {
                    xMap.put(neurons[i], paramsFile.getDoubleCell(neurons[i], "x0"));
                    yMap.put(neurons[i], paramsFile.getDoubleCell(neurons[i], "y0"));
                    sigmaXMap.put(neurons[i], paramsFile.getDoubleCell(neurons[i], "SigmaX"));
                    sigmaYMap.put(neurons[i], paramsFile.getDoubleCell(neurons[i], "SigmaY"));
                    thetaMap.put(neurons[i], paramsFile.getDoubleCell(neurons[i], "Theta"));
                }
                
                if (showEIGaussians || showDisplacements) {
                    xEIMap.put(neurons[i], paramsFile.getDoubleCell(neurons[i], "EIx0"));
                    yEIMap.put(neurons[i], paramsFile.getDoubleCell(neurons[i], "EIy0"));
                    sigmaXEIMap.put(neurons[i], paramsFile.getDoubleCell(neurons[i], "EISigmaX"));
                    sigmaYEIMap.put(neurons[i], paramsFile.getDoubleCell(neurons[i], "EISigmaY"));
                    thetaEIMap.put(neurons[i], paramsFile.getDoubleCell(neurons[i], "EITheta"));
                }
                
                if (showContours && movieExists) {
                    contourXMap.put(neurons[i], paramsFile.getArrayCell(neurons[i], "contourX"));
                    contourYMap.put(neurons[i], paramsFile.getArrayCell(neurons[i], "contourY"));
                    simpleContourXMap.put(neurons[i], paramsFile.getArrayCell(neurons[i], "simpleContourX"));
                    simpleContourYMap.put(neurons[i], paramsFile.getArrayCell(neurons[i], "simpleContourY"));
                }
            }
            
            PlotPanel p = new PlotPanel("mosaic");
            p.setAntiallias(true);
            
            // if centerX and centerY are not available from file, assume that they match the center of the stas.
            if ((Double.isNaN(centerX) || Double.isNaN(centerY)) && movieExists) {
                centerX = 0.0;
                centerY = 0.0;
                int count = 0;
                for (Integer id : xMap.keySet()) {
                    Double x0 = xMap.get(id);
                    Double y0 = yMap.get(id);
                    if (x0 != null && y0 != null && !Double.isNaN(x0.doubleValue()) && !Double.isNaN(y0.doubleValue())) {
                        centerX += (x0.doubleValue()*runParams.micronsPerStixelX+runParams.xOffset);
                        centerY +=  (y0.doubleValue()*runParams.micronsPerStixelY+runParams.yOffset);
                        count++;
                    }
                }
                if (count!=0) {
                    centerX/=count;
                    centerY/=count;
                } else {
                    centerX = runParams.micronsPerStixelX*runParams.width/2 + runParams.xOffset;
                    centerY = runParams.micronsPerStixelY*runParams.height/2 + runParams.yOffset;
                }

                p.addBackgroundText("Corners not available; Calibration approximated", PlotPanel.LEFT, PlotPanel.BOTTOM, Font.decode("Arial 10"), Color.black);
            }
            
            boolean arrowValid;
            double[] staCenter = new double[2], eiCenter = new double[2];
            double[] staSigma = new double[2], eiSigma = new double[2];
            MeanVarianceCalculator mvc = new MeanVarianceCalculator(MeanVarianceCalculator.UNBIASED);
            for (int i = 0; i < neurons.length; i++) {
                int id = neurons[i];
                Double x0 = xMap.get(id);
                Double y0 = yMap.get(id);
                Double sigmaX = sigmaXMap.get(id);
                Double sigmaY = sigmaYMap.get(id);
                Double theta = thetaMap.get(id);

                arrowValid = true;
                if (x0 != null && y0 != null && (showDisplacements || showGaussians || showLegend)) {
                    staCenter[0] = x0 + runParams.xOffset/runParams.micronsPerStixelX;
                    staCenter[1] = y0 + runParams.yOffset/runParams.micronsPerStixelY;
                    staSigma[0] = sigmaX;
                    staSigma[1] = sigmaY;
                    if (!Double.isNaN(staCenter[0]) && !Double.isNaN(staCenter[1])) {

                        mvc.add(Math.sqrt(sigmaX*sigmaY));
                        //Must use built in scaling functions for parametric ellipse.
                        //Otherwise, asymmetry between X and Y will bite you.
                        //major and minor axes do not point in the direction of x and y.
                        if (showGaussians) {
                            ParametricEllipse fit = new ParametricEllipse(staCenter[0], staCenter[1], staSigma[0], staSigma[1], -theta,
                                                        runParams.micronsPerStixelX/1000, runParams.micronsPerStixelY/1000); 
                            if (showIDs) fit.setLabel("" + id);
                            p.addData(fit, gaussianFitStyle);
                        }				
                        staCenter[0]*=runParams.micronsPerStixelX/1000;
                        staCenter[1]*=runParams.micronsPerStixelY/1000;
                    } else {
                        arrowValid = false;
                    }
                }

                Double x0EI = xEIMap.get(id);
                Double y0EI = yEIMap.get(id);
                Double sigmaXEI = sigmaXEIMap.get(id);
                Double sigmaYEI = sigmaYEIMap.get(id);

                Double thetaEI = thetaEIMap.get(id);
                if (thetaEI != null && ((flipX ? -1 : 1) * (flipY ? -1 : 1) == -1)) {
                    while(thetaEI > Math.PI) {thetaEI-=2*Math.PI;}
                    while(thetaEI < -Math.PI) {thetaEI+=2*Math.PI;}
                    thetaEI = -thetaEI;
                }

                if (x0EI != null && y0EI != null && (showDisplacements || showEIGaussians)) {

                    eiCenter[0] = (flipX ? -1 : 1) * x0EI.doubleValue();
                    eiCenter[1] = (flipY ? -1 : 1) * y0EI.doubleValue();

                    rotate2D(eiCenter, angle);

                    eiCenter[0] = (eiCenter[0]+centerX)/1000;
                    eiCenter[1] = (eiCenter[1]+centerY)/1000;

                    eiSigma[0] = sigmaXEI/1000.0;
                    eiSigma[1] = sigmaYEI/1000.0;

                    if (!Double.isNaN(eiCenter[0]) && !Double.isNaN(eiCenter[1])) {
                        if (showEIGaussians)
                            p.addData(new ParametricEllipse(eiCenter[0], eiCenter[1], eiSigma[0], eiSigma[1], -thetaEI+angle), gaussianFitStyleEI);
                    } else {
                        arrowValid = false;
                    }
                }

                if (arrowValid == true && showDisplacements)
                    p.addArrow(new Arrow(eiCenter[0], eiCenter[1], staCenter[0], staCenter[1]));
            }

            if (showContours && movieExists) {
                for(Integer id: contourXMap.keySet()) {
                    double[] xVals = contourXMap.get(id);
                    double[] yVals = contourYMap.get(id);
                    double[] xValsSimple = simpleContourXMap.get(id);
                    double[] yValsSimple = simpleContourYMap.get(id);

                    p.addDataPlotter(new PolygonPlotter());

                    if (xVals != null && xValsSimple != null) {
                        Polygon2D polygon = new Polygon2D(xVals, yVals, xVals.length);
                        Polygon2D simplePolygon = new Polygon2D(xValsSimple, yValsSimple,
                                xValsSimple.length);

                        Polygon2D[] polygons = polygon.divideIntoSeparatePolygons();

                        Polygon2DAdapter p2d;

                        p2d = new Polygon2DAdapter();
                        simplePolygon.transform(runParams.micronsPerStixelX/1000, 
                                runParams.micronsPerStixelY/1000, runParams.xOffset/1000, runParams.yOffset/1000);
                        p2d.setPolygon(simplePolygon);
                        p.addData(p2d,simpleContour);
                        for (int j = 0; j < polygons.length; j++) {

                            p2d = new Polygon2DAdapter();
                            polygons[j].transform(runParams.micronsPerStixelX/1000, 
                                    runParams.micronsPerStixelY/1000, runParams.xOffset/1000, runParams.yOffset/1000);
                            p2d.setPolygon(polygons[j]);
                            p.addData(p2d, fullContour);
                        }
                    }
                }
            }


            String className = InteractiveTree.pathToString(classPath).replace("/", " ").replace("All", "");

            p.addToLegend(viewer.experimentName + File.separator + viewer.datasetName);
            p.addToLegend(neurons.length + " " + className + " cells");
            if (neurons.length != 0 && movieExists) {
                double c = 2 * Math.sqrt(runParams.micronsPerStixelX*runParams.micronsPerStixelY);
                p.addToLegend(
                        "RF Diam: " +
                        StringUtil.format(c * mvc.getMean(), 1) + "\u00B1" +
                        StringUtil.format(c * mvc.getMeanVariance(), 1) + "\u03bcm");	
            }
            int kSpikesPerNeuron = (int) round(totalSpikes / neurons.length);
            p.addToLegend("Spikes/cell: " + kSpikesPerNeuron);
            p.setAxisVisible(showAxes);
            p.setLegendVisible(showLegend);
            p.setLabels("x (mm)", "y (mm)");
            if (movieExists) {
                p.setRange(runParams.xOffset/1000,
                        (runParams.micronsPerStixelX*runParams.width+runParams.xOffset)/1000,
                        runParams.yOffset/1000,
                        (runParams.micronsPerStixelY*runParams.height+runParams.yOffset)/1000);
            } else{
                double[] bounds = viewer.electrodeMap.getBounds();
                p.setRange(bounds[0]/1000 - .25, bounds[1]/1000 +.25, bounds[2]/1000- .25, bounds[3]/1000 + .25);
            }

            if (movieExists) {
                final double pixelWidth = runParams.micronsPerStixelX;
                final double pixelHeight = runParams.micronsPerStixelY;
                final double xOffset = runParams.xOffset;
                final double yOffset = runParams.yOffset;
                final TreePath finalClassPath = classPath;
                p.addSelectionAction(new SelectionAction("Print IDs") {
                    public void selectionPerformed(JComponent source, Selection selection) {
                        SelectionRegion r = selection.getSelection();

                        for (int i = 0; i < neurons.length; i++) {
                            //add the timecourse data
                            double x = (pixelWidth * paramsFile.getDoubleCell(neurons[i], "x0") + xOffset)/1000;
                            double y = (pixelHeight * paramsFile.getDoubleCell(neurons[i], "y0") + yOffset)/1000;
                            //  System.out.println(x + "   " + y);
                            if (r.contains(x, y)) {
                                System.out.println(neurons[i]);
                            }
                        }
                    }
                });

                p.addSelectionAction(new SelectionAction("New Class") {
                    public void selectionPerformed(JComponent source, Selection selection) {
                        SelectionRegion r = selection.getSelection();
                        IntegerList idList = new IntegerList();

                        for (int i = 0; i < neurons.length; i++) {
                            //add the timecourse data
                            double x = (pixelWidth * paramsFile.getDoubleCell(neurons[i], "x0") + xOffset)/1000;
                            double y = (pixelHeight * paramsFile.getDoubleCell(neurons[i], "y0") + yOffset)/1000;
                            if (r.contains(x, y)) {
                                idList.add(neurons[i]);
                            }
                        }

                        String className = "nc" + viewer.newClassIndex++;
                        viewer.classificationHappened(
                                idList.toArray(), finalClassPath.pathByAddingChild(new
                                        DefaultMutableTreeNode(className, true)));
                    }
                });
            }


            return p;  
        } catch (Exception ex) {
            ex.printStackTrace();
            return new JLabel("Could not make mosaic plot.");

        }
    }

    public void rotate2D(double[] point, double angle) {
        double newX = point[0]*Math.cos(angle)-point[1]*Math.sin(angle);
        double newY = point[0]*Math.sin(angle)+point[1]*Math.cos(angle);
        point[0] = newX;
        point[1] = newY;
    }

}
