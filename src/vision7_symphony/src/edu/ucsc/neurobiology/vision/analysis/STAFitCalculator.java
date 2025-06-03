package edu.ucsc.neurobiology.vision.analysis;

import edu.ucsc.neurobiology.vision.math.fitting.Gaussian2DFunction;
import edu.ucsc.neurobiology.vision.plot.Polygon2D;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class STAFitCalculator implements ParametersCalculator {

    public boolean fitContours;
    public double colorToFit;
    public double contourCutValue;
    public int downsampleFactor;
    
    private enum ColorToFit {RED, GREEN, BLUE};

    public void init(String rootPath, String mainFilePath) {}

    public String getName() {
        return "Gaussian Fit";
    }


    public String[][] getParameterTypes() {
        return new String[][] { 
                {"x0", "Double"}
                , {"y0", "Double"}
                , {"SigmaX", "Double"}
                , {"SigmaY", "Double"}
                , {"Theta", "Double"}
                , {"gAmp", "Double"}

                , {"contourX", "DoubleArray"}
                , {"contourY", "DoubleArray"}
                , {"contourArea", "Double"}
                , {"simpleContourX", "DoubleArray"}
                , {"simpleContourY", "DoubleArray"}
                , {"simpleContourArea", "Double"}

        };
    }


    public Object[] getParameterValues(ParameterCalculatorContext c) {
        Object[] gaussianParams;
        Gaussian2DFunction gFit = null;
        if (c.currentSTA != null) {
            gFit = c.currentSTA.fit(downsampleFactor);
        }

        if (gFit != null) {
            gFit.scaleUp(downsampleFactor);
            gaussianParams = new Object[] {
                    new Double(gFit.getX0()),
                    new Double(gFit.getY0()),
                    new Double(gFit.getSigmaX()),
                    new Double(gFit.getSigmaY()),
                    new Double(gFit.getTheta()),
                    new Double(gFit.getAmplitude())
            };
        } else {
            gaussianParams = new Object[] {
                    new Double(Double.NaN),
                    new Double(Double.NaN),
                    new Double(Double.NaN),
                    new Double(Double.NaN),
                    new Double(Double.NaN),
                    new Double(Double.NaN),
            };
        }

        Polygon2D poly;
        Polygon2D[] polys;
        Object[] contourParams;
        if (fitContours && c.currentSTA != null && c.currentSTA.getMainFrameIndex() != -1) {
            poly = c.currentSTA.getContour(contourCutValue, (int) colorToFit, true);
            polys = poly.divideIntoSeparatePolygons();
            double totalPolyArea = 0.0;
            for (int j = 0; j < polys.length; j++) {
                totalPolyArea += polys[j].getArea();
            }
            Polygon2D simple = poly.simplifyPolygon();
            double simpleArea = simple.getArea();

            contourParams = new Object[] {
                    poly.xPoints,
                    poly.yPoints,
                    new Double(totalPolyArea),
                    simple.xPoints,
                    simple.yPoints,
                    new Double(simpleArea)
            };
        } else {
            contourParams = new Object[] {
                    new double[0],
                    new double[0],
                    new Double(Double.NaN),
                    new double[0],
                    new double[0],
                    new Double(Double.NaN)
            };
        }
        Object[] params = new Object[gaussianParams.length + contourParams.length];
        for (int i = 0; i < gaussianParams.length; i++) {
            params[i] = gaussianParams[i];
        }
        for (int i = 0; i < contourParams.length; i++) {
            params[gaussianParams.length + i] = contourParams[i];
        }
        return params;
    }
}
