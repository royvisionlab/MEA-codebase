package edu.ucsc.neurobiology.vision.analysis;

import edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMap;
import edu.ucsc.neurobiology.vision.math.FitFailedException;
import edu.ucsc.neurobiology.vision.math.fitting.Fitter;
import edu.ucsc.neurobiology.vision.math.fitting.Gaussian2DFunction;
import edu.ucsc.neurobiology.vision.plot.PlotUtil;


/**
 *
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class EIFitCalculator
implements ParametersCalculator {

    private final static boolean DEBUG = false;


    public void init(String rootPath, String mainFilePath) {}


    public String getName() {
        return "EI Fit";
    }


    public String[][] getParameterTypes() {
        return new String[][] { 
                {"EIx0", "Double"}//microns
                , {"EIy0", "Double"}//microns
                , {"EISigmaX", "Double"}//microns
                , {"EISigmaY", "Double"}//microns
                , {"EITheta", "Double"}//radians
                , {"EIgAmp", "Double"}

        };
    }


    public Object[] getParameterValues(ParameterCalculatorContext c) {
        Gaussian2DFunction gFit = null;
        if (c.currentEI != null) {
            gFit = fitEI(c.currentEI, c.map);
        }

        Object[] gaussianParams;
        if (gFit != null) {
            return gaussianParams = new Object[] {
                    new Double(gFit.getX0()),
                    new Double(gFit.getY0()),
                    new Double(gFit.getSigmaX()),
                    new Double(gFit.getSigmaY()),
                    new Double(gFit.getTheta()),
                    new Double(gFit.getAmplitude())
            };
        } else {
            return gaussianParams = new Object[] {
                    new Double(Double.NaN),
                    new Double(Double.NaN),
                    new Double(Double.NaN),
                    new Double(Double.NaN),
                    new Double(Double.NaN),
                    new Double(Double.NaN),
            };
        }
    }
    
    
    /**
     * Fits a Gaussian to an EI
     * @param float[image or error][electrode][time point] electrophysiological image
     * @param ElectrodeMap
     * @return Gaussian2DFunction
     */
    public Gaussian2DFunction fitEI(float[][][] ei, ElectrodeMap map) {


        //Convert EI to fitter friendly format.
        int nElectrodes = ei[0].length;	
        double[][] x = new double[2][nElectrodes];
        double[] y = new double[nElectrodes];
        double[] sigma = new double[nElectrodes];


        //for each electrode, find max
        int maxElectrode = -1;
        double maxOfMaxes = -1;
        for(int electrode=0; electrode<nElectrodes; electrode++) {
            double max = -1;
//			double maxTime = -1;
            double errorAtMax = -1;
            for(int time=0; time<ei[0][electrode].length; time++) {
                if(Math.abs(ei[0][electrode][time]) > max) {
                    max = Math.abs(ei[0][electrode][time]);
                    errorAtMax = ei[1][electrode][time];
//					maxTime = time;
                }	
            }

            x[0][electrode]=map.getXPosition(electrode);
            x[1][electrode]=map.getYPosition(electrode);
            y[electrode]  = max;
            sigma[electrode] = errorAtMax;

            if(max > maxOfMaxes) {
                maxOfMaxes = max;
                maxElectrode = electrode;
            }

        }




        Fitter fitter = new Fitter();
        fitter.setLastChiSquaredVariation(1e-3);
        fitter.setMaxIterations(1000);
        double[] params = new double[7];
        Gaussian2DFunction gBest = null;
        double minChi2 = Double.POSITIVE_INFINITY;
        
        //If maxElectrode is still -1, fitting will fail. Return -1
        if(maxElectrode == -1) {
         return null;
        }

        //Fit is highly unstable.  Try a range of sigmas to make getting a good fit likely.
        for (double i = 10; i <= 100; i += 10) {
            for (double j = 10; j <= 100; j += 10) {

                params[0] = maxOfMaxes; // A
                params[1] = map.getXPosition(maxElectrode);
                params[2] = map.getYPosition(maxElectrode); // y0
                params[3] = i; // sigX
                params[4] = j; // sigY
                params[5] = 0; // Theta, radians
                params[6] = 0; //B


                Gaussian2DFunction g2d = new Gaussian2DFunction();
                g2d.setParameters(params);
                try {
                    fitter.fit(g2d, x, y, sigma, nElectrodes);
                    if (DEBUG) {
                        System.out.println(g2d);
                    }
                    double chi2 = g2d.getChiSquared();
                    if (chi2 < minChi2) {
                        minChi2 = chi2;
                        gBest = g2d;
                    }
                } catch (FitFailedException e) {
                    if (DEBUG) {
                        System.out.println(e.getMessage() + "\n");
                    }
                }
            }
        }

        // validity cuts
        if (gBest == null) {
            return null;
        } else {
            double[] bounds = map.getBounds();
            gBest.normalizeVariables();
            if (gBest.getX0() < bounds[0] - 100||
                gBest.getX0() > bounds[1] + 100||
                gBest.getY0() < bounds[2] - 100||
                gBest.getY0() > bounds[3] + 100 ||
                Math.abs(gBest.getSigmaX()) > (bounds[1] - bounds[0]) ||
                Math.abs(gBest.getSigmaY()) > bounds[3] - bounds[2]) {

                return null;
            }
        }
        
        
        return gBest;
    }



    public static void main(String[] args) {
        int width = 64;
        int height = 64;
        int nPoints = width*height;
        double[][] x = new double[2][nPoints];

        double[] y = new double[nPoints];


        double[] sigma = new double[nPoints];

        int n=0;
//		double average = 0;
        for(int i=0; i<width; i++) {
            for(int j=0; j<height; j++) {
                x[0][n] = i+.5;
                x[1][n] = height-j-1+.5;
                y[n] = Math.pow(Math.E, -((x[0][n]-width/2)*(x[0][n]-width/2)+(x[1][n]-height/2)*(x[1][n])-height/2)/100)/16;
                sigma[n] = Math.sqrt(y[n]);
//				average+=y[n];
                n++;
            }
        }
//		average /= y.length;
        for (int i = 0; i < y.length; i++) {
            //           y[i] -= average;
        }
        PlotUtil.showArray("", y);


        Fitter fitter = new Fitter();
        fitter.setLastChiSquaredVariation(1e-3);
        fitter.setMaxIterations(1000);
        double[] params = new double[7];

        Gaussian2DFunction gBest = null;
        double minChi2 = Double.POSITIVE_INFINITY;

        for(double i=1; i<=3; i+=0.5) {
            for(double j=1; j<=3; j+=.5) {
                params[0] = 1;//maxOfMaxes; // A
                params[1] = 32;//map.getXPosition(maxElectrode);
                params[2] = 32;//map.getYPosition(maxElectrode); // y0
                params[3] = i; // sigX
                params[4] = j; // sigY
                params[5] = 0; // Theta, radians
                params[6] = 0; //B


                Gaussian2DFunction g2d = new Gaussian2DFunction();
                g2d.setParameters(params);
                try {
                    fitter.fit(g2d, x, y, sigma, nPoints);
                    if (DEBUG) {
                        System.out.println(g2d);
                    }
                    double chi2 = g2d.getChiSquared();
                    if (chi2 < minChi2) {
                        minChi2 = chi2;
                        gBest = g2d;
                    }
                } catch (FitFailedException e) {
                    if (DEBUG) {
                        System.out.println(e.getMessage() + "\n");
                    }
                }
            }


        }

        System.out.println(gBest);
    }
}
