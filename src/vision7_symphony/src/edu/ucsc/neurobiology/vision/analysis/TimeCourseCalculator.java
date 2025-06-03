package edu.ucsc.neurobiology.vision.analysis;

import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.math.fitting.*;
import edu.ucsc.neurobiology.vision.plot.*;

import java.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class TimeCourseCalculator implements ParametersCalculator {
    // These normally get assigned from config.xml immediately after this is created.
    public double significance = 3;
    public int nTemporalSubfilters = 2;
    public double latestZero = -40;
    public double stimulusVariance = 0.0;

    public void init(String rootPath, String mainFilePath) {}

    public String getName() {
        return "Time Course";
    }

    public String[][] getParameterTypes() {
        return new String[][] { {"RedTimeCourse", "DoubleArray"}
        , {"GreenTimeCourse", "DoubleArray"}
        , {"BlueTimeCourse", "DoubleArray"}

        , {"t1", "Double"}
        , {"t2", "Double"}
        , {"t3", "Double"}

        , {"a1", "Double"}
        , {"a2", "Double"}
        , {"a3", "Double"}

        , {"n1", "Double"}
        , {"n2", "Double"}
        , {"n3", "Double"}
        , {"tOffset", "Double"}

        , {"dot", "Double"}
        , {"dot2", "Double"}
        , {"srm", "Double"}
        , {"rl", "Double"}
        , {"amp1", "Double"}
        , {"amp2", "Double"}
        , {"amp3", "Double"}

        , {"blueness", "Double"}

        , {"RedVTimeCourse", "DoubleArray"}
        , {"GreenVTimeCourse", "DoubleArray"}
        , {"BlueVTimeCourse", "DoubleArray"}

        , {"t1V", "Double"}
        , {"t2V", "Double"}
        , {"t3V", "Double"}

        , {"a1V", "Double"}
        , {"a2V", "Double"}
        , {"a3V", "Double"}

        , {"n1V", "Double"}
        , {"n2V", "Double"}
        , {"n3V", "Double"}
        , {"tVOffset", "Double"}

        , {"dotV", "Double"}
        , {"dot2V", "Double"}
        , {"srmV", "Double"}
        , {"rlV", "Double"}
        , {"amp1V", "Double"}
        , {"amp2V", "Double"}
        , {"amp3V", "Double"}

        , {"bluenessV", "Double"}
        };
    }

    public Object[] getParameterValues(ParameterCalculatorContext c) {
        double[][] timeFilters = null, timeFiltersError = null, varianceTimeFilters = null;
        double refreshTime = 0.0;
        if (c.currentSTA != null) {
            timeFilters = c.currentSTA.getTimeFilters(significance);
            timeFiltersError = c.currentSTA.getTimeFiltersError();
            refreshTime = c.currentSTA.getRefreshTime();
            
            if (c.stv != null)
                varianceTimeFilters = c.currentSTA.getComparisonTimeFilters(significance, c.stv);
        }
        
        Object[] timeCourseParams  = getTimeCourseParams(timeFilters,         timeFiltersError, refreshTime, 0.0);
        Object[] timeCourseEParams = getTimeCourseParams(varianceTimeFilters, timeFiltersError, refreshTime, stimulusVariance);

        Object[] toReturn = new Object[timeCourseParams.length + timeCourseEParams.length];
        for (int i = 0; i < timeCourseParams.length; i++) {
            toReturn[i] = timeCourseParams[i]; 
        }
        for (int i = timeCourseParams.length; i < toReturn.length; i++) {
            toReturn[i] = timeCourseEParams[i-timeCourseParams.length];
        }

        return toReturn;
    }

    public Object[] getTimeCourseParams(double[][] timeFilters, double[][] timeFiltersError, double refreshTime, double offset) {
        double t1 = Double.NaN;
        double t2 = Double.NaN;
        double t3 = Double.NaN;
        double a1 = Double.NaN;
        double a2 = Double.NaN;
        double a3 = Double.NaN;
        double n1 = Double.NaN;
        double n2 = Double.NaN;
        double n3 = Double.NaN;
        double dot = Double.NaN;
        double dot2 = Double.NaN;
        double srm = Double.NaN;
        double rl = Double.NaN;
        double amp1 = offset;
        double amp2 = offset;
        double amp3 = offset;


        if (timeFilters != null) {
            if (nTemporalSubfilters != 0) {
                int[] n1n2n3 = new int[3];
                try {
                    STATimeFunction1 f = null;
                    if (nTemporalSubfilters == 3) {
                        f = fitSTATimeFilter(timeFilters[1],
                                timeFiltersError[1], false, 0,
                                30, 0, 16, 8, 16, n1n2n3,
                                refreshTime, offset);

                        // Sort based on t values.
                        if (f.getT1() < f.getT2()) {
                            f = new STATimeFunction1(
                                    f.getA2(), f.getT2(), f.getA1(), f.getT1(), f.getA3(),
                                    f.getT3(), f.getN2(), f.getN1(), f.getN3(), f.offset);
                        }
                        if (f.getT2() < f.getT3()) {
                            f = new STATimeFunction1(
                                    f.getA1(), f.getT1(), f.getA3(), f.getT3(), f.getA2(),
                                    f.getT2(), f.getN1(), f.getN3(), f.getN2(), f.offset);
                        }
                        if (f.getT1() < f.getT2()) {
                            f = new STATimeFunction1(
                                    f.getA2(), f.getT2(), f.getA1(), f.getT1(), f.getA3(),
                                    f.getT3(), f.getN2(), f.getN1(), f.getN3(), f.offset);
                        }
                    } else if (nTemporalSubfilters == 2) {
                        f = fitSTATimeFilter(timeFilters[1],
                                timeFiltersError[1], false, 1,
                                30, 1, 20, n1n2n3, refreshTime, offset);

                        //sort based on t values
                        if (f.getT1() < f.getT2()) {
                            f = new STATimeFunction1(
                                    f.getA2(), f.getT2(), f.getA1(), f.getT1(), n1n2n3[1],
                                    n1n2n3[0], f.offset);
                        }
                    } else {
                        f = fitSTATimeFilter(timeFilters[1],
                                timeFiltersError[1], false, 0,
                                30, n1n2n3, refreshTime, offset);
                    }

                    t1 = f.getT1();
                    t2 = f.getT2();
                    t3 = f.getT3();
                    a1 = f.getA1();
                    a2 = f.getA2();
                    a3 = f.getA3();
                    n1 = f.getN1();
                    n2 = f.getN2();
                    n3 = f.getN3();

                    f.offset = 0; //to make parameter calculations more correct
                    // calculate extra TC params
                    double x1 = -refreshTime * timeFilters[1].length;

                    double i1 = Integrator.gauss8(f, x1, 0, 0.01);
                    double i2 = Integrator.gauss8(new AbsAdapterFunction(f), x1, 0, 0.01);
                    dot = 1 - Math.abs(i1) / i2;

                    double onCell = f.getValueAt(latestZero);
                    onCell = onCell / Math.abs(onCell);
                    dot2 = 1 -  onCell * i1 / i2;
                    srm = i1/i2;

                    rl = findResponceLatency(f, x1, 0, 10000);
                    f.offset = offset;

                    double[] peakAmps = findPeakAmps(f, x1, latestZero, .01);

                    
                    for (int i = peakAmps.length-1; i >= 0; i--) {
                        if (i == peakAmps.length-1) {
                            amp1 = peakAmps[i];
                        } else if (i == peakAmps.length-2) {
                            amp2 = peakAmps[i];
                        } else if(i == peakAmps.length-3) {
                            amp3 = peakAmps[i];
                        } else {
                            break;
                        }
                    }


                } catch (Exception e) {
//					System.out.println("Timecourse fit failed on neuron " + c.id + ".");
                }
            }
            
            double blueness;
            if (timeFilters != null) {
                double greenSensitivity = 0.0;
                double blueSensitivity = 0.0;
                for (int i = timeFilters[1].length / 2; i < timeFilters[1].length; i++) {
                    greenSensitivity += (Math.abs(timeFilters[1][i])-offset);
                    blueSensitivity  += (Math.abs(timeFilters[2][i])-offset);
                }
                blueness = blueSensitivity / greenSensitivity;
            } else {
                blueness = 0.0;
            }

                return new Object[] {
                    (timeFilters != null) ? timeFilters[0] : null,
                    (timeFilters != null) ? timeFilters[1] : null,
                    (timeFilters != null) ? timeFilters[2] : null,

                    new Double(t1),
                    new Double(t2),
                    new Double(t3),
                    new Double(a1),
                    new Double(a2),
                    new Double(a3),
                    new Double(n1),
                    new Double(n2),
                    new Double(n3),
                    new Double(offset),

                    new Double(dot),
                    new Double(dot2),
                    new Double(srm),
                    new Double(rl),
                    new Double(amp1),
                    new Double(amp2),
                    new Double(amp3),

                    new Double(blueness)

                };
        } else {
            return new Object[] {
                null, null, null,
                new Double(Double.NaN), new Double(Double.NaN), new Double(Double.NaN), new Double(Double.NaN),
                new Double(Double.NaN), new Double(Double.NaN), new Double(Double.NaN), new Double(Double.NaN),
                new Double(Double.NaN), new Double(Double.NaN), new Double(Double.NaN), new Double(Double.NaN),
                new Double(Double.NaN), new Double(Double.NaN), new Double(Double.NaN), new Double(Double.NaN),
                new Double(Double.NaN), new Double(Double.NaN)};
        }
    }





    public Object[] getParameterValuesOld(ParameterCalculatorContext c) {
        // get the time courses
        double[][] timeFilters = null, timeFiltersError = null;
        double[][] varianceTimeFilters = null;
        
        if (c.currentSTA != null) {
            timeFilters = c.currentSTA.getTimeFilters(significance);
            timeFiltersError = c.currentSTA.getTimeFiltersError();
            varianceTimeFilters =
                (c.stv != null) ?
                        c.currentSTA.getComparisonTimeFilters(significance, c.stv) : null;
        }

        double t1 = Double.NaN;
        double t2 = Double.NaN;
        double t3 = Double.NaN;
        double a1 = Double.NaN;
        double a2 = Double.NaN;
        double a3 = Double.NaN;
        double n1 = Double.NaN;
        double n2 = Double.NaN;
        double n3 = Double.NaN;
        double dot = Double.NaN;
        double dot2 = Double.NaN;
        double srm = Double.NaN;
        double rl = Double.NaN;



        if (nTemporalSubfilters != 0) {
            int[] n1n2n3 = new int[3];
            try {
                STATimeFunction1 f = null;
                if (nTemporalSubfilters == 3) {
                    f = fitSTATimeFilter(timeFilters[1],
                            timeFiltersError[1], false, 0,
                            30, 0, 16, 8, 16, n1n2n3,
                            c.currentSTA.getRefreshTime(), 0.0);

                    //Sort based on t values.
                    if (f.getT1() < f.getT2()) {
                        f = new STATimeFunction1(
                                f.getA2(), f.getT2(), f.getA1(), f.getT1(), f.getA3(),
                                f.getT3(), f.getN2(), f.getN1(), f.getN3(), f.offset);
                    }
                    if (f.getT2() < f.getT3()) {
                        f = new STATimeFunction1(
                                f.getA1(), f.getT1(), f.getA3(), f.getT3(), f.getA2(),
                                f.getT2(), f.getN1(), f.getN3(), f.getN2(), f.offset);
                    }
                    if (f.getT1() < f.getT2()) {
                        f = new STATimeFunction1(
                                f.getA2(), f.getT2(), f.getA1(), f.getT1(), f.getA3(),
                                f.getT3(), f.getN2(), f.getN1(), f.getN3(), f.offset);
                    }
                } else if (nTemporalSubfilters == 2) {
                    f = fitSTATimeFilter(timeFilters[1],
                            timeFiltersError[1], false, 1,
                            30, 1, 20, n1n2n3, c.currentSTA.getRefreshTime(), 0.0);

                    //sort based on t values
                    if (f.getT1() < f.getT2()) {
                        f = new STATimeFunction1(
                                f.getA2(), f.getT2(), f.getA1(), f.getT1(), n1n2n3[1],
                                n1n2n3[0], f.offset);
                    }
                } else {
                    f = fitSTATimeFilter(timeFilters[1],
                            timeFiltersError[1], false, 0,
                            30, n1n2n3, c.currentSTA.getRefreshTime(), 0.0);
                }

                t1 = f.getT1();
                t2 = f.getT2();
                t3 = f.getT3();
                a1 = f.getA1();
                a2 = f.getA2();
                a3 = f.getA3();
                n1 = f.getN1();
                n2 = f.getN2();
                n3 = f.getN3();

                // calculate extra TC params
                double x1 = -c.currentSTA.getRefreshTime() * timeFilters[1].length;

                double i1 = Integrator.gauss8(f, x1, 0, 0.01);
                double i2 = Integrator.gauss8(new AbsAdapterFunction(f), x1, 0, 0.01);
                dot = 1 - Math.abs(i1) / i2;

                double onCell = f.getValueAt(latestZero);
                onCell = onCell / Math.abs(onCell);
                dot2 = 1 -  onCell * i1 / i2;
                srm = i1/i2;

                rl = findResponceLatency(f, x1, 0, 10000);


                // calculate the amplitudes

            } catch (Exception e) {
//				System.out.println("Timecourse fit failed on neuron " + c.id + ".");
            }
        }
        double blueness;
        if (timeFilters != null) {
            double greenSensitivity = 0.0;
            double blueSensitivity = 0.0;
            for (int i = timeFilters[1].length / 2; i < timeFilters[1].length;
            i++) {
                greenSensitivity += Math.abs(timeFilters[1][i]);
                blueSensitivity += Math.abs(timeFilters[2][i]);

            }
            blueness = blueSensitivity / greenSensitivity;
        } else {
            blueness = 0.0;
        }

        return new Object[] {
                (timeFilters != null) ? timeFilters[0] : null,
                (timeFilters != null) ? timeFilters[1] : null,
                (timeFilters != null) ? timeFilters[2] : null,

                (varianceTimeFilters != null) ? varianceTimeFilters[0] : null,
                (varianceTimeFilters != null) ? varianceTimeFilters[1] : null,
                (varianceTimeFilters != null) ? varianceTimeFilters[2] : null,

                new Double(t1),
                new Double(t2),
                new Double(t3),
                new Double(a1),
                new Double(a2),
                new Double(a3),
                new Double(n1),
                new Double(n2),
                new Double(n3),

                new Double(dot),
                new Double(dot2),
                new Double(srm),
                new Double(rl),

                new Double(blueness)
        };
    }


    public static double findResponceLatency(FunctionData f, double a, double b,
            int nSteps) throws
            CannotEvaluateException {

        double dx = (b - a) / nSteps;
        double s = 0;
        double x = b - dx / 2;
        double rl = 0.0;
        double maxS = 0.0;

        for (int j = 0; j < nSteps; j++, x -= dx) {
            s += f.getValueAt(x);
            if (Math.abs(s) > maxS) {
                rl = x;
                maxS = Math.abs(s);
            }
        }

        return rl;
    }


    //Fit for a single filter
    public static STATimeFunction1 fitSTATimeFilter(double[] timeCourse,
            double[] timeCourseError, boolean isSTV, int I1, int I2, int[] n1n2n3,
            double refreshTime, double offset) throws Exception {
        return fitSTATimeFilter(timeCourse, timeCourseError, isSTV, I1, I2, 0, 0, 0, 0, 1,
                n1n2n3, refreshTime, offset);
    }


//	Fit for a double filter
    public static STATimeFunction1 fitSTATimeFilter(double[] timeCourse,
            double[] timeCourseError, boolean isSTV, int I1, int I2, int J1, int J2,
            int[] n1n2n3, double refreshTime, double offset) throws Exception {
        return fitSTATimeFilter(timeCourse, timeCourseError, isSTV, I1, I2, J1, J2, 0, 0,
                2, n1n2n3, refreshTime, offset);
    }


//	Fit for a triple filter
    public static STATimeFunction1 fitSTATimeFilter(double[] timeCourse,
            double[] timeCourseError, boolean isSTV, int I1, int I2, int J1, int J2, int K1,
            int K2, int[] n1n2n3, double refreshTime, double offset) throws Exception {
        return fitSTATimeFilter(timeCourse, timeCourseError, isSTV, I1, I2, J1, J2, K1,
                K2, 3, n1n2n3, refreshTime, offset);
    }


//	Leave private.  Use public functions, above.
    private static STATimeFunction1 fitSTATimeFilter(double[] timeCourse,
            double[] timeCourseError, boolean isSTV, int I1, int I2, int J1,
            int J2, int K1, int K2, int nFilters, int[] n1n2n3, double refreshTime, double offset) throws
            Exception {

        // extract the time-courses
        ScatterPlot sp = new ScatterPlot("");
        double[][] x = new double[1][timeCourse.length];

        //Code for STV
//		tCourse = sta.getComparisonTimeFilters(5.0, stv)[1];
//		MathUtil.sub(tCourse, MathUtil.mean(tCourse));
        //      MathUtil.multiply(timeCourseError,.4);
        ///      

        for (int f = 0; f < timeCourse.length; f++) {

            sp.add( - (timeCourse.length - f - 1) * refreshTime, timeCourse[f],
                    timeCourseError[f]);
            x[0][f] = - (timeCourse.length - f - 1) * refreshTime;
        }

        int max = MathUtil.maxAbsIndex(timeCourse);

        STATimeFunction1 bestF = null;
        double minChi2 = Double.MAX_VALUE;
        int _n1 = -1, _n2 = -1, _n3 = -1;

        for (int n1 = I1; n1 <= I2; n1++) {
            for (int n2 = J1; n2 <= J2; n2++) {
                for (int n3 = K1; n3 <= K2; n3++) {
                    try {
                        STATimeFunction1 t;
                        if (nFilters == 3) {
                            t = new STATimeFunction1(
                                    timeCourse[max], -40, -timeCourse[max] / 2, -70,
                                    timeCourse[max] / 4, -100,
                                    n1, n2, n3, offset);
                        } else if (nFilters == 2) {
                            t = new STATimeFunction1(
                                    timeCourse[max], -40, -timeCourse[max] / 2, -70,
                                    n1, n2, offset);

                        } else {
                            t = new STATimeFunction1(
                                    timeCourse[max], -40,
                                    n1, offset);

                        }

                        Fitter f = new Fitter();
                        f.fit(t, x, timeCourse, timeCourseError);

//						STATimeFunction1 t = fitSTA(sta, n1, n2, false);
                        double chi2 = t.getChiSquared();
                        if (chi2 < minChi2) {
                            minChi2 = chi2;
                            bestF = t;
                            _n1 = n1;
                            _n2 = n2;
                            _n3 = n3;
                        }
                    } catch (FitFailedException ex) {
//						System.err.println("Failed " + n1 + ", " + n2 + "," + n3);
                    }
                }
            }
        }

        if (n1n2n3 != null) {
            n1n2n3[0] = _n1;
            n1n2n3[1] = _n2;
            n1n2n3[2] = _n3;
        }

        return bestF;
    }

    //finds the amplitudes of peaks to accuracy dx.  Assumes that the function is smooth (that
    //is, it does not check for false positives).
    public static double[] findPeakAmps(FunctionData f, double x1, double x2, double dx) throws 
    CannotEvaluateException {
        double deriv = Double.NaN;
        double oldDeriv = (f.getValueAt(x1+dx) - f.getValueAt(x1))/dx;
        ArrayList<Double> peaksList = new ArrayList<Double>();
        for(double x=x1+dx; x < x2-dx; x+=dx) {
            deriv = (f.getValueAt(x+dx) - f.getValueAt(x-dx))/(2*dx);
            //if peak found
            if(deriv*oldDeriv <= 0.0) {
                peaksList.add(new Double(f.getValueAt(x)));
            }
            oldDeriv = deriv;
        }

        double[] peaks = new double[peaksList.size()];
        for(int i=0; i< peaksList.size(); i++) {
            peaks[i] = peaksList.get(i).doubleValue();
        }
        return peaks;
    }

    /*written by Matthew Grivich, returns the max(y-offset) between x1 and x2 on f to
    the accuracy of dx*.  Preserves the sign of y-offset.*/
    public static double findAbsMax(FunctionData f, double x1, double x2, double dx, double offset) throws
    CannotEvaluateException {
        double max = 0;
        for(double x=x1; x<x2; x+=dx) {
            double y = f.getValueAt(x)-offset;
            if(Math.abs(y)> Math.abs(max)) {
                max = y ;
            }
        }
        return max;
    }


}
