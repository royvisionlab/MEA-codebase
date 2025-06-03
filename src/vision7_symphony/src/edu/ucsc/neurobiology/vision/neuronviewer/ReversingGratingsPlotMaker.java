package edu.ucsc.neurobiology.vision.neuronviewer;

import java.awt.Component;
import java.io.IOException;
import java.util.HashMap;

import javax.swing.tree.TreePath;

import edu.ucsc.neurobiology.vision.Config;
import edu.ucsc.neurobiology.vision.analysis.ReversingGratings;
import edu.ucsc.neurobiology.vision.io.ParametersFile;
import edu.ucsc.neurobiology.vision.math.MathUtil;
import edu.ucsc.neurobiology.vision.math.MeanVarianceCalculator;
import edu.ucsc.neurobiology.vision.math.fitting.DOG1DFourierFunction;
import edu.ucsc.neurobiology.vision.math.fitting.Fitter;
import edu.ucsc.neurobiology.vision.plot.AxisType;
import edu.ucsc.neurobiology.vision.plot.PlotPanel;
import edu.ucsc.neurobiology.vision.plot.ScatterPlot;
import edu.ucsc.neurobiology.vision.util.IntegerList;
import edu.ucsc.neurobiology.vision.util.StringUtil;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @owner Matthew Grivich, UCSC
 */
public class ReversingGratingsPlotMaker
extends PlotMaker {

    public boolean showAverage = false;
    public boolean showF1 = true;
    public boolean showF2 = true;
    public boolean normalize = true;
    public double y1 = 1;
    public double y2 = 100;
    public AxisType xAxisType = AxisType.LOG10;
    public AxisType yAxisType = AxisType.LOG10;
    public enum PlotFunction {
        OVERLAYED, AVERAGE, FIT} ;
        public PlotFunction plotFunction = PlotFunction.OVERLAYED;


        public ReversingGratingsPlotMaker() {
            super("Reversing Gratings", CLASS_AND_NEURON_PLOT);
        }


        public void initialize(NeuronViewer viewer) {
            super.initialize(viewer);
        }


        public Component makePlot(IntegerList list, int plotType, TreePath classPath) {
            if (plotType == CLASS_PLOT && plotFunction.equals(PlotFunction.AVERAGE) &&
                    list.size() > 0) {
                try {
                    PlotPanel p = getAverageHarmonicsPanel(list.toArray(), normalize, paramsFile);
                    p.setYRange(y1, y2);
                    p.getXAxis().setAxisType(xAxisType);
                    p.getYAxis().setAxisType(yAxisType);
                    return p;
                } catch (IOException ex) {
                    return null;
                }
            }  else {
                PlotPanel p = getHarmonicsPanel(list.toArray(), paramsFile, showF1,
                        showF2, normalize, viewer.configuration);
                p.setYRange(y1, y2);
                p.getXAxis().setAxisType(xAxisType);
                p.getYAxis().setAxisType(yAxisType);
                return p;
            }
        }


        public static PlotPanel getHarmonicsPanel(int[] id, ParametersFile paramsFile,
                boolean showF1, boolean showF2, boolean normalize, Config config) {
            if (id.length == 0) {
                return new PlotPanel("F1F2");
            }


            double[] freqs = paramsFile.getArrayCell(id[0], "reversingFrequencies");
            PlotPanel p = new PlotPanel("F1F2");

            // calculate overall scaling factor
            double overallsum = 0;
            int nn = 0;

            // find the number of temporal frequencies
            int nTemporalPeriods = -1;
            for (int tPeriod = 0; tPeriod < 3; tPeriod++) {
                String f1Name = "T" + (tPeriod + 1) + "reversingF1";
                String f2Name = "T" + (tPeriod + 1) + "reversingF2";
                if (paramsFile.hasParameter(f1Name) && paramsFile.hasParameter(f2Name)) {
                    nTemporalPeriods = tPeriod + 1;
                }
            }

            // make the F graphs

            for (int i = 0; i < id.length; i++) {
                double norm = 0.0;

                for (int tPeriod = 0; tPeriod < nTemporalPeriods; tPeriod++) {
                    double[] _f1 = paramsFile.getArrayCell(id[
                                                              i], "T" + (tPeriod + 1) + "reversingF1");
                    double[] _f2 = paramsFile.getArrayCell(id[
                                                              i], "T" + (tPeriod + 1) + "reversingF2");
                    if(_f1 != null)
                        norm+= MathUtil.sumSquare(_f1);
                    if(_f2 != null)
                        norm+=MathUtil.sumSquare(_f2);
                }
                norm = Math.sqrt(norm);

                for (int tPeriod = 0; tPeriod < nTemporalPeriods; tPeriod++) {
                    double[] _f1 = paramsFile.getArrayCell(id[
                                                              i], "T" + (tPeriod + 1) + "reversingF1");
                    double[] _f2 = paramsFile.getArrayCell(id[
                                                              i], "T" + (tPeriod + 1) + "reversingF2");


                    if (_f1 != null && showF1) {

                        if(norm >  0.0 && normalize) {
                            MathUtil.divide(_f1,  norm);
                        }
                        p.addData(new ScatterPlot(freqs, _f1, new double[_f1.length]),
                                ReversingGratings.style1Thick[tPeriod]);
                    } 

                    if (_f2 != null && showF2) {

                        if(norm >  0.0 && normalize) {
                            MathUtil.divide(_f2,  norm);
                        }
                        p.addData(new ScatterPlot(freqs, _f2, new double[_f2.length]),
                                ReversingGratings.style2Thick[tPeriod]);
                    } 
                }




            }



//			p.setRange(0, 3, 0, 150);
//			p.setLabels(xLabel, yLabel);
//			p.setAxesType(AxisType.LOG10, AxisType.LOG10);
//			p.autoscale();

            p.loadStyle("f1f2", config);
            p.autoscale();
            return p;
        }


        public static DOG1DFourierFunction fit(
                double[] freqs, double[] f, double[] fErr,
                PlotPanel p1, String name, HashMap params) {

            DOG1DFourierFunction ff = null;
            for (double s = 0.01; s < 1; s += 0.005) {
                DOG1DFourierFunction fun;
                double max = MathUtil.max(f);
//				double s = 0.100;
//				if (name.equals("F1")) {
//				fun = new DOG1DFourierFunction(max, s, 0);
//				fun.setParameterState(2, false);
//				} else {
                fun = new DOG1DFourierFunction(max, -max / 5, s, s * 10, f[0]);
//				}

                try {
                    Fitter.fit1D(fun, freqs, f, fErr, freqs.length);
                    if (ff != null) {
                        if (fun.getChiSquared() <
                                ff.getChiSquared()
                        /*&&  fun.getS2() > 1.5 * fun.getS1()*/) {
                            ff = fun;
                        }
                    } else {
                        ff = fun;
                    }
                } catch (Exception ex) {
                }
            }
            /*
                    // get the frequency of maximum
                    Num wMax;
                    if (ff.getA2() < 0) {
//            double v = -a2 * s2 * s2 / (a1 * s1 * s1);
//            w = sqrt(2 * log(v) / (s2 * s2 - s1 * s1));
                        wMax = ff.a2().mul(ff.s2().sqr()).inv();
                        wMax = wMax.div(ff.a1()).div(ff.s1().sqr());

                        Num h = Num.zero.add(ff.s2().sqr()).sub(ff.s1().sqr());
                        wMax = wMax.ln().mul(Num.two).div(h).sqrt();
                    } else {
                        wMax = Num.zero;
                    }

                    // get the amaplitude at max
                    Num fMax = ff.getValueAt(wMax);

                    // get the 3db cutoff frequency and error
             DOG1DFourierFunction f1 = new DOG1DFourierFunction(0, 0, 0, 0, 0);
             DOG1DFourierFunction f2 = new DOG1DFourierFunction(0, 0, 0, 0, 0);
                    double error = Double.NEGATIVE_INFINITY;
                    int param = -1;
                    for (int i = 0; i < ff.getParametersCount(); i++) {
                        f1.setParameters(ff.getParameters());
                        f2.setParameters(ff.getParameters());
             f1.setParameter(i, ff.getParameter(i) + ff.getParameterError(i));
             f2.setParameter(i, ff.getParameter(i) - ff.getParameterError(i));

                        Num cutoff1 = get3dbCutoff(f1, wMax.x);
                        Num cutoff2 = get3dbCutoff(f2, wMax.x);
                        if (Math.abs(cutoff1.x - cutoff2.x) > error) {
                            error = Math.abs(cutoff1.x - cutoff2.x);
                            param = i;
                        }
                    }

                    Num wCutoff = get3dbCutoff(ff, wMax.x);
                    wCutoff = new Num(wCutoff.x, error);
             */
            p1.addToLegend(
                    "a1=" + ff.a1().toString(1) + ", " +
                    "a2=" + ff.a2().toString(1) + ", " +
                    "s1=" + ff.s1().toString(3) + ", " +
                    "s2=" + ff.s2().toString(3) + ", " +
                    "b=" + ff.b().toString(3) + ", " +
                    "\u03C7\u00B2=" + StringUtil.format(ff.getChiSquared(), 1) + ", "
            );
            /*
                    p1.addToLegend("wc = " + wCutoff.toString(1) + " @ " +
                                   ff.getParameterNames()[param]);

                    params.put(name + "wMax", wMax);
                    params.put(name + "wCutoff", wCutoff);
                    params.put(name + "fMax", fMax);
             */
            return ff;
        }


        public static PlotPanel getAverageHarmonicsPanel(
                int[] id, boolean normalize, ParametersFile pFile) throws
                IOException {

            // find the number of temporal frequencies
            int nTemporalPeriods = -1;
            for (int tPeriod = 0; tPeriod < 3; tPeriod++) {
                String f1Name = "T" + (tPeriod + 1) + "reversingF1";
                String f2Name = "T" + (tPeriod + 1) + "reversingF2";
                if (pFile.hasParameter(f1Name) && pFile.hasParameter(f2Name)) {
                    nTemporalPeriods = tPeriod + 1;
                }
            }


//			calculate overall scaling factor
            double overallsum = 0.0;
            int nn = 0;
            if(normalize) {
                for (int i = 0; i < id.length; i++) {
                    for (int tPeriod = 0; tPeriod < nTemporalPeriods; tPeriod++) {
                        double[] _f1 = pFile.getArrayCell(id[
                                                             i], "T" + (tPeriod + 1) + "reversingF1");
                        double[] _f2 = pFile.getArrayCell(id[
                                                             i], "T" + (tPeriod + 1) + "reversingF2");
                        if (_f1 != null && _f2 != null) {
                            overallsum += MathUtil.sum(_f1) + MathUtil.sum(_f2);
                            nn++;
                        }
                    }
                }
                overallsum /= nn;
            }

            PlotPanel p = new PlotPanel("F1F2");
            double[] _freqs = pFile.getArrayCell(id[0], "reversingFrequencies");
            final int nFreq = _freqs.length;
            for (int tPeriod = 0; tPeriod < 3; tPeriod++) {

                MeanVarianceCalculator[] mvcF1 = new MeanVarianceCalculator[nFreq];
                MeanVarianceCalculator[] mvcF2 = new MeanVarianceCalculator[nFreq];
                for (int periodIndex = 0; periodIndex < nFreq; periodIndex++) {
                    mvcF1[periodIndex] = new MeanVarianceCalculator();
                    mvcF2[periodIndex] = new MeanVarianceCalculator();
                }

                for (int i = 0; i < id.length; i++) {

                    double[] _f1 = pFile.getArrayCell(id[
                                                         i], "T" + (tPeriod + 1) + "reversingF1");
                    double[] _f2 = pFile.getArrayCell(id[
                                                         i], "T" + (tPeriod + 1) + "reversingF2");

                    if (_f1 != null && _f2 != null) {
                        double sum = MathUtil.sum(_f1) + MathUtil.sum(_f2);
                        if (sum > 0.0 && normalize) {
                            MathUtil.multiply(_f1, overallsum / sum);
                            MathUtil.multiply(_f2, overallsum / sum);
                        }

                        for (int j = 0; j < _f1.length; j++) {
                            mvcF1[j].add(_f1[j]);
                            mvcF2[j].add(_f2[j]);
                        }
                    }
                }

                double[] f1 = new double[nFreq];
                double[] f2 = new double[nFreq];
                double[] f1Err = new double[nFreq];
                double[] f2Err = new double[nFreq];
                for (int sPeriod = 0; sPeriod < nFreq; sPeriod++) {
                    f1[sPeriod] = mvcF1[sPeriod].getMean();
                    f1Err[sPeriod] = mvcF1[sPeriod].getMeanVariance();
                    f2[sPeriod] = mvcF2[sPeriod].getMean();
                    f2Err[sPeriod] = mvcF1[sPeriod].getMeanVariance();
                }

                if (!Double.isNaN(mvcF1[0].getMean())) {

                    p.addData(new ScatterPlot(_freqs, f1, f1Err),
                            ReversingGratings.style1Thick[tPeriod]);
                    p.addData(new ScatterPlot(_freqs, f2, f2Err),
                            ReversingGratings.style2Thick[tPeriod]);
                }

            }
            p.setLabels("Spatial Frequency (cyc/mm)", "Spike Rate (spikes/sec)");
            p.autoscale();

            return p;

        }


        public Component[] makePlots(IntegerList list, int plotType,
                TreePath classPath) {
            return null;
        }

}
