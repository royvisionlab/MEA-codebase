package edu.ucsc.neurobiology.vision.neuronviewer;

import static edu.ucsc.neurobiology.vision.plot.SymbolType.NONE;
import static java.awt.Color.black;

import java.awt.Component;
import java.io.IOException;

import javax.swing.tree.TreePath;

import edu.ucsc.neurobiology.vision.analysis.AutocorrelationCalculator;
import edu.ucsc.neurobiology.vision.math.MathUtil;
import edu.ucsc.neurobiology.vision.math.MeanVarianceCalculator;
import edu.ucsc.neurobiology.vision.math.Num;
import edu.ucsc.neurobiology.vision.plot.PlotPanel;
import edu.ucsc.neurobiology.vision.plot.ScatterPlot;
import edu.ucsc.neurobiology.vision.plot.ScatterPlotStyle;
import edu.ucsc.neurobiology.vision.util.IntegerList;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ACFPlotMaker
    extends PlotMaker {

  
    public boolean normalize = true;

    private ScatterPlotStyle autocorrelationStyle = new ScatterPlotStyle(
        "Autocorrelation", NONE, 1, black, true, black, .25f);


    public ACFPlotMaker() {
        super("ACF", CLASS_AND_NEURON_PLOT);
    }


    public Component makePlot(IntegerList list, int plotType, TreePath classPath) {
        int[] neurons = list.toArray();


        PlotPanel acfsPanel = new PlotPanel("acf");
        acfsPanel.setLabels("Time difference (ms)", "# of pairs");

        
        for (int index = 0; index < neurons.length; index++) {
            double[] acf = paramsFile.getArrayCell(neurons[index], "Auto");
            double binning = Double.NaN;
            try {
                binning = paramsFile.getDoubleCell(neurons[index], "acfBinning");
            //Value used to be hard coded.  Older parameters files do not have the value.
            } catch(NullPointerException e) {
                binning = .5;
            }
            double[] xVals = new double[acf.length];
            double[] noError = new double[acf.length];
            for(int i=0; i<xVals.length; i++) {
                xVals[i] = binning*i;
            }
            if (acf != null) {
                if (normalize) {
                    MathUtil.divide(acf, Math.sqrt(MathUtil.sumSquare(acf)));
                }
                acfsPanel.addData(new ScatterPlot(xVals, acf, noError), autocorrelationStyle);
            }
       }

        if (neuronFile != null) {
            if (neurons.length == 1) {
                try {
                    acfsPanel.addToLegend(
                        "c: " + AutocorrelationCalculator.getContamination(neurons[0],
                        neuronFile));
                    acfsPanel.addToLegend("El: " + neuronFile.getElectrode(neurons[0]));
                    acfsPanel.addToLegend("N: " + neuronFile.getSpikeCount(neurons[0]));
                    acfsPanel.setLegendVisible(true);
                } catch (IOException e) {
                    e.printStackTrace();
                }
            } else {
                try {

                    MeanVarianceCalculator mvc = new MeanVarianceCalculator();
                    for (int i = 0; i < neurons.length; i++) {
                        mvc.add(neuronFile.getSpikeCount(neurons[i]));
                    }
                    Num n = new Num(mvc.getMean(), mvc.getMeanVariance());
                    acfsPanel.addToLegend("Spikes/neuron: " + n.toString(0));
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
        acfsPanel.autoscale();

        return acfsPanel;
    }

}
