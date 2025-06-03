package edu.ucsc.neurobiology.vision.neuronviewer;

import java.util.*;

import edu.ucsc.neurobiology.vision.plot.*;


/**
 * A plot maker that allows for other plots to be added on top of the standard plot.
 * The standard plot must be PlotPabel based and the implementer must add the
 * additionalPlots with additionalStyles manually to his standard plot.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public abstract class CompoundPlotMaker
    extends PlotMaker {

    protected ArrayList<PlotData> additionalPlots = new ArrayList<PlotData>();
    protected ArrayList<PlotStyle> additionalStyles = new ArrayList<PlotStyle>();


    public CompoundPlotMaker(String name, int plotType) {
        super(name, plotType);
    }


    /**
     * Used to register a plot with a style to this PlotMaker.
     *
     * @param plot PlotData
     * @param style PlotStyle
     */
    public void addAdditionalPlot(PlotData plot, PlotStyle style) {
//        System.err.println(style.getDescription() + " on " + this.name);

        for (int i = 0; i < additionalStyles.size(); i++) {
            if (additionalStyles.get(i).getDescription().equals(style.getDescription())) {
                throw new Error("You are adding a second style named: " +
                                style.getDescription());
            }
        }

        additionalPlots.add(plot);
        additionalStyles.add(style);
    }


    /**
     * Used to remove a plot with a given styleName from this PlotMaker.
     *
     * @param styleName String
     */
    public void removePlot(String styleName) {
        for (int i = 0; i < additionalStyles.size(); i++) {
            if (additionalStyles.get(i).getDescription().equals(styleName)) {
                additionalPlots.remove(i);
                additionalStyles.remove(i);
            }
        }
    }


}
