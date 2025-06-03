package edu.ucsc.neurobiology.vision.plot;

import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class GrowableHistogram
    implements HistogramData {
    private String description;
    private double min;
    private double binInterval;
    private DoubleList towers;


    public GrowableHistogram(String description, double min, double binInterval) {
        this.description = description;
        this.min = min;
        this.binInterval = binInterval;

        towers = new DoubleList();
    }


    public double getBinInterval() {
        return binInterval;
    }


    public void addBin(double value) {
        towers.add(value);
    }


    public String getDescription() {
        return description;
    }


    public int getBinCount() {
        return towers.size();
    }


    public double getMin() {
        return min;
    }


    public double getMax() {
        return min + towers.size() * binInterval;
    }


    public void setDescription(String description) {
        this.description = description;
    }


    public double getBin(int i) {
        return towers.get(i);
    }
}
