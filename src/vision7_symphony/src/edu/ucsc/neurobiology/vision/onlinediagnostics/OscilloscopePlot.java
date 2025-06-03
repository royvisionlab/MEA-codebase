package edu.ucsc.neurobiology.vision.onlinediagnostics;

import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.plot.*;


/**
 * Used by the OscilloscopeWindow class to do the real plot drawing.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class OscilloscopePlot
    implements ScatterPlotData, ChangeableData, SampleListener {

    private short[] xList, yList, yListNew;
    private String description;
    private DrawingControl drawingControl;
    private int electrode = 1;
    private int nSamples;
    private int dataBlockIndex = 0;
    private int currentSample = 0;


    public OscilloscopePlot(String description, int nSamples, int nPoints) {
        this.description = description;
        this.xList = new short[nPoints];
        this.yList = new short[nPoints];
        this.yListNew = new short[nPoints];

        for (int i = 0; i < nSamples; i++) {
            xList[i] = (short) i;
        }

        this.nSamples = nSamples;
    }


    public void setElectrode(int electrode) {
        this.electrode = electrode;
    }


    public boolean hasXErrors() {
        return false;
    }


    public boolean hasYErrors() {
        return false;
    }


    public int getElectrode() {
        return electrode;
    }


    public String getDescription() {
        return description;
    }


    public int getPointCount() {
        return xList.length;
    }


    public void getDataPoint(int i, double[] point) {
        point[0] = xList[i];
        point[1] = yList[i];
        point[2] = 0;
    }


    public final void processSample(short[] sample) {
        yListNew[currentSample] = sample[electrode];
        currentSample++;

        if (currentSample >= nSamples) {
            short[] temp = yList;
            yList = yListNew;
            yListNew = temp;

            if (drawingControl != null) {
                drawingControl.updateNeeded(this);
            }
            dataBlockIndex++;
            currentSample = 0;
        }
    }


    public void finishSampleProcessing() {
    }


    public void setDrawingControl(DrawingControl drawingControl) {
        this.drawingControl = drawingControl;
    }

}
