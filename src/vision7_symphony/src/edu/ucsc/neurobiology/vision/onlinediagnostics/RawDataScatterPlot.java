package edu.ucsc.neurobiology.vision.onlinediagnostics;

import java.io.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * A scatter plot implementation that interfaces directly with the RawDataFile and is
 * used to display the raw data. When selecting araw data file in the Data Manager this
 * viewer is used.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class RawDataScatterPlot
    implements ScatterPlotData, ScaleChangeListener, ChangeableData {

    private String fileNameNoExtension;
    private RawDataHeader header;
    private RawDataFile rawDataFile;
    private int electrode, nElectrodes;
    private int startSample, nSamples;
    private short[] shortBuffer;
    private DrawingControl control;
    private double samplingFrequency;


    /**
     * Create a RawDataScatterPlot which allows access to the raw
     * neurobiology data contained in the given file.
     */
    public RawDataScatterPlot(RawDataFile rawDataFile, int initialElectrode) throws
        IOException {
        this.electrode = initialElectrode;
        this.rawDataFile = rawDataFile;
        this.header = rawDataFile.getHeader();
        this.nElectrodes = header.getNumberOfElectrodes();
        this.samplingFrequency = header.getSamplingFrequency();
        fileNameNoExtension = StringUtil.removeExtension(rawDataFile.getName());
    }


    public int getElectrodesCount() {
        return header.getNumberOfElectrodes();
    }


    public void setDrawingControl(DrawingControl control) {
        this.control = control;
    }


    public String getFileName() {
        return fileNameNoExtension;
    }


    public void xScaleChanged(double xMin, double xMax) {
        startSample = (int) (xMin*samplingFrequency);
        nSamples = (int) (xMax*samplingFrequency - startSample);
        shortBuffer = new short[nSamples];

        try {
            rawDataFile.getData(electrode, startSample, shortBuffer);
        } catch (IOException e) {
            Vision.reportException(e);
        }
    }


    public void yScaleChanged(double yMin, double yMax) {
    }


    public void setElectrode(int electrode) {
        if (this.electrode == electrode) {
            return;
        } else {
            this.electrode = electrode;
            try {
                rawDataFile.getData(electrode, startSample, shortBuffer);
                
               
            } catch (IOException e) {
                Vision.reportException(e);
            }
            control.updateNeeded(this);
        }
    }


    public boolean hasXErrors() {
        return false;
    }


    public boolean hasYErrors() {
        return false;
    }


    public int getPointCount() {
        return nSamples;
    }


    public void getDataPoint(int i, double[] point) {
        point[0] = (i + startSample)/samplingFrequency;
        point[1] = shortBuffer[i];
        point[2] = 0;
    }


    public String getDescription() {
        return "Raw Data";
    }

}
