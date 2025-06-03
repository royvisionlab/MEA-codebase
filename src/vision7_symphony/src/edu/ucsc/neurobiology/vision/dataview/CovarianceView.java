package edu.ucsc.neurobiology.vision.dataview;

import java.beans.*;
import java.io.*;

import javax.swing.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.plot.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class CovarianceView {
    private JComponent controller;


    private CovarianceFile covFile;



    // the PlotPanel instance on which the spikes will be plotted
    private PlotPanel covView;

    private int electrode = 33;
    private int nElectrodes;


    /**
     *
     */
    public CovarianceView(CovarianceFile covFile) {
        this.covFile = covFile;

        nElectrodes = covFile.getMaxElectrode();


        Vision app = Vision.getInstance();

        covView = new PlotPanel();
        try {
        updateCovPanel();
        } catch(IOException ex) {
            ex.printStackTrace();
        }
        app.createFrame(covView, getController(), null, covFile.getFileName());


        /*	PlotUtil.showData("Spatial Power Spectrum",
                new DoubleHistogram2D("", 0, size, 0, size, expandSymmetricMatrix(covMatrix, size)),
                new HistogramStyle());
         */}

    double[][] expandSymmetricMatrix(float[] matrix, int size) {

        double[][] expanded = new double[size][size];

        int index = 0;
        for (int i = 0; i < size; i++) {
            for (int j = i; j < size; j++) {
                expanded[i][j] = matrix[index];
                expanded[j][i] = matrix[index];
                index++;
            }
        }
        return expanded;
    }

    private void updateCovPanel() throws IOException {
        covView.removeAllData();
        float[] covMatrix = null;
        double[][] expandedCovMatrix = null;
        int size;
        covMatrix = covFile.getCovarianceMatrix(electrode);
        if(covMatrix != null) {
            size = (int) Math.round( (Math.sqrt(1 + 8 * covMatrix.length) - 1) / 2); 
            expandedCovMatrix = expandSymmetricMatrix(covMatrix, size);
        } else {
            size = 1;
            System.out.println("No covariance matrix for electrode " + electrode + ".");
            expandedCovMatrix = new double[1][1];
        }
        
        covView.addData(new DoubleHistogram2D("", 0, size, 0, size, expandedCovMatrix) , new HistogramStyle());
        covView.autoscale();
    }

    private JComponent getController() {
        ParametersTable t = new ParametersTable();

        t.addParameter(
                new IntegerParameter("Electrode", null, null, electrode, 0, nElectrodes),
                new PropertyChangeListener() {
                    public void propertyChange(PropertyChangeEvent event) {
                        electrode = ( (IntegerParameter) event.getNewValue()).getValue();
                        try {
                            updateCovPanel();
                        } catch (IOException e) {
                            Vision.reportException(e);
                        }
                    }
                }
        );

        return new JScrollPane(t);
    }
}

