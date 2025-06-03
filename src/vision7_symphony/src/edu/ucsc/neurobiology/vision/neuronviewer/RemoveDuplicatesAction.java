package edu.ucsc.neurobiology.vision.neuronviewer;

import java.io.*;
import java.util.*;

import javax.swing.tree.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.stimulus.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class RemoveDuplicatesAction
extends CalculationAction {

    public RemoveDuplicatesAction() {
        super("Remove Duplicates", CalculationAction.CLASS_ACTION);
    }
    
    public void doAction(final IntegerList list, final TreePath classPath) {
        ParametersTable table = configuration.showDialog(
                "Remove Duplicates", "Remove Duplicates", viewer.mainFrame);
        if(table==null) {
            return;
        }

        double minRedChiSq = table.getDoubleParameter("minRedChiSq");
        double maxSeparation = table.getDoubleParameter("maxSeparation");
        double significance = table.getDoubleParameter("significance");
        double dataType = table.getEnumeratorParameter("dataType");
        
        int[] neurons;
        
        try {

            if(dataType == 0) {

                neurons = removeDuplicatesBySTA(list, maxSeparation,
                        minRedChiSq, significance);
            } else {

                neurons = removeDuplicatesByEI(list, maxSeparation,
                        minRedChiSq, significance);
            }
    
            if (neurons != null && neurons.length != 0) {
                viewer.classificationHappened(
                        neurons, classPath.pathByAddingChild(new
                                DefaultMutableTreeNode("Duplicates", true)));
            }
        } catch (IOException e) {
            Vision.reportException(e);
        }	

    }


    public int[] removeDuplicatesByEI(IntegerList list, double maxSeparation,
            double minRedChiSq, double significance) throws
            IOException {

        Vision.getInstance().sendMessage("Removing Duplicates...");

        ArrayList<Integer> neuronsInClass = new ArrayList<Integer>();
        IntegerList duplicateNeurons = new IntegerList();
        for (int i = 0; i < list.size(); i++) {
            neuronsInClass.add(list.get(i));
        }

        for (int i = neuronsInClass.size() - 1; i >= 0; i--) {
            Vision.getInstance().setProgress(
                    100 * (neuronsInClass.size() - i) / neuronsInClass.size());

            int currentNeuron = neuronsInClass.get(i);
            double currentSpikes = viewer.paramsFile.getDoubleCell(currentNeuron,
                    "nSpikes");
            double xCurrent = viewer.paramsFile.getDoubleCell(currentNeuron, "EIx0");
            double yCurrent = viewer.paramsFile.getDoubleCell(currentNeuron, "EIy0");
            if(Double.isNaN(xCurrent)) {
                float[][][] ei = viewer.imgFile.getImage(currentNeuron);
                int maxElectrode = viewer.imgFile.getMaxElectrode(ei);
                xCurrent = viewer.electrodeMap.getXPosition(maxElectrode);
                yCurrent = viewer.electrodeMap.getYPosition(maxElectrode);
            }
    //		STA currentSTA = viewer.staFile.getSTA(currentNeuron);

            for (int j = i - 1; j >= 0; j--) {
                int comparisonNeuron = (Integer) neuronsInClass.get(j);
                double comparisonSpikes = viewer.paramsFile.getDoubleCell(
                        comparisonNeuron,
                        "nSpikes");
                double xComparison = viewer.paramsFile.getDoubleCell(comparisonNeuron,
                "EIx0");
                double yComparison = viewer.paramsFile.getDoubleCell(comparisonNeuron,
                "EIy0");
                if(Double.isNaN(xComparison)) {
                    float[][][] ei = viewer.imgFile.getImage(comparisonNeuron);
                    int maxElectrode = viewer.imgFile.getMaxElectrode(ei);
                    xComparison = viewer.electrodeMap.getXPosition(maxElectrode);
                    yComparison = viewer.electrodeMap.getYPosition(maxElectrode);
                }
                double distanceSquared =
                    (xCurrent - xComparison) * (xCurrent - xComparison) +
                    (yCurrent - yComparison) * (yCurrent - yComparison);

                if (distanceSquared < maxSeparation * maxSeparation
                        || Double.isNaN(distanceSquared)) {
                    double redChiSq = viewer.imgFile.compareEIs(currentNeuron, comparisonNeuron, significance);

                    if (redChiSq < minRedChiSq) {
                        if (currentSpikes < comparisonSpikes) {
                            duplicateNeurons.add(neuronsInClass.get(i));
                            neuronsInClass.remove(neuronsInClass.get(i));

                            System.out.println(currentNeuron + " similar to " +
                                    comparisonNeuron + " reduced chi^2: " + redChiSq);
                        } else {
                            duplicateNeurons.add(neuronsInClass.get(j));
                            neuronsInClass.remove(neuronsInClass.get(j));
                            System.out.println(comparisonNeuron + " similar to " +
                                    currentNeuron + " reduced chi^2: " + redChiSq);
                        }
                        break;
                    }
                }
            }
        }
        System.out.println();
        Vision.getInstance().sendMessage("Done.");
        Vision.getInstance().setProgress(0);

        return duplicateNeurons.toArray();
    }

    public int[] removeDuplicatesBySTA(IntegerList list, double maxSeparation,
            double minRedChiSq, double significance) throws
            IOException {

        Vision.getInstance().sendMessage("Removing Duplicates...");

        ArrayList<Integer> neuronsInClass = new ArrayList<Integer>();
        IntegerList duplicateNeurons = new IntegerList();
        for (int i = 0; i < list.size(); i++) {
            neuronsInClass.add(list.get(i));
        }

        for (int i = neuronsInClass.size() - 1; i >= 0; i--) {
            Vision.getInstance().setProgress(
                    100 * (neuronsInClass.size() - i) / neuronsInClass.size());

            int currentNeuron = neuronsInClass.get(i);
            double currentSpikes = viewer.paramsFile.getDoubleCell(currentNeuron,
            "nSpikes");
            double xCurrent = viewer.paramsFile.getDoubleCell(currentNeuron, "x0");
            double yCurrent = viewer.paramsFile.getDoubleCell(currentNeuron, "y0");
            STA currentSTA = viewer.staCollection.getSTA(currentNeuron);

            for (int j = i - 1; j >= 0; j--) {
                int comparisonNeuron = (Integer) neuronsInClass.get(j);
                double comparisonSpikes = viewer.paramsFile.getDoubleCell(
                        comparisonNeuron,
                "nSpikes");
                double xComparison = viewer.paramsFile.getDoubleCell(comparisonNeuron,
                "x0");
                double yComparison = viewer.paramsFile.getDoubleCell(comparisonNeuron,
                "y0");
                double distanceSquared =
                    (xCurrent - xComparison) * (xCurrent - xComparison) +
                    (yCurrent - yComparison) * (yCurrent - yComparison);

                if (distanceSquared < maxSeparation * maxSeparation
                        || Double.isNaN(distanceSquared)) {

                    double redChiSq = currentSTA.compareSTAs(viewer.staCollection.getSTA(
                            comparisonNeuron), significance);

                    if (redChiSq < minRedChiSq) {
                        if (currentSpikes < comparisonSpikes) {
                            duplicateNeurons.add(neuronsInClass.get(i));
                            neuronsInClass.remove(neuronsInClass.get(i));

                            System.out.println(currentNeuron + " similar to " +
                                    comparisonNeuron + " reduced chi^2: " + redChiSq);
                        } else {
                            duplicateNeurons.add(neuronsInClass.get(j));
                            neuronsInClass.remove(neuronsInClass.get(j));
                            System.out.println(comparisonNeuron + " similar to " +
                                    currentNeuron + " reduced chi^2: " + redChiSq);
                        }
                        break;
                    }
                }
            }
        }

        System.out.println();
        Vision.getInstance().sendMessage("Done.");
        Vision.getInstance().setProgress(0);

        return duplicateNeurons.toArray();
    }



}
