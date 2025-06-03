package edu.ucsc.neurobiology.vision.neuronviewer;

import java.io.*;

import javax.swing.tree.*;

import edu.ucsc.neurobiology.vision.stimulus.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class CalculateOverallTimeCourse
    extends CalculationAction {

    public CalculateOverallTimeCourse() {
        super("Calculate Overall Time Course", CalculationAction.NEURON_ACTION);
    }


    public void doAction(final IntegerList list, final TreePath classTreePath) {
        int id = list.get(0);

        try {
            STA sta = viewer.staCollection.getSTA(id);
            int w = sta.getWidth();
            int h = sta.getHeight();
            int staDepth = sta.size();

            double[][] timeFilter = new double[3][staDepth];
            float[] intensity = new float[3];

            for (int fIndex = 0; fIndex < staDepth; fIndex++) {
                ImageFrame f = sta.getFrame(fIndex);
                for (int x = 0; x < w; x++) {
                    for (int y = 0; y < h; y++) {
                        f.getPixel(x, y, intensity);
                        for (int cIndex = 0; cIndex < 3; cIndex++) {
                            timeFilter[cIndex][fIndex] += intensity[cIndex];
                        }
                    }
                }
            }

            viewer.paramsFile.setCell(id, "RedTimeCourse", timeFilter[0]);
            viewer.paramsFile.setCell(id, "GreenTimeCourse", timeFilter[1]);
            viewer.paramsFile.setCell(id, "BlueTimeCourse", timeFilter[2]);
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }


}
