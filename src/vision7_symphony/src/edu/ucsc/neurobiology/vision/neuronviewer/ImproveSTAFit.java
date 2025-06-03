package edu.ucsc.neurobiology.vision.neuronviewer;

import java.io.*;

import javax.swing.tree.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.math.fitting.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.stimulus.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ImproveSTAFit
    extends CalculationAction {

    public ImproveSTAFit() {
        super("Improve STA Fit", CalculationAction.NEURON_ACTION);
    }


    public void doAction(final IntegerList list, final TreePath classTreePath) {
        ParametersTable table = Vision.getConfig().showDialog(
            "Improve STA Fit", "Improve STA Fit", viewer);
        if (table == null) {
            return;
        }
        int frame1 = table.getIntParameter("Frame 1");
        int frame2 = table.getIntParameter("Frame 2");
        double x0 = table.getDoubleParameter("x0");
        double y0 = table.getDoubleParameter("y0");

        try {
            int id = list.get(0);
            STA sta = viewer.staCollection.getSTA(id);
            double stixelWidth = sta.getStixelWidth();
            double stixelHeight = sta.getStixelHeight();
            int w = sta.getWidth();
            int h = sta.getHeight();

            STAFrame frame = new STAFrame(w, h, stixelWidth, stixelHeight);
            for (int k = frame1; k <= frame2; k++) {
                ImageFrame f = sta.getFrame(k);
                for (int i = 0; i < w; i++) {
                    for (int j = 0; j < h; j++) {
                        for (int c = 0; c < 3; c++) {
                            frame.addToPixel(i, j, c, f.getPixel(i, j, c));
                            frame.setPixelError(
                                i, j, c,
                                frame.getPixelError(i, j, c) + f.getPixelError(i, j, c));
                        }
                    }
                }
            }

            Gaussian2DFunction g = frame.improveFit(1, x0, y0);
            System.err.println(g);
            ParametricEllipse e = new ParametricEllipse(
                g.getX0() * stixelWidth, g.getY0() * stixelHeight, g.getSigmaX() * stixelWidth,
                g.getSigmaY() * stixelHeight, -g.getTheta());

            STA sta1 = new STA(new STAFrame[] {frame}, 1);
            PlotPanel p = STAPlotMaker.makeSTAPanel(sta1, false, 2, 3, false, true, false,
                viewer.globalsFile, e);
            PlotUtil.showData("", p);

            viewer.paramsFile.setCell(id, "x0", g.getX0());
            viewer.paramsFile.setCell(id, "y0", g.getY0());
            viewer.paramsFile.setCell(id, "SigmaX", g.getSigmaX());
            viewer.paramsFile.setCell(id, "SigmaY", g.getSigmaY());
            viewer.paramsFile.setCell(id, "Theta", g.getTheta());
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

}