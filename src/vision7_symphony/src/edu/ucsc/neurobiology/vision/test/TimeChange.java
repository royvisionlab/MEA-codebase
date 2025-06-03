package edu.ucsc.neurobiology.vision.test;

import static java.io.File.separator;

import java.io.File;
import java.io.IOException;

import edu.ucsc.neurobiology.vision.io.ParametersFile;
import edu.ucsc.neurobiology.vision.io.STAFile;
import edu.ucsc.neurobiology.vision.plot.PlotPanel;
import edu.ucsc.neurobiology.vision.plot.PlotUtil;
import edu.ucsc.neurobiology.vision.plot.ScatterPlot;
import edu.ucsc.neurobiology.vision.util.DoubleList;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class TimeChange {

    public static void main(String[] args) {

        String[] datasetNames = new String[] {
                                "f:\\Good Data\\2005-04-14-0",
                                "f:\\Good Data\\2005-07-07-0",
                                "f:\\Good Data\\2005-07-07-6",
                                "f:\\Good Data\\2005-06-07-4",
                                "f:\\Good Data\\2005-06-07-5",
                                "f:\\Good Data\\2005-09-27-5",
//
//                                "f:\\Good Data\\2005-04-19-2", // bad data
//                                "f:\\Good Data\\2005-04-14-3", // bad data

                                "f:\\Good Data\\2005-04-19-3",
                                "f:\\Good Data\\2005-04-21-4",
                                "f:\\Good Data\\2005-04-06-0",
                                "f:\\Good Data\\2005-04-06-4",
                                "f:\\Good Data\\2005-04-26-2", // some huge cells exist
                                "f:\\Good Data\\2005-05-26-2",

                                "f:\\Good Data\\2005-04-26-0",
                                "f:\\Good Data\\2005-09-09-1",
        };

        String[] classes = {
                           "All/ON/Midget",
                           "All/OFF/Midget",
                           "All/ON/Parasol",
                           "All/OFF/Parasol",
                           "All/ON/BTY",
                           "All/OFF/BTY"
        };

        PlotPanel[] p = new PlotPanel[classes.length];
        DoubleList[] valueList = new DoubleList[classes.length];
        DoubleList[] xList = new DoubleList[classes.length];
        for (int i = 0; i < p.length; i++) {
            p[i] = new PlotPanel();
            valueList[i] = new DoubleList();
            xList[i] = new DoubleList();
        }

        for (int j = 0; j < datasetNames.length; j++) {
            System.err.println(datasetNames[j]);
            File dataset = new File(datasetNames[j]);
            File[] folders = dataset.listFiles();

            for (int k = 0; k < classes.length; k++) {
                valueList[k].clear();
                xList[k].clear();
            }

            for (int i = 0; i < folders.length; i++) {
                String name = folders[i].getName();
                if (name.startsWith("data")) {
                    File staFile = new File(
                        datasetNames[j] + separator + name + separator + name + ".sta");
                    File paramsFile = new File(
                        datasetNames[j] + separator + name + separator + name + ".params");

                    if (!staFile.exists() || !staFile.isFile() || !staFile.canRead() ||
                        !paramsFile.exists() || !paramsFile.isFile() ||
                        !paramsFile.canRead()) {
//                    System.err.println("No STA or Params in " + name);
                        continue;
                    }

                    STAFile sta = null;
                    ParametersFile params = null;

                    try {
                        sta = new STAFile(staFile);
                    } catch (IOException ex) {
//                    System.err.println("STA file does not open in " + name);
                        continue;
                    }

                    try {
                        params = new ParametersFile(paramsFile);
                    } catch (IOException ex) {
//                    System.err.println("Params file does not open in " + name);
                        continue;
                    }

                    if (sta.getWidth() != 32 || sta.getHeight() != 16) {
//                    System.err.println("STA is not 20x20 in " + name);
                        continue;
                    }

                    for (int k = 0; k < classes.length; k++) {
//                        double v = params.getAverage("-s2/s1", classes[k]);
                        double v = params.getAverage("-zero1", classes[k]);

                        if (!Double.isNaN(v)) {
                            valueList[k].add(v);
                            xList[k].add(Integer.parseInt(name.substring(4)));
//                            System.err.println("DOT in " + name + " = " +
//                                               StringUtil.format(dot, 2));
                        }
                    }
                }
            }

            for (int k = 0; k < classes.length; k++) {
                if (xList[k].size() >= 2) {
                    p[k].addData(new ScatterPlot(xList[k].toArray(),
                                                 valueList[k].toArray(), null),
                                 "SQUARE 3 black, SOLID 1 black");
                }
            }
        }

        for (int k = 0; k < classes.length; k++) {
//            p[k].setRange( -1, 20, 0, 1.5);
            p[k].setRange( -1, 20, 40, 160);
            PlotUtil.showData(classes[k], p[k]);
        }
    }
}
