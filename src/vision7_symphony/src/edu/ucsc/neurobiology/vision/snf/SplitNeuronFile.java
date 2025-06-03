package edu.ucsc.neurobiology.vision.snf;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * Calculation used to break a neuron file in pieces after few datasets have
 * been analyzed together.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class SplitNeuronFile
    extends AbstractCalculation {


    public void startCalculation() {
        Vision app = Vision.getInstance();
    }


    public void setParameters(HashMap<String,String> parameters) {
    }


    public static void splitTimes(int[] times, int[] splitT, IntegerList[] t) {
        for (int i = 0; i < t.length; i++) {
            t[i].clear();
        }

        for (int i = 0; i < times.length; i++) {
            for (int j = splitT.length - 1; j >= 0; j--) {
                if (times[i] >= splitT[j]) {
                    t[j].add(times[i] - splitT[j]);
                    break;
                }
            }
        }

//        for (int i = 0; i < times.length; i++) {
//            for (int j = splitT.length - 1; j >= 0; j--) {
//                if (times[i] >= splitT[j]) {
//                    t[j].add(times[i]);
//                    break;
//                }
//            }
//        }
    }


    private static void splitFolder(String masterFolder) throws IOException {
        String experimentFolder = new File(masterFolder).getParentFile().getAbsolutePath();
        String masterNeuronFile =
            masterFolder + File.separator + new File(masterFolder).getName() + ".neurons";
        String floderStructureFile =
            masterFolder + File.separator + "folder-structure.txt";

//        System.out.println("experimentFolder " + experimentFolder);
//        System.out.println("masterNeuronFile " + masterNeuronFile);

        // red the folder structure
        LineNumberReader r = new LineNumberReader(new FileReader(floderStructureFile));
        ArrayList<String> names = new ArrayList<String>();
        IntegerList nSamples = new IntegerList();
        String line;
        while ( (line = r.readLine()) != null) {
            int i = line.indexOf(":");
            if (i == -1) {
                throw new IOException("Line " + r.getLineNumber() + " - missing :");
            }

            try {
                String name = line.substring(0, i).trim();
                int n = Integer.parseInt(line.substring(i + 1).trim());
                names.add(name);
                nSamples.add(n);
//                System.out.println(name + ": " + n);
            } catch (NumberFormatException e) {
                throw new IOException(
                    "Line " + r.getLineNumber() + " - wrong number format");
            }
        }

        int[] splitT = new int[nSamples.size()];
        IntegerList[] times = new IntegerList[nSamples.size()];
        for (int i = 0; i < nSamples.size(); i++) {
            times[i] = new IntegerList();
        }
        for (int i = 1; i < nSamples.size(); i++) {
            splitT[i] = splitT[i - 1] + nSamples.get(i - 1);
        }
        IOUtil.printArray(splitT);

        NeuronFile masterNF = new NeuronFile(masterNeuronFile);
        splitTimes(masterNF.getTTLTimes(), splitT, times);

        System.out.println("TTL 1: " + times[0].get(0) + ", " + times[0].get(1));
        System.out.println("TTL 2: " + times[1].get(0) + ", " + times[1].get(1));

        // create the output neuron files
        NeuronFile[] nf = new NeuronFile[nSamples.size()];
        for (int i = 0; i < nSamples.size(); i++) {
            // create the dataset folder
            new File(experimentFolder + File.separator + names.get(i)).mkdir();

            // create the neuron file
            String name =
                experimentFolder + File.separator + names.get(i) + File.separator +
                names.get(i) + VisionParams.NEURON_FILE_EXTENSION;
            //TOBEDONE
//            nf[i] = new NeuronFile(NeuronFile.INT_VERSION, name, 2000, nSamples.get(i),
//                                   times[i].toArray());
        }

        // split the times
        int[] id = masterNF.getIDList();
        for (int i = 0; i < id.length; i++) {
            System.out.println("Neuron " + i + "/" + id.length);
            splitTimes(masterNF.getSpikeTimes(id[i]), splitT, times);
            int electrode = masterNF.getElectrode(id[i]);
            for (int j = 0; j < nSamples.size(); j++) {
                nf[j].addNeuron(electrode, id[i], times[j].toArray(), times[j].size());
            }
        }

        // close neuron files
        masterNF.close();
        for (int i = 0; i < nSamples.size(); i++) {
            nf[i].close();
        }

    }


    public static void main(String[] args) throws IOException {
        Config config = Vision.getInstance().getConfig();
        String group = "SplitNeuronFile";
        ParametersTable t = config.getParameterGroup(group);
        String folder = null;

        if (args.length == t.getParametersCount()) {
            folder = args[0];
        } else {
            System.out.println(
                "Incorrect number of command line arguments: " + t.getParametersCount() +
                " required: \"datasetFolder\"");
            System.out.println("Chose parameters from the GUI...");
            t = config.showDialog(group, group, null);
            if (t == null) {
                System.out.println(
                    "You did not provide the required input. The program will now exis.");
                System.exit(1);
            } else {
                folder = t.getParameter("datasetFolder").valueAsString();
            }
        }

        splitFolder(folder);

        System.exit(1);
    }

}
