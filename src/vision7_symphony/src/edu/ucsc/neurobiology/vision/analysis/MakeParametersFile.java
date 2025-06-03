package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;
import java.lang.reflect.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.neuronviewer.*;
import edu.ucsc.neurobiology.vision.util.*;
import edu.ucsc.neurobiology.vision.electrodemap.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class MakeParametersFile
    extends AbstractCalculation {

    private NeuronFile neuronFile;
    private STACollection staCollection, steFile;
    private PhysiologicalImagingFile eiFile;
    private ElectrodeMap electrodeMap;
    private ParametersFile paramsFile;
    private int[] neurons;
    private HashMap<String, String> parameters;

    public String mainFilePath;
    public String rootPath;
    public int nThreads;

    public void setParameters(HashMap<String, String> _parameters) {
        parameters = _parameters;

        File f = new File(parameters.get("MainFilePath")).getAbsoluteFile();
        mainFilePath = f.getAbsolutePath();
        rootPath = f.getParentFile().getAbsolutePath();
        nThreads = new Integer(parameters.get("nThreads"));
    }


    public void startCalculation() throws Exception {
        File mainFile = new File(mainFilePath);
        String dataSetName = mainFile.getName();
        String fullDataSetName = mainFilePath + File.separator + dataSetName;
        File file;

        this.neuronFile = new NeuronFile(fullDataSetName +
                                         VisionParams.NEURON_FILE_EXTENSION);
        this.neurons = neuronFile.getIDList();

        staCollection = STACollection.Factory.fromFilename(fullDataSetName + VisionParams.STA_FILE_EXTENSION);

        file = new File(fullDataSetName + ".stv");
        if (IOUtil.isValidFile(file)) {
            this.steFile = new STAFile(file.getAbsolutePath());
        }
        
        file = new File(fullDataSetName + ".ei");
        if (IOUtil.isValidFile(file)) {
            this.eiFile = new PhysiologicalImagingFile(file.getAbsolutePath());
            electrodeMap = ElectrodeMapFactory.getElectrodeMap(eiFile.getArrayID());
        }
        

        // Delete old params file
        File pFile = new File(fullDataSetName + VisionParams.PARAMS_FILE_EXTENSION);
        pFile.delete();
        
        ArrayList<ParametersCalculator>[] calculators = new ArrayList[nThreads];
        for (int thread = 0; thread < nThreads; thread++) {
            calculators[thread] = createCalculators();
        }
        

        // get the names and types
        ArrayList<String[]> types = new ArrayList<String[]> ();
        types.add(new String[] {"ID", "Double"});
        types.add(new String[] {"classID", "String"});
        for (ParametersCalculator c : calculators[0]) {
            String[][] t = c.getParameterTypes();
            for (int i = 0; i < t.length; i++) {
                types.add(t[i]);
            }
        }

        String[][] paramNamesAndTypes = new String[types.size()][];
        for (int i = 0; i < paramNamesAndTypes.length; i++) {
            paramNamesAndTypes[i] = types.get(i);
        }

        // create params file
        paramsFile = new ParametersFile(pFile.getAbsolutePath(), paramNamesAndTypes, 10000);

        Vision vision = Vision.getInstance();
        vision.startProgressBar();
        ComputeThread[] threads = new ComputeThread[nThreads];
        for (int thread = 0; thread < nThreads; thread++) {
            // FIXME: Use the ParallelUtil method
            threads[thread] = 
                new ComputeThread(thread*neurons.length/nThreads
                        ,(thread+1)*neurons.length/nThreads, calculators[thread]);
            threads[thread].start();
        }

        boolean alive = true;                                        
        while (alive) {
            Thread.sleep(100);
           
            int count = 0;
            alive = false;
            for (int thread = 0; thread < nThreads; thread++) {     	
                alive = alive || threads[thread].isAlive();
                count += threads[thread].nIndex - threads[thread].index1;
            }
            vision.setProgress(100 * (count) / (neurons.length));            
        }

        // print report
        long totalTime = 0L;
        for (int thread = 0; thread < nThreads; thread++) {
            totalTime += MathUtil.sum(threads[thread].time);
        }
        for (int i = 0; i < threads[0].time.length; i++) {
            long calcTypeTime = 0L;
            for (int thread = 0; thread < nThreads; thread++) {
                calcTypeTime += threads[thread].time[i];
            }

             System.err.println(
                   threads[0].calculatorList.get(i).getName() + ": " + ((100*calcTypeTime)/totalTime) + " %");
        }
//        double t1 = MathUtil.sum(thread1.time);
//        double t2 = MathUtil.sum(thread2.time);
//        for (int i = 0; i < thread1.time.length; i++) {
//            int percent = (int) (
//                100.0 * (thread1.time[i] / t1 + thread2.time[i] / t2) / 2.0);
//            System.err.println(
//                thread1.calculatorList.get(i).getName() + ": " + percent + " %");
//        }
        
        vision.endProgressBar();
        paramsFile.close(true);
        Vision.getInstance().getCalculationManager().calculationDone();
    }


    private ArrayList<ParametersCalculator> createCalculators() {
        // initialize calculators
        ArrayList<ParametersCalculator> calculatorList = new ArrayList<ParametersCalculator>();
        String baseName = "edu.ucsc.neurobiology.vision.analysis.";

        for (String calculatorName : parameters.keySet()) {
            if (calculatorName.endsWith("Calculator") && parameters.get(calculatorName).equals("true")) {
                try {
                    ParametersCalculator c =
                        (ParametersCalculator) Class.forName(baseName + calculatorName).newInstance();

                    // set the parameters
                    for (String name1 : parameters.keySet()) {
                        if (name1.startsWith(calculatorName + ".")) {
                            String fieldName = name1.substring(name1.indexOf(".") + 1);
                            String v = parameters.get(name1);

                            Field field = c.getClass().getField(fieldName);

                            if (field.getType().equals(int.class)) {
                                field.set(c, Integer.valueOf(v));
                            } else if (field.getType().equals(float.class)) {
                                field.set(c, Float.valueOf(v));
                            } else if (field.getType().equals(double.class)) {
                                field.set(c, Double.valueOf(v));
                            } else if (field.getType().equals(boolean.class)) {
                                field.set(c, Boolean.valueOf(v));
                            } else if (field.getType().equals(String.class)) {
                                field.set(c, v);
                            }
                        }
                    }

                    c.init(rootPath, mainFilePath);

                    calculatorList.add(c);
                } catch (Exception ex) {
                    System.err.println(
                        "Couldn't make " + calculatorName + ": " + ex.getMessage() +
                        ". See vision-output.txt");
                    ex.printStackTrace(Vision.getInstance().getFileStdOut());
                    continue;
                }
            }
        }

        return calculatorList;
    }


    class ComputeThread
        extends Thread {

        int index1, index2, nIndex;
        ArrayList<ParametersCalculator> calculatorList;
        long[] time;


        public ComputeThread(
            int index1, int index2, ArrayList<ParametersCalculator> calculatorList) throws
            IOException {

            super("MakeParametersFile ComputeThread");

            this.index1 = index1;
            this.index2 = index2;
            this.calculatorList = calculatorList;
            time = new long[calculatorList.size()];
        }


        public void run() {
            ParameterCalculatorContext context = new ParameterCalculatorContext(
                paramsFile, neuronFile, staCollection, steFile, eiFile, electrodeMap);

            for (nIndex = index1; nIndex < index2; nIndex++) {
                try {
                    context.update(neurons[nIndex]);
                } catch (IOException ex) {
                    ex.printStackTrace();
                }

                // save the values
                int index = 0;
                context.row[index++] = new Double(context.id);
                context.row[index++] = InteractiveTree.ROOT_NAME;

//                ArrayList<String[]> types = new ArrayList<String[]> ();
                for (int i = 0; i < calculatorList.size(); i++) {
                    Object[] v = null;
                    long t1 = System.nanoTime();
                    try {
                        v = calculatorList.get(i).getParameterValues(context);
                    } catch (IOException e) {
                        e.printStackTrace();
                        break;
                    }
                    long t2 = System.nanoTime();
                    time[i] += t2 - t1;

                    for (int j = 0; j < v.length; j++) {
                        context.row[index++] = v[j];
                    }
                }

                paramsFile.addRow(context.row);
            }
        }
    }

}
