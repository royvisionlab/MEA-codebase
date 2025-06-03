package edu.ucsc.neurobiology.vision.neuronviewer;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import javax.swing.JOptionPane;
import javax.swing.tree.TreePath;

import edu.ucsc.neurobiology.vision.Vision;
import edu.ucsc.neurobiology.vision.analysis.NDScatterClassifier;
import edu.ucsc.neurobiology.vision.analysis.PCACalculator;
import edu.ucsc.neurobiology.vision.io.ParametersFile;
import edu.ucsc.neurobiology.vision.math.CannotEvaluateException;
import edu.ucsc.neurobiology.vision.math.MathUtil;
import edu.ucsc.neurobiology.vision.math.TooManyIterationsException;
import edu.ucsc.neurobiology.vision.parameters.IntegerParameter;
import edu.ucsc.neurobiology.vision.parameters.ParametersTable;
import edu.ucsc.neurobiology.vision.parameters.StringParameter;
import edu.ucsc.neurobiology.vision.util.IntegerList;


/**
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class NDClassifyAction
extends CalculationAction {

    public NDClassifyAction() {
        super("ND Classify", CalculationAction.CLASS_ACTION);
    }


    public void doAction(final IntegerList list, final TreePath classTreePath) {


        ParametersTable table = configuration.showDialog(
                "NDClassification", "ND Classification", viewer.mainFrame);
        if (table == null) {
            return;
        }

        try {
            String classPath = InteractiveTree.pathToString(classTreePath);
            String cutString = makeCutString(classPath); 
            if(cutString.matches("")) return;

            HashMap<Integer, String> neuronsMap = paramsFile.evaluate("classID", cutString);


            int[] neurons = new int[neuronsMap.size()];
            String[] classNames = new String[neuronsMap.size()];

            int i=0;
            for(Integer key: neuronsMap.keySet()) {
                neurons[i] = key.intValue();
                classNames[i] = (String) neuronsMap.get(key);
                i+=1;
            }

            Object[] paramsInfo = getNDParams(cutString,
                    neurons, table, paramsFile);

            NDScatterClassifier nDScatterClassifier =
                new NDScatterClassifier(viewer, classTreePath, neurons, classNames, 
                        (double[][]) paramsInfo[0], (String[]) paramsInfo[1]);

        } catch (CannotEvaluateException ex) {
            Vision.reportException(ex);
        }

    }

    //Get the cutstring.  Includes all subclasses.
    private String makeCutString(String classPath) {
        String cutString = "";
        try {
            boolean found = false;
            HashMap<Integer, String> classLabels = paramsFile.evaluate("classID", "");
            ArrayList<String> classesList = new ArrayList<String>();
            for (Integer key : classLabels.keySet()) {
                String currentClass = classLabels.get(key);

                for (int i = 0; i < classesList.size(); i++) {
                    if (currentClass.equals(classesList.get(i))) {
                        found = true;
                    }
                }
                if (!found && currentClass.startsWith(classPath)) {
                    classesList.add(currentClass);
                }
                found = false;

            }

            if (classesList.size() == 0) {
                JOptionPane.showMessageDialog(null, "No neurons in class.",
                        "Error", JOptionPane.ERROR_MESSAGE);
                return "";
            }

            cutString = "classID==\"" + (String) classesList.get(0) + "\"";
            for (int i = 1; i < classesList.size(); i++) {
                cutString =  "classID==\"" + (String) classesList.get(i) + "\"" + "|" + cutString;
            }
        } catch (CannotEvaluateException ex) {
            Vision.reportException(ex);
        }
        return cutString;
    }


    public static Object[] getNDParams(String cutString, int[] neuronsInClass,
            ParametersTable table, ParametersFile paramsFile) throws
            CannotEvaluateException {

        //Parameters to do PCA calculation on.
        ArrayList<HashMap> params = new ArrayList<HashMap>();
        //Names of parameters
        ArrayList<String> dimensionNamesList = new ArrayList<String>();

        //Classes of each neuron
        //ArrayList classNames;
        //className.add(paramsFile.evalute("CLASSNAME", paramsFile.getIDList()

        for (int j = 0; j < 4; j++) {
            String expression = null;
            try {
                expression = ( (StringParameter) table.getParameter(
                        "Value Expression " + j)).getValue();
            } catch (NullPointerException e) {}
            if (!expression.equals("")) {
                params.add(paramsFile.evaluate(expression, paramsFile.getIDList()));
                dimensionNamesList.add(new String(expression));
            }
        }

        for (int j = 0; j < 6; j++) {
            String expression = null;
            int nComponents = 0;
            try {
                expression = ( (StringParameter) table.getParameter(
                        "PCA Expression " + j)).getValue();
                nComponents = ( (IntegerParameter) table.getParameter(
                        "N Components " + j)).getValue();

            } catch (NullPointerException e) {
                if (expression == null) {
                    break;
                }
            }

            if (!expression.equals("")) {

                HashMap<Integer, double[]> vMap = paramsFile.evaluate(expression, cutString);

                //Remove mean
                for (Iterator<Integer> iter = vMap.keySet().iterator(); iter.hasNext(); ) {
                    Object key = iter.next();
                    double[] x = (double[]) vMap.get(key);
                    if (x != null) {
                        double mean = MathUtil.mean(x);
                        for (int i = 0; i < x.length; i++) {
                            x[i] -= mean;
                        }
                    }
                }

                try {
                    PCACalculator pca = new PCACalculator(vMap);
                    for (int i = 0; i < nComponents; i++) {
                        params.add(pca.getPCAComponent(i));
                    }
                } catch (TooManyIterationsException e) {}

                for (int k = 0; k < nComponents; k++) {
                    dimensionNamesList.add(new String(expression + " " + k));
                }

            }
        }

        String[] dimensionNames = new String[dimensionNamesList.size()];
        for (int i = 0; i < dimensionNamesList.size(); i++) {
            dimensionNames[i] = (String) dimensionNamesList.get(i);
        }

        // convert params list to double[parameter][neuron]
        double[][] dataIn = new double[params.size()][neuronsInClass.length];
        for (int j = 0; j < params.size(); j++) {
            HashMap currentHashMap = params.get(j);
            for (int i = 0; i < neuronsInClass.length; i++) {
                dataIn[j][i] = ( (Double) currentHashMap.get(
                        new Integer(neuronsInClass[i]))).doubleValue();
            }
        }

        Object[] paramsInfo = { (Object) dataIn, (Object) dimensionNames};
        return paramsInfo;
    }


}
