package edu.ucsc.neurobiology.vision.neuronviewer;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import javax.swing.tree.TreePath;

import cern.colt.matrix.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import edu.ucsc.neurobiology.vision.Vision;
import edu.ucsc.neurobiology.vision.io.ParametersFile;
import edu.ucsc.neurobiology.vision.math.CannotEvaluateException;
import edu.ucsc.neurobiology.vision.math.CovarianceMatrix;
import edu.ucsc.neurobiology.vision.math.MeanVarianceCalculator;
import edu.ucsc.neurobiology.vision.math.PCA;
import edu.ucsc.neurobiology.vision.math.TooManyIterationsException;
import edu.ucsc.neurobiology.vision.parameters.ParametersTable;
import edu.ucsc.neurobiology.vision.plot.DoubleHistogram;
import edu.ucsc.neurobiology.vision.plot.PlotUtil;
import edu.ucsc.neurobiology.vision.util.IntegerList;
import edu.ucsc.neurobiology.vision.util.StringUtil;
import edu.ucsc.neurobiology.vision.util.Table;


/**
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class CalculateClassSeparationAction
    extends CalculationAction {

    public CalculateClassSeparationAction() {
        super("Calculate Class Separation", CalculationAction.CLASS_ACTION);
    }


    public void doAction(final IntegerList list, final TreePath classTreePath) {
        try {
            calculateClassSeparation(classTreePath);
        } catch (CannotEvaluateException ex) {
            Vision.reportException(ex);
        }
    }


    public double[] calculateClassSeparation(TreePath classTreePath) throws
        CannotEvaluateException {
        ArrayList<Double> separationList = new ArrayList<Double> ();
        DoubleHistogram hist = new DoubleHistogram("", 0, 100, 3);

        String classPath = InteractiveTree.pathToString(classTreePath);
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

        String[] classes = new String[classesList.size()];
        for (int i = 0; i < classesList.size(); i++) {
            classes[i] = (String) classesList.get(i);
        }

        ParametersTable table = configuration.showDialog(
            "NDClassification", "ND Classification", viewer.mainFrame);
        if (table == null) {
            return null;
        }
//          ParametersTable table = configuration.getParameterGroup("NDClassification");
        Table table2 = new Table(classes.length, classes.length + 1);
        table2.setTitle(0, "");

        //data[parameter][neuron]
        double[][] data = null;
        for (int i = 0; i < classes.length; i++) {
            //set column headings
            table2.setTitle(i + 1, classes[i]);

            //set row headings
            table2.setCell(i, 0, classes[classes.length - 1 - i]);
            //set diagonal values to zero.
            table2.setCell(classes.length - 1 - i, i + 1, 0.0);

            for (int j = i + 1; j < classes.length; j++) {
                Object[] neuronsByClassWithNames = getNeuronsInClasses(new String[] {
                    classes[i], classes[j]}, paramsFile);
                int[] neuronsInClasses = (int[]) neuronsByClassWithNames[0];
                String[] classNamesByNeuron = (String[]) neuronsByClassWithNames[1];
                String str = "classID==\"" + classes[i] + "\"" + "|" +
                             "classID==\"" + classes[j] + "\"";
                String[] parts = classes[i].split("/");
                boolean weak = false;
                for (int k = 0; k < parts.length; k++) {

                    if (parts[k].equals("out") || parts[k].equals("Axons") ||
                        parts[k].equals("Duplicates") || parts[k].equals("contaminated")
                        || parts[k].equals("mismapped") || parts[k].equals("Weak")
                        || parts[k].equals("DS") || parts[k].equals("lost")
                        || parts[k].equals("minimal classes")) {
                        weak = true;
                        break;
                    }
                }
                if (weak == true) {
                    continue;
                }

                parts = classes[j].split("/");
                weak = false;
                for (int k = 0; k < parts.length; k++) {
                    if (parts[k].equals("out") || parts[k].equals("Axons") ||
                        parts[k].equals("Duplicates") || parts[k].equals("contaminated")
                        || parts[k].equals("mismapped") || parts[k].equals("Weak")
                        || parts[k].equals("DS") || parts[k].equals("lost") ||
                        parts[k].equals("minimal classes")) {
                        weak = true;

                        break;
                    }
                }

                if (weak == true) {
                    continue;
                }

//             if(str.equals("classID==\"" + "All/Off/DS/1" + "\"" + "|" +
//                          "classID==\"" + "All/On/Large/A" + "\"")) {
//                System.out.println(str);
                Object[] dataObject = new Object[2];
                dataObject = NDClassifyAction.getNDParams(str, neuronsInClasses, table,
                    paramsFile);

                //data[param][neuron]
                data = (double[][]) dataObject[0];

                //Calculate mean and covariance for each cluster, pairwise.
                //covCalculator[class]
                //"class" = 3 is both classes

                CovarianceMatrix[] covCalculator = new CovarianceMatrix[3];
                for (int ii = 0; ii < 3; ii++) {
                    covCalculator[ii] = new CovarianceMatrix(data.length);
                }

                for (int neuron = 0; neuron < neuronsInClasses.length; neuron++) {
                    boolean badNeuron = false;
                    //if neuron is in class zero
                    if (classes[i] == classNamesByNeuron[neuron]) {
                        double[] neuronData = new double[data.length];

                        for (int param = 0; param < data.length; param++) {
                            neuronData[param] = data[param][neuron];
                            if (Double.isNaN(data[param][neuron]) ||
                                data[param][neuron] == 0.0) {
                                badNeuron = true;
                            }
                        }
                        if (!badNeuron) {
                            covCalculator[0].addData(neuronData);
                            covCalculator[2].addData(neuronData);
                        }
                        //else if neuron is in class one
                    } else if (classes[j] == classNamesByNeuron[neuron]) {
                        double[] neuronData = new double[data.length];

                        for (int param = 0; param < data.length; param++) {
                            neuronData[param] = data[param][neuron];
                            if (Double.isNaN(data[param][neuron]) ||
                                data[param][neuron] == 0.0) {
                                badNeuron = true;
                            }
                        }
                        if (!badNeuron) {
                            covCalculator[1].addData(neuronData);
                            covCalculator[2].addData(neuronData);
                        }
                        //else, bug in code.
                    } else {
                        System.out.println(classes[i]);
                        System.out.println(classes[j]);
                        System.out.println(classNamesByNeuron[neuron]);
                        throw new IllegalStateException("Error");
                    }
                }

                Algebra alg = new Algebra();
                DenseDoubleMatrix2D[] covMatrices = new DenseDoubleMatrix2D[3];
                DenseDoubleMatrix1D[] meanMatrices = new DenseDoubleMatrix1D[2];
                covMatrices[0] = new DenseDoubleMatrix2D(data.length, data.length);
                covMatrices[1] = new DenseDoubleMatrix2D(data.length, data.length);
                covMatrices[2] = new DenseDoubleMatrix2D(data.length, data.length);
                meanMatrices[0] = new DenseDoubleMatrix1D(data.length);
                meanMatrices[1] = new DenseDoubleMatrix1D(data.length);

                DenseDoubleMatrix1D direction = new DenseDoubleMatrix1D(data.length);

//                float[] covar = covCalculator[2].getCovariance();
//             int kk = 0;
//            for (int ii = 0; ii < data.length; ii++) {
//                   meanMatrices[1].setQuick(ii, means[1][ii]);

//                for (int jj = ii; jj < data.length; jj++) {
//                    covMatrices[2].setQuick(ii, jj, covar[kk]);
//                    covMatrices[2].setQuick(jj, ii, covar[kk]);
//                    kk++;
//                }
//            }

                //find whitener that diagonalizes full covariance matrix, and
                //scales eigenvalues to all be the same.
                PCA pca = new PCA(covCalculator[2].getCovariance());
                try {
                    pca.doPCA();
                } catch (TooManyIterationsException e) {
                    e.printStackTrace();
                }
                DenseDoubleMatrix2D whitener = new DenseDoubleMatrix2D(data.length,
                    data.length);
                DenseDoubleMatrix2D covEigen = new DenseDoubleMatrix2D(pca.
                    getEigenVectors());

                for (int ii = 0; ii < data.length; ii++) {
                    whitener.setQuick(ii, ii, 1 / Math.sqrt(pca.getEigenValue(ii)));
                }

                whitener = (DenseDoubleMatrix2D) alg.mult(whitener, covEigen);

                float[] covar = covCalculator[0].getCovariance();
                float[][] means = new float[2][];

                //Convert matrices to Algebra ready form
                means[0] = covCalculator[0].getMeans();
                int kk = 0;
                for (int ii = 0; ii < data.length; ii++) {
                    meanMatrices[0].setQuick(ii, means[0][ii]);
                    for (int jj = ii; jj < data.length; jj++) {
                        covMatrices[0].setQuick(ii, jj, covar[kk]);
                        covMatrices[0].setQuick(jj, ii, covar[kk]);
                        kk++;
                    }
                }

                covar = covCalculator[1].getCovariance();
                means[1] = covCalculator[1].getMeans();
                kk = 0;
                for (int ii = 0; ii < data.length; ii++) {
                    meanMatrices[1].setQuick(ii, means[1][ii]);
                    direction.setQuick(ii, means[1][ii] - means[0][ii]);
                    for (int jj = ii; jj < data.length; jj++) {
                        covMatrices[1].setQuick(ii, jj, covar[kk]);
                        covMatrices[1].setQuick(jj, ii, covar[kk]);
                        kk++;
                    }
                }

//                 covar = covCalculator[2].getCovariance();
//                kk = 0;
//               for (int ii = 0; ii < data.length; ii++) {
//                   meanMatrices[1].setQuick(ii, means[1][ii]);

//                   for (int jj = ii; jj < data.length; jj++) {
//                       covMatrices[2].setQuick(ii, jj, covar[kk]);
//                       covMatrices[2].setQuick(jj, ii, covar[kk]);
//                       kk++;
//                   }
//               }

                //normalize the direction vector
                double normalizer = Math.sqrt(alg.mult(direction, direction));
                for (int ii = 0; ii < data.length; ii++) {
                    direction.setQuick(ii, direction.getQuick(ii) / normalizer);
                }

                //convert all objects to white space
                direction = (DenseDoubleMatrix1D) alg.mult(whitener, direction);
                meanMatrices[0] = (DenseDoubleMatrix1D) alg.mult(whitener, meanMatrices[0]);
                meanMatrices[1] = (DenseDoubleMatrix1D) alg.mult(whitener, meanMatrices[1]);
                covMatrices[0] = (DenseDoubleMatrix2D) alg.mult(whitener,
                    alg.mult(covMatrices[0], alg.transpose(whitener)));
                covMatrices[1] = (DenseDoubleMatrix2D) alg.mult(whitener,
                    alg.mult(covMatrices[1], alg.transpose(whitener)));
//            covMatrices[2] = (DenseDoubleMatrix2D) alg.mult(whitener,alg.mult(covMatrices[2],
//alg.transpose(whitener)));

                double distance = alg.mult(meanMatrices[0], direction) -
                                  alg.mult(meanMatrices[1], direction);
                double scaler = Math.sqrt(alg.mult(alg.mult(covMatrices[0], direction),
                    direction)
                                          +
                                          alg.mult(alg.mult(covMatrices[1], direction),
                    direction));

//                System.out.println("Covariance 0");
//             System.out.println(covMatrices[0].getQuick(0, 0) + "  " +
//                                covMatrices[0].getQuick(1, 0) + "  " +
//                                covMatrices[0].getQuick(0, 1) + "  " +
//                                covMatrices[0].getQuick(1, 1));
//
//              System.out.println("Covariance 1");
//              System.out.println(covMatrices[1].getQuick(0, 0) + "  " +
//                                 covMatrices[1].getQuick(1, 0) + "  " +
//                                 covMatrices[1].getQuick(0, 1) + "  " +
//                                 covMatrices[1].getQuick(1, 1));
//
//              System.out.println("Covariance 2");
//             System.out.println(covMatrices[2].getQuick(0, 0) + "  " +
//                                covMatrices[2].getQuick(1, 0) + "  " +
//                                covMatrices[2].getQuick(0, 1) + "  " +
//                                covMatrices[2].getQuick(1, 1));



//                System.out.println("Mean 0");
//               System.out.println(meanMatrices[0].getQuick(0) + "  " +
//                                  meanMatrices[0].getQuick(1));
//
//                System.out.println("Mean 1");
//                System.out.println(meanMatrices[1].getQuick(0) + "  " +
//                                   meanMatrices[1].getQuick(1));
//
//                System.out.println("Direction");
//                     System.out.println(direction.getQuick(0) + "  " +
//                                        direction.getQuick(1));
//
//            System.out.println(alg.mult(meanMatrices[0], direction));
//
//            System.out.println(alg.mult(meanMatrices[1], direction));
//           System.out.println(distance);
//           System.out.println("****");
//           System.out.println(alg.mult(alg.mult(covMatrices[0], direction), direction));
//           System.out.println(alg.mult(alg.mult(covMatrices[1], direction), direction));
                double separation = Math.abs(distance / scaler);
                table2.setCell(classes.length - j - 1, i + 1,
                               StringUtil.format(separation, 2));
                hist.fill(separation, 1.0);
                separationList.add(new Double(separation));
            }
        }
        table2.draw(System.err);
        try {
            table2.draw(new PrintStream(new File(viewer.filePathRoot + ".table")));
        } catch (IOException e) {
            e.printStackTrace();
        }
        PlotUtil.showData("", hist);

        double[] separationArray = new double[separationList.size()];
        for (int i = 0; i < separationArray.length; i++) {
            separationArray[i] = ( (Double) separationList.get(i)).doubleValue();
        }
        return separationArray;
    }


    @SuppressWarnings("unused")
    private void calculateClassSeparationParamByParam() throws CannotEvaluateException {
        boolean found = false;

        HashMap<Integer, String> classLabels = paramsFile.evaluate("classID", "");

        ArrayList<String> classesList = new ArrayList<String>();
        for (Iterator<Integer> iter = classLabels.keySet().iterator(); iter.hasNext(); ) {
            Integer key = (Integer) iter.next();
            String currentClass = (String) classLabels.get(key);
            for (int i = 0; i < classesList.size(); i++) {
                if (currentClass.equals( (String) classesList.get(i))) {
                    found = true;
                }
            }
            if (!found) {
                classesList.add(currentClass);
            }
            found = false;
        }
        String[] classes = new String[classesList.size()];
        for (int i = 0; i < classesList.size(); i++) {
            classes[i] = (String) classesList.get(i);
        }

        configuration.showDialog("NDClassification", "ND Classification", viewer.mainFrame);

        //data[parameter][neuron]
        double[][] data = null;
        for (int i = 0; i < classes.length; i++) {
            for (int j = i + 1; j < classes.length; j++) {
//         data[i] = (double[][]) getNDParams(classes[i], getNeuronsInClass(classes[i], paramsFile), table)[0];

                Object[] neuronsByClassWithNames = getNeuronsInClasses(new String[] {
                    classes[i], classes[j]}
                    , paramsFile);
                int[] neuronsInClasses = (int[]) neuronsByClassWithNames[0];
                String[] classNamesByNeuron = (String[]) neuronsByClassWithNames[1];
                String str = "classID==\"" + classes[i] + "\"" + "|" +
                             "classID==\"" + classes[j] + "\"";
//             if(str.equals("classID==\"" + "All/On/TwoPeaks/Slow" + "\"" + "|" +
//                          "classID==\"" + "All/On/Sustained/A" + "\"")) {
//
                System.out.println(str);
                Object[] dataObject = new Object[2];
//                dataObject = getNDParams(str,
//                                         neuronsInClasses, table);
                data = (double[][]) dataObject[0];

                //Calculate mean and variance for each cluster, pairwise.
                //mvc[class][parameter]
                MeanVarianceCalculator[][] mvc = new MeanVarianceCalculator[2][data.
                                                 length];
                for (int ii = 0; ii < 2; ii++) {
                    for (int jj = 0; jj < data.length; jj++) {
                        mvc[ii][jj] = new MeanVarianceCalculator();
                    }
                }

                for (int neuron = 0; neuron < neuronsInClasses.length; neuron++) {
                    if (classes[i] == classNamesByNeuron[neuron]) {
                        for (int param = 0; param < mvc[0].length; param++) {
                            if (!Double.isNaN(data[param][neuron]) &&
                                data[param][neuron] != 0.0) {
                                mvc[0][param].add(data[param][neuron]);
                            }
                        }
                    } else if (classes[j] == classNamesByNeuron[neuron]) {
                        for (int param = 0; param < mvc[0].length; param++) {
                            if (!Double.isNaN(data[param][neuron]) &&
                                data[param][neuron] != 0.0) {
                                mvc[1][param].add(data[param][neuron]);
                            }
                        }
                    } else {
                        System.out.println(classes[i]);
                        System.out.println(classes[j]);
                        System.out.println(classNamesByNeuron[neuron]);
                        throw new IllegalStateException("Error");
                    }
                }

                double maxSeparation = 0.0;
                for (int param = 0; param < mvc[0].length; param++) {
                    maxSeparation = Math.max(Math.abs( (mvc[0][param].getMean() -
                        mvc[1][param].getMean()) / (
                            mvc[0][param].getStandardDeviation() +
                            mvc[1][param].getStandardDeviation())),
                                             maxSeparation);
                }
                System.out.println(maxSeparation);
//            }
            }
        }

        return;
    }


    public static Object[] getNeuronsInClasses(String[] classNames, ParametersFile params) {
        HashMap<Integer,? extends Object> classLabels = params.getClassIDs();

        IntegerList idList = new IntegerList();
        ArrayList<String> nameList = new ArrayList<String>();
        for (Iterator<Integer> iter = classLabels.keySet().iterator(); iter.hasNext(); ) {
            Integer key = iter.next();
            String classType = (String) classLabels.get(key);
            for (int i = 0; i < classNames.length; i++) {
                if (classType.equals(classNames[i])) {
                    idList.add(key.intValue());
                    nameList.add(classNames[i]);
                }
            }
        }
        String[] classNamesByNeuron = new String[nameList.size()];
        for (int i = 0; i < nameList.size(); i++) {
            classNamesByNeuron[i] = (String) nameList.get(i);
        }
        return new Object[] {idList.toArray(), classNamesByNeuron};
    }


}
