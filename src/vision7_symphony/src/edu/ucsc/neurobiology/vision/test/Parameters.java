package edu.ucsc.neurobiology.vision.test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import edu.ucsc.neurobiology.vision.io.ParametersFile;
import edu.ucsc.neurobiology.vision.math.MeanVarianceCalculator;
import edu.ucsc.neurobiology.vision.util.StringUtil;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class Parameters {
    static String[][] columnNames;
    static ArrayList<Expr> expressions;
    static String normalizingClass = "All/OFF/Parasol";
    static {
        expressions = new ArrayList<Expr>();
        expressions.add(new Expr("RedTimeCourse", "DoubleArray"));
        expressions.add(new Expr("GreenTimeCourse", "DoubleArray"));
        expressions.add(new Expr("BlueTimeCourse", "DoubleArray"));
//        expressions.add(new Expr("t1", "Double"));
//        expressions.add(new Expr("t2", "Double"));
//        expressions.add(new Expr("t3", "Double"));
//        expressions.add(new Expr("a1", "Double"));
//        expressions.add(new Expr("a2", "Double"));
//        expressions.add(new Expr("a3", "Double"));
//        expressions.add(new Expr("n1", "Double"));
//        expressions.add(new Expr("n2", "Double"));
//        expressions.add(new Expr("n3", "Double"));
//        expressions.add(new Expr("dot", "Double", false));
//        expressions.add(new Expr("rl", "Double", false));

        expressions.add(new Expr("size", "2*116*((SigmaX*SigmaY)^0.5)", "Double", false));
        expressions.add(new Expr("Auto", "DoubleArray"));
        expressions.add(new Expr("nSpikes", "Double"));
        expressions.add(new Expr("reversingFrequencies", "DoubleArray"));
        expressions.add(new Expr("T1reversingF1", "DoubleArray"));
        expressions.add(new Expr("T1reversingF2", "DoubleArray"));
        expressions.add(new Expr("n", "max(T1reversingF2/T1reversingF1)", "Double"));

//        expressions.add(new Expr("flashResponse", "DoubleArray"));
//        expressions.add(new Expr("flashBinSize", "Double"));


        columnNames = new String[3 + expressions.size()][2];
        columnNames[0][0] = "ID";
        columnNames[0][1] = "Double";
        columnNames[1][0] = "classID";
        columnNames[1][1] = "String";
        columnNames[2][0] = "comment";
        columnNames[2][1] = "String";
        for (int i = 0; i < expressions.size(); i++) {
            columnNames[3 + i][0] = expressions.get(i).name;
            columnNames[3 + i][1] = expressions.get(i).type;
        }
    }


    static class Expr {
        String name, expression, type;
        boolean normalize;

        public Expr(String expression, String type) {
            this.name = expression;
            this.expression = expression;
            this.type = type;
            this.normalize = false;
        }


        public Expr(String name, String expression, String type) {
            this.name = name;
            this.expression = expression;
            this.type = type;
            this.normalize = false;
        }


        public Expr(String expression, String type, boolean normalize) {
            this(expression, type);
            this.normalize = normalize;
        }


        public Expr(String name, String expression, String type, boolean normalize) {
            this(name, expression, type);
            this.normalize = normalize;
        }
    }


    public static String toString(MeanVarianceCalculator mvc, int fDigits) {
        return StringUtil.format(mvc.getMean(), fDigits) + "\u00B1" +
            StringUtil.format(mvc.getMeanVariance(), fDigits);
    }


    static public void makeAverageParamsFile(String[] pName, String paramsFileName,
                                             String ...className) throws Exception {
        new File(paramsFileName).delete();
        ParametersFile newParamsFile = new ParametersFile(paramsFileName,
            columnNames, 1000);

        int id = 0;
        for (int prepIndex = 0; prepIndex < pName.length; prepIndex++) {
            ParametersFile pf = null;
            try {
                pf = new ParametersFile(pName[prepIndex]);
            } catch (IOException ex1) {
                System.err.println("could not open " + pName[prepIndex]);
                continue;
            }

            for (int cIndex = 0; cIndex < className.length; cIndex++) {
                int n = pf.getNeuronsInClass(className[cIndex]).length;
                if (n != 0) {
                    System.err.println(
                        n + " " + className[cIndex] + " in " + pName[prepIndex]);

                    Object[] row = new Object[columnNames.length];
                    row[0] = new Double(id);
                    row[1] = className[cIndex];
                    row[2] = pName[prepIndex];

                    for (int i = 0; i < expressions.size(); i++) {
                        Expr expr = expressions.get(i);
                        if (expr.type.equals("Double")) {
                            row[3 + i] = pf.getAverage(expr.expression, className[cIndex]);
                        } else if (expr.type.equals("DoubleArray")) {
                            row[3 + i] = pf.getArrayAverage(expr.expression,
                                className[cIndex]);
                        } else {
                            throw new Error("Unknown type");
                        }
                    }

                    newParamsFile.addRow(row);
                }
                id++;
            }
        }

        newParamsFile.close(true);
    }


    static public void makeCompoundParamsFile(String[] pName, String paramsFileName,
                                              String ...className) throws Exception {
        new File(paramsFileName).delete();
        ParametersFile newParamsFile = new ParametersFile(paramsFileName,
            columnNames, 10000);

        int newID = 0;
        all:for (int prepIndex = 0; prepIndex < pName.length; prepIndex++) {
            ParametersFile pf = new ParametersFile(pName[prepIndex]);

            double[] norm = new double[columnNames.length];
            for (int i = 0; i < expressions.size(); i++) {
                Expr expr = expressions.get(i);
                if (expr.normalize && expr.type.equals("Double")) {
                    Object v = pf.getAverage(expr.expression, normalizingClass);
                    if ( ( (Double) v).isNaN()) {
                        System.err.println("Skipping over: " + pName[prepIndex]);
                        continue all;
                    }
                    norm[3 + i] = (Double) v;
                }
            }

            for (int cIndex = 0; cIndex < className.length; cIndex++) {
                // save the values
                int[] idList = pf.getNeuronsInClass(className[cIndex]);

                for (int k = 0; k < idList.length; k++) {
                    Object[] row = new Object[columnNames.length];
                    row[0] = new Double(newID);
                    row[1] = className[cIndex];
                    row[2] = pName[prepIndex];

                    for (int i = 0; i < expressions.size(); i++) {
                        Expr expr = expressions.get(i);
                        Object v = pf.evaluate(expr.expression, idList[k]);

                        if (v instanceof Double && expr.normalize) {
                            row[3 + i] = (Double) v / norm[3 + i];
                        } else {
                            row[3 + i] = v;
                        }
                    }

                    newParamsFile.addRow(row);
                    newID++;
                }
            }
        }

        newParamsFile.close(true);
    }


    public static void main(String[] args) throws Exception {
        // OFF Small Amacrine preps
        String[] amacrine = new String[] {
                            "f:\\Good Data\\2005-04-26-0\\data009\\data009.params",
//                          "f:\\Good Data\\2005-09-09-2\\data002\\data002.params",
//                          "f:\\Good Data\\2005-05-26-2\\data003\\data003.params",
                            "f:\\Good Data\\2005-04-06-0\\data002\\data002.params",
                            "f:\\Good Data\\2005-09-27-5\\data004\\data004.params",
        };

        String[] files1 = new String[] {
                          // OFF BTY preps
                          "f:\\Good Data\\2005-04-26-0\\data009\\data009.params",
                          "f:\\Good Data\\2005-09-09-1\\data002\\data002.params",
                          "f:\\Good Data\\2005-04-14-0\\data002\\data002.params",

                          "f:\\Good Data\\2005-07-07-0\\data002\\data002.params",
                          "f:\\Good Data\\2005-07-07-6\\data000\\data000.params",
                          "f:\\Good Data\\2005-06-07-4\\data000\\data000.params",
                          "f:\\Good Data\\2005-06-07-5\\data000\\data000.params",
//                          "f:\\Good Data\\2005-09-27-5\\data004\\data004.params",
                          "f:\\Good Data\\2005-04-19-3\\data002\\data002.params",
                          "f:\\Good Data\\2005-04-14-3\\data003\\data003.params",
                          "f:\\Good Data\\2005-04-21-4\\data002\\data002.params",
                          "f:\\Good Data\\2005-04-06-0\\data002\\data002.params",
                          "f:\\Good Data\\2005-04-06-4\\data002\\data002.params",
                          "f:\\Good Data\\2005-04-26-2\\data002\\data002.params",
                          "f:\\Good Data\\2005-05-26-2\\data003\\data003.params",

//                          // Jeff's prep
//                          "F:\\Good Data\\2005-08-08-0\\data001\\data001.params"

                          // bad preps
//                          "f:\\Good Data\\2005-04-19-2\\data002\\data002.params",
        };

        String[] classes = {
//                           "All/ON/Midget",
//                           "All/OFF/Midget",
//                           "All/ON/Parasol",
//                           "All/OFF/Parasol",
//                           "All/ON/BTY",
//                           "All/OFF/BTY",
                           "All/OFF/Amacrine"
        };

        makeCompoundParamsFile(files1, "f:\\Good Data\\all.params", classes);
//        makeAverageParamsFile(files1, "f:\\Good Data\\average.params", classes);
    }
}
