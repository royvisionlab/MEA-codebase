package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;


/**
 * This interface must be implemented by a class that wants to compute parameters to be
 * stored in the ParametersFile. The class has to also be listed in the
 * "Make Parameters File" group in config.xml.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface ParametersCalculator {

    /**
     * Returns the name of the calculation.
     * @return String
     */
    String getName();


    /**
     * Called by MakeParametersFile to initialize the calculator.
     *
     * @param rootPath String  the path to the folder in which the main dataset folder is
     * @param mainFilePath String the path to the main dataset
     * @throws Exception
     */
    void init(String rootPath, String mainFilePath) throws Exception;


    /**
     * Must return the param names and types, like this.
     * <pre>
        return new String[][] {
              {"RedTimeCourse", "DoubleArray"}
            , {"t1",            "Double"}
            , {"name",          "String"}
        };
     * </pre>
     *
     * @return String[][]
     */
    String[][] getParameterTypes();


    /**
     * Must return the calculated params in the same order as returned by getParameterTypes().
     * The object types must be:
     * <pre>
     * DoubleArray - double[]
     * Double - Double
     * String - String
     * </pre>
     *
     * @param c ParameterCalculatorContext the context of the calculation
     * @return Object[]
     * @throws IOException
     */
    Object[] getParameterValues(ParameterCalculatorContext c) throws IOException;
}
