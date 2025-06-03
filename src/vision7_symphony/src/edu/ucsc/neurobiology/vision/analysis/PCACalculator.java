package edu.ucsc.neurobiology.vision.analysis;

import java.util.*;

import edu.ucsc.neurobiology.vision.math.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class PCACalculator {
    private HashMap<Integer, double[]> dataMap;
    private PCA pca;


    public PCACalculator(HashMap<Integer, double[]> dataMap) throws
        TooManyIterationsException {

        this.dataMap = dataMap;

        boolean first = true;
        CovarianceMatrix c = null;
        for (Integer key : dataMap.keySet()) {
            double[] data = dataMap.get(key);
            if (data != null) {
                if (first) {
                    c = new CovarianceMatrix(data.length);
                    first = false;
                }

                boolean isNaN = false;
                for (int j = 0; j < data.length; j++) {
                    if (Double.isNaN(data[j])) {
                        isNaN = true;
                    }
                }

                if (isNaN) {
                    System.out.println("PCACalculator ---- NAN");
                } else {
                    c.addData(data);
                }
            }
        }

        pca = new PCA(c.getCovariance());
        pca.doPCA();
    }


    public HashMap<Integer,Double> getPCAComponent(int dimension) {
        HashMap<Integer,Double> pcaData = new HashMap<Integer,Double>();
        /*
                double min = Double.POSITIVE_INFINITY;
//        double max = Double.NEGATIVE_INFINITY;
                for (Integer key : dataMap.keySet()) {
                    double[] data = dataMap.get(key);
                    double x = (data != null) ? pca.project(data, dimension) : 0;
//            if (x > max) {
//                max = x;
//            }
                    if (x < min) {
                        min = x;
                    }
                }
         */
//        double s = -min / Math.abs(min);
//        double p = Math.round(Math.abs(Math.log10(Math.abs(min))));
//        double factor = Math.pow(10, s * p);

        for (Integer key : dataMap.keySet()) {
            double[] data = dataMap.get(key);
            if (data != null) {
                pcaData.put(key, pca.project(data, dimension) /* * factor*/);
            } else {
                pcaData.put(key, new Double(0));
            }
        }

        return pcaData;
    }

}
