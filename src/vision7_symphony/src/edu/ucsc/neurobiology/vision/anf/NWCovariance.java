package edu.ucsc.neurobiology.vision.anf;

import java.io.*;
import java.util.*;

import cern.colt.matrix.*;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.plot.DoubleHistogram2D;
import edu.ucsc.neurobiology.vision.plot.HistogramStyle;
import edu.ucsc.neurobiology.vision.plot.PlotUtil;


/**
 * Compute noise whitened covariance matrices. This code is integrated into the default
 * automated neuron finding (anf) system. The process of whitening occurs immediately after
 * standard spike finding has completed. It uses the noise matrix calculated during spike finding 
 * to whiten the standard covariance matrix. PCA is then performed on the whitened covariance matrix that is returned from
 * this class. 
 * 
 * @author tamachado@salk.edu (7/29/08)
 * @author mgrivich@salk.edu (3/29/10), moved noise covariance calculation to spike finding for a substantial performance boost.
 */
public class NWCovariance
extends AbstractCalculation  {

    private ElectrodeUsage electrodeUsage;

    private int   nElectrodes;//electrodes in this dataset
    private int   nThreads;   //threads for BLAS operations to use
    private int   nlPoints;   //spike window size to left of spike
    private int   nrPoints;   //spike window size to right of spike

    private boolean[] badElectrodes;
    private boolean[][] adjacentMap;
    private int[] neighborCount;
    private String rawPath;
    private String outputPath;
    private Vision app;
    private RawDataHeader header;
    private DoubleMatrix2D[] reducedWhiteMatrix;

    //Extension of standard covariance file
    public static final String UNWHITENED_EXT = ".cov";

    //Extension of whitened covariance file
    public static final String WHITENED_EXT = ".wcov";

    public void startCalculation() throws Exception {
        app.startProgressBar();
        //Create Noise Whitened Projections
        try {
            //Read the raw data header
            RawDataWrapper rFile = new RawDataWrapper(rawPath);
            header = rFile.getHeader();
            nElectrodes = header.getNumberOfElectrodes();

            computeCovarianceMatrices();
            computeNWCovarianceFile();
        } catch (IOException e) {
            e.printStackTrace();
            throw new RuntimeException(e);
        }
    }


    /**
     * This is computes the whitened covariance matrices from noise data. 
     * 
     * NOTE: Because of the way the spikeAligner in vision works, each spike vector
     * is essentially shortened by two. Look at the spike finding and spike aligner
     * classes in vision for more information.
     * 
     * Consequently, the size of the covariance matrix in vision is:
     * 
     * covMatrix[electrode] = ((vector length - 2) * neighbors)^2
     * 
     * 
     * @throws IOException upon failure reading the spikes files or writing covariance files
     */
    private void computeCovarianceMatrices() throws IOException {
        neighborCount = new int[nElectrodes];
        badElectrodes = new boolean[nElectrodes];       
        
        //Parse spikes covariance file
        CovarianceFile oldCovFile = new CovarianceFile(outputPath + new File(outputPath).getName() +  UNWHITENED_EXT);
        SpikeFile spikeFile  = new SpikeFile(outputPath + new File(outputPath).getName() + ".spikes");
        
        VisionHeader oldCovHeader = oldCovFile.getHeader();
        this.nlPoints = oldCovHeader.nlPoints;
        this.nrPoints = oldCovHeader.nrPoints;
        this.electrodeUsage = oldCovHeader.electrodeUsage;
        oldCovFile.close();
        
        File outputDir = new File(outputPath);
        
        //Compute the covariance matrices
        String outputName =
            outputDir.getAbsolutePath() + File.separator + outputDir.getName();
        
        String noiseCovarianceFileName = outputName + ".ncov";
        String covarianceFileName = outputName + ".cov";
        
        try {
            VisionHeader sHeader = (new CovarianceFile(covarianceFileName)).getHeader();
            ElectrodeMap map = ElectrodeMapFactory.getElectrodeMap(sHeader.arrayID);
            int[][] adjacent = new int[nElectrodes][];
            adjacentMap = new boolean[nElectrodes][];
            
            for(int electrode = 1; electrode <nElectrodes; electrode++) {
                if(spikeFile.getSpikesCount(electrode) == 0) badElectrodes[electrode] = true;
            }

            //Build matrix of spike values:
            //- We use a double matrix because we're going to call linear algebraic
            //  methods on this matrix that will need to do floating point arithmetic
            //- Matrix is implemented in row major order, so it's 
            //  faster for each row to be a spike vector
            //- There exists a noise covariance matrix for each electrode
            for (int electrode = 1; electrode < nElectrodes; electrode++) {
                
                adjacent[electrode] = map.getAdjacentsTo(electrode, electrodeUsage.ordinal());
                adjacentMap[electrode] = new boolean[adjacent[electrode].length];
                neighborCount[electrode] = 0;

                //Ignore bad/disconnected adjacent electrodes-- they might make the covariance matrix singular.
                //
                //The standard covariance matrix calculation code does not throw away these neighbors
                //consequently, we must resize our matrices to be [nPoints][nNeighbors] in the final
                //step of whitening when we combine the noise covariance matrices with the
                //standard covariance matrices. This is why we save out a list of the ignored neighbors
                //here (adjacentMap). Instead of saving the electrode number, we save a bit vector of 
                //maxNeighborLength to show which areas of the covariance matrix must be filled with zeros
                for (int neighbor = 0; neighbor < adjacent[electrode].length; neighbor++) {

                    //See if neighbor is bad. It must be bad if the first selected spike time is zero
                    //since that implies that the spike window begins before the first sample. The random
                    //spike selection code does not allow such a spike to be chosen.

                    int n = adjacent[electrode][neighbor]; 
                    if (badElectrodes[n] == true) {
                        adjacentMap[electrode][neighbor] = false;
                    } else {
                        adjacentMap[electrode][neighbor] = true;
                        neighborCount[electrode]++;
                    }
                }
            }

            //Take the covariance and inverse of each noise matrix:
            CovarianceFile covarianceFile = new CovarianceFile(noiseCovarianceFileName);

            reducedWhiteMatrix = new DoubleMatrix2D[nElectrodes];
            for (int electrode = 1; electrode < nElectrodes; electrode++) {
                if (badElectrodes[electrode]) continue;
                app.setProgress((int) (50.0 * electrode /(nElectrodes - 1)));
                
                DoubleMatrix2D covMatrix = visionCovarianceToColt(covarianceFile.getCovarianceMatrix(electrode),electrode);
                int size = covMatrix.rows();  //= covMatrix.columns();         
            
                //Remove bad channels.  If you whiten a covariance matrix
                //with zeros, you get infinities.  These will be added back in later.
                int windowLength = nlPoints + nrPoints + 1 -2;
                int reducedSize = windowLength*neighborCount[electrode];
                DenseDoubleMatrix2D reducedCovMatrix = new DenseDoubleMatrix2D(reducedSize, reducedSize);
                int skippedColumns = 0, k=0;
                int skippedRows = 0;
                for (int r = 0; r < size; r++) {
                    skippedColumns = 0;
                    k = r/windowLength;
                    if (adjacentMap[electrode][k] == true) {
                        for (int c = 0; c < size; c++) {
                            k = c/windowLength;
                            if (adjacentMap[electrode][k] == true) {
                                reducedCovMatrix.setQuick(r-skippedRows, c - skippedColumns, covMatrix.get(r, c));
                            } else {
                                skippedColumns++;
                            }
                        }
                    } else {
                        skippedRows++;
                    }
                }
            

                try {
                    reducedWhiteMatrix[electrode] = whitenMatrix(reducedCovMatrix);
                } catch(Exception e) {
                    badElectrodes[electrode] = true;
                }
            }

        } catch(IllegalArgumentException e) {
            e.printStackTrace();
            throw new IOException("Error during noise matrix calculation!");
        }
        spikeFile.close();
    }

    /**
     * This is a helper method for computeCovarianceMatrices() that actually computes
     * the term inverse( sqrtm ( cov(noise) ) ), where sqrtm is the square root of a matrix
     * computed by eigendecomposition. covMatrix is the covariance matrix of the noise.
     * 
     * 
     * @throws IllegalArgumentException if one of the linear algebraic operations fails
     */
    private DoubleMatrix2D whitenMatrix(DoubleMatrix2D covMatrix) {
        Algebra a = new Algebra();

        EigenvalueDecomposition ev = new EigenvalueDecomposition(covMatrix);

        //We need to compute the term inverse( sqrtm ( cov(noise) ) )
        //this is equivalent to inverse( V * sqrt(D) * inverse(V))
        //where V and D come from the eigendecomposition of the covariance matrix
        DoubleMatrix2D diagMatrix = ev.getD();

        //We are going to take the sqrt of the eigenvalues. Since this operation
        //does not use Colt, we need to manually round everything to the tolerance
        //used in doing the linear algebra. If we don't do this, we could end up
        //with very small negative numbers (i.e. -1e-15)! This would obviously cause problems,
        //although the problem should only appear if the matrix is already very ill-conditioned
        Double tolerance = a.property().tolerance();

        for (int i = 0; i < diagMatrix.rows(); i++) {
            diagMatrix.setQuick(i, i, Math.sqrt(Math.rint(diagMatrix.getQuick(i, i)/tolerance)*tolerance));
        }

        DoubleMatrix2D temp = a.mult(ev.getV(), diagMatrix);

        return a.inverse(a.mult(temp, a.inverse(ev.getV())));
    }

    /**
     * 
     * This is the final step in the whitening process: combining the two covariance
     * matrices into a single covariance matrix to be serialized.
     * 
     * The whitening transformation used is as follows:
     * 
     * cov(new) = inverse( sqrtm ( cov(noise) ) ) * cov(data) * inverse( sqrtm ( cov(noise) ) )
     *  
     * The first n eigenvectors of the new covariance matrix are whitened principal components.
     * We represent the data in this new space and write out a projections file to the output path. 
     * Linear algebraic operations are done using the Colt API (which is also used by Vision).
     * 
     * @throws IOException upon failure reading the covariance files or writing the prj file
     */
    private void computeNWCovarianceFile() throws IOException {

        if (reducedWhiteMatrix == null) {
            System.out.println("Warning: Whitened Noise Matrices Not In Memory, Recalculating...");
            computeCovarianceMatrices();
        }

        String unwhitenedCovPath = outputPath + new File(outputPath).getName() +  UNWHITENED_EXT;
        String whitenedCovPath   = outputPath + new File(outputPath).getName() +  WHITENED_EXT;
        String spikesPath        = outputPath + new File(outputPath).getName() + ".spikes";

        File fcovn = new File(unwhitenedCovPath);
        File fspk  = new File(spikesPath);

        if (!fspk.canRead() || !fcovn.canRead() ) {
            throw new IOException("Cannot Read Spike Covariance Matrix!");
        }

        SpikeFile spikesFile = new SpikeFile(spikesPath);

        //Open the standard spikes covariance file
        CovarianceFile covariance = new CovarianceFile(unwhitenedCovPath);

        //Create a new covariance file to store the whitened covariance matrices
        VisionHeader covHeader = covariance.getHeader();
        covHeader.covarianceType = VisionHeader.NOISE_WHITENED;
        CovarianceFile whiteCovariance = new CovarianceFile(whitenedCovPath, covHeader);

        



        try {
            for (int electrode = 1; electrode < nElectrodes; electrode++) {

                app.setProgress((int) (.5 + 50.0 * electrode/(nElectrodes - 1)));

                //Parse spike covariance matrix
                float[] covMatrix = covariance.getCovarianceMatrix(electrode);

                //Skip bad electrodes
                if (reducedWhiteMatrix[electrode] == null || covMatrix == null || badElectrodes[electrode] == true)
                    continue;

                //Make it into a DoubleMatrix2D (which is 4x the size of the normal covariance matrix)
                //the cov file only stores half of the matrix (since it's symmetric) and it uses floats.
                //The number of rows in the whitened noise covariance matrix should match the number
                //of rows in the normal covariance matrix stored in the .cov file. If it doesn't
                //the program will throw an exception.

                //Solving the quadratic equation:
                int size = (int) Math.round( (Math.sqrt(1 + 8 * covMatrix.length) - 1) / 2);               
                DenseDoubleMatrix2D doubleCovMatrix = new DenseDoubleMatrix2D(size, size);

                int index = 0;
                for (int i = 0; i < size; i++) {
                    for (int j = i; j < size; j++) {
                        doubleCovMatrix.setQuick(i, j, covMatrix[index]);
                        doubleCovMatrix.setQuick(j, i, covMatrix[index]);
                        index++;
                    }
                }

                //Add zeros back into the whitened noise covariance matrix if necessary.

                DenseDoubleMatrix2D bigWhiteMatrix = new DenseDoubleMatrix2D(size, size);

                int windowLength = nlPoints + nrPoints + 1 - 2;
                int vectorLength = windowLength * adjacentMap[electrode].length;
                int skippedColumns = 0 , skippedRows = 0, k;

                for (int r = 0; r < vectorLength; r++) {
                    k = r/windowLength;
                    if(adjacentMap[electrode][k] == true) {
                        skippedColumns = 0;
                        for (int c = 0; c < vectorLength; c++) {
                            //Figure out which neighbor's values we are putting into our matrix
                            k = c/windowLength;
                            //If it's a good neighbor, copy stuff over
                            if (adjacentMap[electrode][k] == true) {
                                bigWhiteMatrix.setQuick(r, c, reducedWhiteMatrix[electrode].get(r-skippedRows, c-skippedColumns));
                                //If it isn't then fill zeros
                            } else {
                                bigWhiteMatrix.setQuick(r, c, 0);
                                skippedColumns++;
                            }
                        }
                    } else {
                        for(int c=0; c<vectorLength; c++) {
                            bigWhiteMatrix.setQuick(r,c,0);
                        }
                        skippedRows++;
                    }
                }
                


                //Mark this as null so the garbage collector will kill it
                covMatrix = null;

                //Multiply it by the whitened covariance matrices (N * C * N)
                //we use smpBlas.dgemm to do multithreaded multiplication using nThreads (from params.xml)
                DoubleMatrix2D temp = new DenseDoubleMatrix2D(size, size);

                Algebra a = new Algebra();
                temp = a.mult(bigWhiteMatrix, doubleCovMatrix);
                doubleCovMatrix = (DenseDoubleMatrix2D) a.mult(temp, bigWhiteMatrix);

                //Mark this as null so the garbage collector will kill it
                temp = null;

                //doubleCovMatrix is the new covariance matrix that we write out to the .wcov file
                float[] white = new float[(size * (size + 1)) / 2];
                double[][] whiteD = doubleCovMatrix.toArray();

                //We must save the covariance matrix in reduced form
                index = 0;
                for (int i = 0; i < size; i++) {
                    for (int j = i; j < size; j++) {
                        white[index] = (float) whiteD[i][j];
                        index++;
                    }
                }

                whiteD = null;  

                //Serialize the whitened covariance matrix to DATASET.wcov
                //remember that the standard covariance matrix is named DATASET.cov
                whiteCovariance.saveCovarianceMatrix(electrode, white);
            }

            //Close whitened covariance here because the PCA computation
            //will need to read from it.
            whiteCovariance.close();
            app.endProgressBar();

            //Print out a list of bad electrodes
            List<Integer> badElectrodesList = new ArrayList<Integer>();
            for (int i = 1; i < nElectrodes; i++) {
                if (badElectrodes[i]) 
                    badElectrodesList.add(i);
            }
            if (!badElectrodesList.isEmpty()) {
                System.out.print("Electrode(s) ");
                System.out.print(badElectrodesList.toString());
                System.out.print(" appear to be bad (had no spikes or singlular covariance matrices).\n\n");
            }
            
            Vision.getInstance().getCalculationManager().calculationDone();
            
        } catch (Exception e) {
            e.printStackTrace();
            throw new IOException("Error writing projections file!");
        } finally {
            if (covariance != null)      covariance.close();
            if (spikesFile != null)      spikesFile.close();
            if (whiteCovariance != null) whiteCovariance.close();
        }
    }

    public void setParameters(HashMap<String, String> parameters) {
        rawPath = (String) parameters.get("Raw Data Path");
        outputPath = (String) parameters.get("Dataset Folder") + File.separator;
        nThreads = Integer.parseInt((String) parameters.get("Number of Threads"));

        // Tell BLAS (some linear algebra functions) to use multiple CPUs.
        SmpBlas.allocateBlas(nThreads, SeqBlas.seqBlas);

        app = Vision.getInstance();
    }

    DenseDoubleMatrix2D visionCovarianceToColt(float[] covMatrix, int electrode) {
        int size = (int) Math.round( (Math.sqrt(1 + 8 * covMatrix.length) - 1) / 2);               
        DenseDoubleMatrix2D doubleCovMatrix = new DenseDoubleMatrix2D(size, size);

        int index = 0;
        for (int i = 0; i < size; i++) {
            for (int j = i; j < size; j++) {
                doubleCovMatrix.setQuick(i, j, covMatrix[index]);
                doubleCovMatrix.setQuick(j, i, covMatrix[index]);
                index++;
            }
        }
        return doubleCovMatrix;
    }

    public static void displayMatrix(DoubleMatrix2D matrix, String label) {
        PlotUtil.showData(label,
                new DoubleHistogram2D("", 0, matrix.size(), 0, matrix.size(), matrix.toArray()),
                new HistogramStyle());
    }
}
