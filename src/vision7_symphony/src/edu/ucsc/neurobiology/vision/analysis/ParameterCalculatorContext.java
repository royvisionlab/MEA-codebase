package edu.ucsc.neurobiology.vision.analysis;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.stimulus.*;
import edu.ucsc.neurobiology.vision.electrodemap.*;


/**
 * Provides a context to ParametersCalculator instances. Specifically the NeuronFile, the
 * current neuron id, the current STA, STV and ei are available. Also a method for accessing
 * parameters already calculated by other ParametersCalculator instances can be accessed
 * by getParameter(String name);
 * 
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ParameterCalculatorContext {
    public final NeuronFile neuronFile;
    public int id;
    public STA currentSTA;
    public STA stv;
    public float[][][] currentEI;
    public final Object[] row;

    private final String[][] paramNamesAndTypes;
    private final STACollection staFile, stvFile;
    private final PhysiologicalImagingFile eiFile;
    public final ElectrodeMap map;
    private final int nParameters;


    public ParameterCalculatorContext(ParametersFile paramsFile, NeuronFile neuronFile,
                                      STACollection staFile, STACollection stvFile, PhysiologicalImagingFile eiFile, ElectrodeMap map) {
        this.neuronFile = neuronFile;
        this.staFile = staFile;
        this.stvFile = stvFile;
        this.eiFile = eiFile;
        this.map = map;

        paramNamesAndTypes = paramsFile.getColumnNamesAndTypes();
        this.nParameters = paramNamesAndTypes.length;

        row = new Object[nParameters];
    }


    /**
     * SHOULD BE CALLED ONLY BY MakeParametersFile.
     *
     * @param id int
     * @throws IOException
     */
    public void update(int id) throws IOException {
        this.id = id;

        // get the STA
        currentSTA = (staFile != null) ? staFile.getSTA(id) : null;
        stv = (stvFile != null) ? stvFile.getSTA(id) : null;
        currentEI = (eiFile != null) ? eiFile.getImage(id) : null;

        // reset the row
        Arrays.fill(row, null);
    }


    /**
     * Call this to get a parameter previously calculataed.
     *
     * @param name String
     * @return Object
     */
    public Object getParameter(String name) {
        int index = -1;
        for (int i = 0; i < nParameters; i++) {
            if (paramNamesAndTypes[i][0].equals(name)) {
                index = i;
            }
        }

        if (index == -1) {
            throw new Error("No such parameter: " + name);
        }

        if (row[index] == null) {
            throw new Error("Parameter " + name + " is not yet calculated");
        }

        return row[index];
    }
}
