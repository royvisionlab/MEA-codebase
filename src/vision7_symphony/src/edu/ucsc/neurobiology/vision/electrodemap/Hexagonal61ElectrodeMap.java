package edu.ucsc.neurobiology.vision.electrodemap;


/**
 *
 * @author Matthew Grivich, The Salk Institute
 */
public class Hexagonal61ElectrodeMap
extends ElectrodeMap {

    private static int nElectrodes = 65;
    private final static float[] xCoord = new float[nElectrodes];
    private final static float[] yCoord = new float[nElectrodes];

    static {
        yCoord[0] = -8;
        xCoord[0] = -8;
        yCoord[1] = 2;
        xCoord[1] = 3;
        yCoord[2] = 4;
        xCoord[2] = 6;
        yCoord[3] = 4;
        xCoord[3] = 4;
        yCoord[4] = 2;
        xCoord[4] = 1;
        yCoord[5] = 6;
        xCoord[5] = 5;
        yCoord[6] = 4;
        xCoord[6] = 2;
        yCoord[7] = 6;
        xCoord[7] = 3;
        yCoord[8] = 8;
        xCoord[8] = 4;
        yCoord[9] = 8;
        xCoord[9] = -8;
        yCoord[10] = 8;
        xCoord[10] = 2;
        yCoord[11] = 6;
        xCoord[11] = 1;
        yCoord[12] = 4;
        xCoord[12] = 0;
        yCoord[13] = 8;
        xCoord[13] = 0;
        yCoord[14] = 6;
        xCoord[14] = -1;
        yCoord[15] = 2;
        xCoord[15] = -1;
        yCoord[16] = 8;
        xCoord[16] = -2;
        yCoord[17] = 4;
        xCoord[17] = -2;
        yCoord[18] = 6;
        xCoord[18] = -3;
        yCoord[19] = 8;
        xCoord[19] = -4;
        yCoord[20] = 6;
        xCoord[20] = -5;
        yCoord[21] = 4;
        xCoord[21] = -4;
        yCoord[22] = 2;
        xCoord[22] = -3;
        yCoord[23] = 4;
        xCoord[23] = -6;
        yCoord[24] = 2;
        xCoord[24] = -5;
        yCoord[25] = 8;
        xCoord[25] = 8;
        yCoord[26] = 0;
        xCoord[26] = -2;
        yCoord[27] = 2;
        xCoord[27] = -7;
        yCoord[28] = 0;
        xCoord[28] = -4;
        yCoord[29] = 0;
        xCoord[29] = -6;
        yCoord[30] = 0;
        xCoord[30] = -8;
        yCoord[31] = -2;
        xCoord[31] = -7;
        yCoord[32] = -2;
        xCoord[32] = -5;
        yCoord[33] = -2;
        xCoord[33] = -3;
        yCoord[34] = -4;
        xCoord[34] = -6;
        yCoord[35] = -4;
        xCoord[35] = -4;
        yCoord[36] = -2;
        xCoord[36] = -1;
        yCoord[37] = -6;
        xCoord[37] = -5;
        yCoord[38] = -4;
        xCoord[38] = -2;
        yCoord[39] = -6;
        xCoord[39] = -3;
        yCoord[40] = -8;
        xCoord[40] = -4;
        yCoord[41] = 0;
        xCoord[41] = 0;
        yCoord[42] = -8;
        xCoord[42] = -2;
        yCoord[43] = -6;
        xCoord[43] = -1;
        yCoord[44] = -4;
        xCoord[44] = 0;
        yCoord[45] = -8;
        xCoord[45] = 0;
        yCoord[46] = -6;
        xCoord[46] = 1;
        yCoord[47] = -2;
        xCoord[47] = 1;
        yCoord[48] = -8;
        xCoord[48] = 2;
        yCoord[49] = -4;
        xCoord[49] = 2;
        yCoord[50] = -6;
        xCoord[50] = 3;
        yCoord[51] = -8;
        xCoord[51] = 4;
        yCoord[52] = -6;
        xCoord[52] = 5;
        yCoord[53] = -4;
        xCoord[53] = 4;
        yCoord[54] = -2;
        xCoord[54] = 3;
        yCoord[55] = -4;
        xCoord[55] = 6;
        yCoord[56] = -2;
        xCoord[56] = 5;
        yCoord[57] = -8;
        xCoord[57] = 8;
        yCoord[58] = 0;
        xCoord[58] = 2;
        yCoord[59] = -2;
        xCoord[59] = 7;
        yCoord[60] = 0;
        xCoord[60] = 4;
        yCoord[61] = 0;
        xCoord[61] = 6;
        yCoord[62] = 0;
        xCoord[62] = 8;
        yCoord[63] = 2;
        xCoord[63] = 7;
        yCoord[64] = 2;
        xCoord[64] = 5;


    }




    /**
     * This constructor contains all the information necessary for the
     * electrode map.  Generally this is called by some utility routine.
     */
    public Hexagonal61ElectrodeMap(double pitch) {
        super(xCoord, yCoord);
        this.pitch = pitch;
        for(int i=0; i< nElectrodes; i++) {
            xPosition[i]*=pitch/2;
            yPosition[i]*=pitch/2;
        }
        //change to standard handedness.  Connector is already down, as it should be.
        flipAndRotate(true, false, 0, 0, 0);
        xAdjacencyScale = Math.sqrt(2);
        yAdjacencyScale = xAdjacencyScale;
    }




}
