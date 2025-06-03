package edu.ucsc.neurobiology.vision.electrodemap;

import java.util.HashMap;


/**
 * This creates electrode maps from various standard
 * configurations. 61 arrays have the ID 0. 512 array IDs start at 501. 519 array IDs start
 * at 1500.
 *
 * @author Charles A. Loomis, University of California, Santa Cruz
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @author Matthew I. Grivich, University of California, Santa Cruz
 */
public class ElectrodeMapFactory {

    private static HashMap<Integer, int[]> disconnectedMap = new HashMap<Integer, int[]>();

    static {
        // 61 arrays
        disconnectedMap.put(new Integer(0), new int[] {
            9, 25, 57 /*, 33, 53*/
            });

        // 512 arrays
        disconnectedMap.put(new Integer(501), new int[] {
                //      48, 49, 68, 69, 70, 71, 72, 73, 74, 121, 122,
                //      125, 126, 254, 255, 256
            });

        disconnectedMap.put(new Integer(502), new int[] {
                //      17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27,
                //      28, 55, 258, 273, 324, 510, 453, 470, 469,
                //      468, 467, 466, 452, 196, 408, 43, 44, 5
            });

        disconnectedMap.put(new Integer(503), new int[] {
                //      148, 149, 150, 151, 152, 153, 154, 155, 156,
                //      157, 243, 244, 423, 424, 425, 484, 485
            });

        disconnectedMap.put(new Integer(504), new int[] {
                // FIXME not sure if 64 is disconnected
                //                      64,
            });

        // 519 arrays
        disconnectedMap.put(new Integer(1500), new int[] {
                //      1, 130, 259, 260, 389, 390, 519
            });

        disconnectedMap.put(new Integer(1501), new int[] { // 2011-01: Usually what we call the old 30 um
                //      1, 130, 259, 260, 389, 390, 519
            });

        disconnectedMap.put(new Integer(1504), new int[] { // 2011-01: We've apparently started using this for the new 30 um
                //      1, 130, 259, 260, 389, 390, 519
            });
                
                
        // 120 um 519
        disconnectedMap.put(new Integer(1505), new int[] { // 2011-03: 120 um gets another one
                //      1, 130, 259, 260, 389, 390, 519
            });
        
        disconnectedMap.put(new Integer(1506), new int[] { // 2011-01: What we seem to be calling the 120 um
                //      1, 130, 259, 260, 389, 390, 519
            });
                
        disconnectedMap.put(new Integer(3500), new int[] { // 2011-01: What the 120 um is supposed to be called
                //      1, 130, 259, 260, 389, 390, 519
            });
                
                
        // Tetrode Array
        disconnectedMap.put(new Integer(2500), new int[] {

            });
        
        // Debug single Electrode array
        disconnectedMap.put(new Integer(9999), new int[] {
            
        });
    }


    public static ElectrodeMap getElectrodeMap(int rawArrayID) {
        int arrayID = rawArrayID & 0xFFFF;
        int part = (rawArrayID >> 16) & 0xFF;
        int nParts = rawArrayID >> 24;

        //              System.err.println(rawArrayID);
        //              System.err.println(arrayID);
        //              System.err.println(part);
        //              System.err.println(nParts);

        ElectrodeMap m;
        if (arrayID < 250) {
            m = create64Map(60);  
        } else if (arrayID >= 250 && arrayID < 499) {
            m = create64Map(30);
        } else if (arrayID >= 500 && arrayID < 1500) {
            m = create512Map();
        } //else if (arrayID == 65545) {  //Probably trash.  Will delete if no bugs reported.  mgrivich 07/26/2007
        // m = create128Map();
        //} 
        else if (arrayID >= 1500 && arrayID < 2500){
            m = create519Map(30);
        } else if (arrayID >=3500 && arrayID < 4000) {
            m = create519Map(120);
        } else if (arrayID == 9999) {
            m = create1Map();
        } else {
            m = createTetrodeMap();
        }

        int[] discElectrodes = (int[]) disconnectedMap.get(new Integer(arrayID));
        if (discElectrodes != null) {
            m.setDisconnected(discElectrodes);
        }

        if (nParts == 0 || nParts == 1) {
            return m;
        } else {
            try {
                // if the map allows parts
                return m.getElectrodeMap(part, nParts);
            } catch (Exception e) {
                // FIXME: Should have some warning here...
                return m;
            }
        }
    }


    private static ElectrodeMap create1Map() {
        int nElectrode = 2;
        float[] xCoord = new float[nElectrode];
        float[] yCoord = new float[nElectrode];

        //Electrode 0 position
        xCoord[0] = -5;
        yCoord[0] = -5;

        //electrode 1 positions
        xCoord[1] = 0;
        yCoord[1] = 0;

        return new ElectrodeMap(xCoord, yCoord /*, adjacent, nearN*/);
    }

    private static ElectrodeMap createTetrodeMap() {
        int nElectrodes = 68;
        float[] xCoord = new float[nElectrodes];
        float[] yCoord = new float[nElectrodes];

        //Electrode 0 position
        //              TTL Electrode
        xCoord[0] = 0;
        yCoord[0] = 0;

        //              Ground Level Electrode
        xCoord[65] = 4;
        yCoord[65] = 0;
        //              High Level Electrode
        xCoord[66] = 8;
        yCoord[66] = 0;
        //              Low Level Electrode
        xCoord[67] = 12;
        yCoord[67] = 0;


        //              Grounded Electrodes
        xCoord[1] = 0;
        yCoord[1] = 4;

        xCoord[2] = 4;
        yCoord[2] = 4;

        xCoord[30] = 8;
        yCoord[30] = 4;

        xCoord[58] = 12;
        yCoord[58] = 4;

        xCoord[59] = 16;
        yCoord[59] = 4;

        xCoord[60] = 0;
        yCoord[60] = 7;

        xCoord[61] = 4;
        yCoord[61] = 7;

        xCoord[62] = 8;
        yCoord[62] = 7;

        xCoord[63] = 12;
        yCoord[63] = 7;

        xCoord[64] = 16;
        yCoord[64] = 7;

        //              EG21 Electrode
        xCoord[3] = 0;
        yCoord[3] = 12;
        //              EG22 Electrode
        xCoord[29] = 4;
        yCoord[29] = 12;
        //              EG41 Electrode
        xCoord[31] = 8;
        yCoord[31] = 12;
        //              EG42 Electrode
        xCoord[57] = 12;
        yCoord[57] = 12;

        //              R2 Electrode
        xCoord[28] = 16;
        yCoord[28] = 12;
        //              R4 Electrode
        xCoord[56] = 20;
        yCoord[56] = 12;

        //              Tetrode Electrodes



        //              T4A1
        xCoord[4] = 0;
        yCoord[4] = 16;
        //              T4A2
        xCoord[6] = 2;
        yCoord[6] = 16;
        //              T4A3
        xCoord[8] = 0;
        yCoord[8] = 18;
        //              T4A4
        xCoord[10] = 2;
        yCoord[10] = 18;

        //              T5A5
        xCoord[12] = 6;
        yCoord[12] = 16;
        //              T5A6
        xCoord[14] = 8;
        yCoord[14] = 16;
        //              T5A7
        xCoord[16] = 6;
        yCoord[16] = 18;
        //              T5A8
        xCoord[18] = 8;
        yCoord[18] = 18;

        //              T6A9
        xCoord[20] = 12;
        yCoord[20] = 16;
        //              T6A10
        xCoord[22] = 14;
        yCoord[22] = 16;
        //              T6A11
        xCoord[24] = 12;
        yCoord[24] = 18;
        //              T6A12
        xCoord[26] = 14;
        yCoord[26] = 18;

        //              T7B1
        xCoord[5] = 18;
        yCoord[5] = 16;
        //              T7B2
        xCoord[7] = 20;
        yCoord[7] = 16;
        //              T7B3
        xCoord[9] = 18;
        yCoord[9] = 18;
        //              T7B4
        xCoord[11] = 20;
        yCoord[11] = 18;



        //              T8B5
        xCoord[13] = 0;
        yCoord[13] = 22;
        //              T8B6
        xCoord[15] = 2;
        yCoord[15] = 22;
        //              T8B7
        xCoord[17] = 0;
        yCoord[17] = 24;
        //              T8B8
        xCoord[19] = 2;
        yCoord[19] = 24;

        //              T9B9
        xCoord[21] = 6;
        yCoord[21] = 22;
        //              T9B10
        xCoord[23] = 8;
        yCoord[23] = 22;
        //              T9B11
        xCoord[25] = 6;
        yCoord[25] = 24;
        //              T9B12
        xCoord[27] = 8;
        yCoord[27] = 24;

        //              T16A1
        xCoord[32] = 12;
        yCoord[32] = 22;
        //              T16A2
        xCoord[34] = 14;
        yCoord[34] = 22;
        //              T16A3
        xCoord[36] = 12;
        yCoord[36] = 24;
        //              T16A4
        xCoord[38] = 14;
        yCoord[38] = 24;

        //              T17A5
        xCoord[40] = 18;
        yCoord[40] = 22;
        //              T17A6
        xCoord[42] = 20;
        yCoord[42] = 22;
        //              T17A7
        xCoord[44] = 18;
        yCoord[44] = 24;
        //              T17A8
        xCoord[46] = 20;
        yCoord[46] = 24;




        //              T18A9
        xCoord[48] = 0;
        yCoord[48] = 28;
        //              T18A10
        xCoord[50] = 2;
        yCoord[50] = 28;
        //              T18A11
        xCoord[52] = 0;
        yCoord[52] = 30;
        //              T18A12
        xCoord[54] = 2;
        yCoord[54] = 30;

        //              T19B1
        xCoord[33] = 6;
        yCoord[33] = 28;
        //              T19B2
        xCoord[35] = 8;
        yCoord[35] = 28;
        //              T19B3
        xCoord[37] = 6;
        yCoord[37] = 30;
        //              T19B4
        xCoord[39] = 8;
        yCoord[39] = 30;

        //              T20B5
        xCoord[41] = 12;
        yCoord[41] = 28;
        //              T20B6
        xCoord[43] = 14;
        yCoord[43] = 28;
        //              T20B7
        xCoord[45] = 12;
        yCoord[45] = 30;
        //              T20B8
        xCoord[47] = 14;
        yCoord[47] = 30;

        //              T21B9
        xCoord[49] = 18;
        yCoord[49] = 28;
        //              T21B10
        xCoord[51] = 20;
        yCoord[51] = 28;
        //              T21B11
        xCoord[53] = 18;
        yCoord[53] = 30;
        //              T21B12
        xCoord[55] = 20;
        yCoord[55] = 30;

        for (int i = 0; i < nElectrodes; i++) {
            xCoord[i] *= 30;
            yCoord[i] *= 30;
        }


        return new ElectrodeMap(xCoord, yCoord);
    }


    private static ElectrodeMap create64Map(double pitch) {
        return new Hexagonal61ElectrodeMap(pitch);
    }


    private static ElectrodeMap create512Map() {
        return new Rectangular512ElectrodeMap();
    }


    private static ElectrodeMap create519Map(double pitch) {
        return new Hexagonal519ElectrodeMap(pitch);
    }



    //128 references probably trash.  Will delete if no bugs reported.  mgrivich 7/26/2007
    //  private static ElectrodeMap create128Map() {
    //  int nElectrode = 128 + 1;
    //  float[] xCoord = new float[nElectrode];
    //  float[] yCoord = new float[nElectrode];


    //  //Electrode 0 position
    //  xCoord[0] = -63 * 15;
    //  yCoord[0] = -34 * 15;

    //  ElectrodeMap m = create512Map();
    //  Point2D.Double p = new Point2D.Double();
    //  for (int i = 0; i < e128.length; i++) {
    //  m.getPosition(e128[i], p);
    //  xCoord[i + 1] = (int) p.x;
    //  yCoord[i + 1] = (int) p.y;
    //  }

    //  return new ElectrodeMap(xCoord, yCoord);
    //  }


    //  public static int[] e128 = {
    //  5,
    //  7,
    //  14,
    //  16,
    //  21,
    //  23,
    //  30,
    //  32,
    //  37,
    //  39,
    //  46,
    //  48,
    //  53,
    //  55,
    //  62,
    //  64,
    //  69,
    //  71,
    //  78,
    //  80,
    //  85,
    //  87,
    //  94,
    //  96,
    //  101,
    //  103,
    //  110,
    //  112,
    //  117,
    //  119,
    //  126,
    //  128,
    //  129,
    //  131,
    //  133,
    //  135,
    //  146,
    //  148,
    //  150,
    //  152,
    //  161,
    //  163,
    //  165,
    //  167,
    //  178,
    //  180,
    //  182,
    //  184,
    //  193,
    //  195,
    //  197,
    //  199,
    //  210,
    //  212,
    //  214,
    //  216,
    //  225,
    //  227,
    //  229,
    //  231,
    //  242,
    //  244,
    //  246,
    //  248,
    //  258,
    //  260,
    //  265,
    //  267,
    //  274,
    //  276,
    //  281,
    //  283,
    //  290,
    //  292,
    //  297,
    //  299,
    //  306,
    //  308,
    //  313,
    //  315,
    //  322,
    //  324,
    //  329,
    //  331,
    //  338,
    //  340,
    //  345,
    //  347,
    //  354,
    //  356,
    //  361,
    //  363,
    //  370,
    //  372,
    //  377,
    //  379,
    //  394,
    //  396,
    //  398,
    //  400,
    //  409,
    //  411,
    //  413,
    //  415,
    //  426,
    //  428,
    //  430,
    //  432,
    //  441,
    //  443,
    //  445,
    //  447,
    //  458,
    //  460,
    //  462,
    //  464,
    //  473,
    //  475,
    //  477,
    //  479,
    //  490,
    //  492,
    //  494,
    //  496,
    //  505,
    //  507,
    //  509,
    //  511,
    //  };
}
