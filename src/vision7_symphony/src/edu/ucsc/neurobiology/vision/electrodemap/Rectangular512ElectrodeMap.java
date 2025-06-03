package edu.ucsc.neurobiology.vision.electrodemap;

import edu.ucsc.neurobiology.vision.plot.*;


/**
 * The implementation of the rectangular 512 map with support for 1, 2, 4, 6, and 8-way cutting.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class Rectangular512ElectrodeMap extends ElectrodeMap {

    private static int nElectrodes = 513;
    private final static float[] xCoord = new float[nElectrodes];
    private final static float[] yCoord = new float[nElectrodes];

    public double micronsPerPixelX, micronsPerPixelY, centerX, centerY, angle;
    public boolean flipX, flipY;

    static {
        //Electrode 0 position
        xCoord[0] = -63;
        yCoord[0] = -34;

        //electrodes 1-129 positions
        for (int i = 0; i < 128; i++) {
            xCoord[i + 1] = (16 - i / 4) * 2 - 1;
            yCoord[i + 1] = - (i % 4) * 8 - 2;
            if ( (i % 8) >= 4) {
                yCoord[i + 1] -= 4;
            }
        }

        //electrodes 257-384 positions
        for (int i = 0; i < 128; i++) {
            xCoord[i + 257] = (i / 4 - 16) * 2 + 1;
            yCoord[i + 257] = - (i % 4) * 8 + 30;
            if ( (i % 8) < 4) {
                yCoord[i + 257] -= 4;
            }
        }

        //electrodes 385-512 positions
        for (int i = 0; i < 128; i++) {
            xCoord[i + 385] = (i % 8) * 4 + 33;
            yCoord[i + 385] = - (i / 8) * 4 + 30;
            if ( (i % 16) < 8) {
                xCoord[i + 385] += 2;
            }
        }

        //electrodes 129-256 positions
        for (int i = 0; i < 128; i++) {
            xCoord[i + 129] = (i % 8) * 4 - 63;
            yCoord[i + 129] = (i / 8) * 4 - 30;
            if ( (i % 16) >= 8) {
                xCoord[i + 129] += 2;
            }
        }


    }




    /**
     * This constructor contains all the information necessary for the
     * electrode map.  Generally this is called by some utility routine.
     */
    public Rectangular512ElectrodeMap() {
        super(xCoord, yCoord);
        super.pitch = 60;
        for (int i = 0; i < nElectrodes; i++) {
            super.xPosition[i] *= (super.pitch/4);
            super.yPosition[i] *= (super.pitch/4);
        }
        
        //No flips or rotations are necessary to put the array in the standard configuration.
        // If rotation changes, must also change adjacency scaling factors!
        xAdjacencyScale = 1;
        yAdjacencyScale = 2 / Math.sqrt(3);
    }


    /**
     *
     * @param part int
     * @param nParts int 1, 2, 4, 6, and 8 parts accepted
     * @return Polygon2D
     */
    public Polygon2D getRegionShape(int part, int nParts) {
        part--;

        final double W = Math.abs(getXPosition(129));  //Width of half the array
        final double H = Math.abs(getYPosition(129));  //Height of half the array
        final double dL = 15; //Half distance between electrodes
        double x1 = 0, x2 = 0, y1 = 0, y2 = 0;

        if (nParts == 1) {
            x1 = -W - dL;
            x2 = W + dL;
            y1 = -H - dL;
            y2 = +H + dL;
        } else if (nParts == 2 || nParts == 4) {
            x1 = -W + 2 * W * (part + 0) / nParts;
            x2 = -W + 2 * W * (part + 1) / nParts;
            y1 = -H - dL;
            y2 = +H + dL;

            if (part == 0) {
                x1 -= dL;
            } else if (part == nParts - 1) {
                x2 += dL;
            }
        } else if (nParts == 6) {
            int n = nParts / 2;
            double[] L = { -W + 0 * 2 * W / n - dL, -W + 1 * 2 * W / n + dL,
                    -W + 2 * 2 * W / n + dL, -W + 3 * 2 * W / n + dL};

            if (part < 3) {
                y1 = -H - dL;
            } else {
                y2 = H + dL;
                part -= 3;
            }

            x1 = L[part + 0];
            x2 = L[part + 1];
        } else if (nParts == 8) {
            x1 = -W + 2 * W * (part%4 + 0) / 4;
            x2 = -W + 2 * W * (part%4 + 1) / 4;
            if(part <4) {
                y1 = -H - dL;
                y2 = 0;
            }
            else {
                y1 = 0;
                y2 = +H +dL;
            }


            if (part == 0 ||part == 4) {
                x1 -= dL;
            } else if (part == 3 || part == 7) {
                x2 += dL;
            }

        }

        else {
            throw new Error("wrong 'nParts' " + nParts);
        }
        Polygon2D p = new Polygon2D();
        p.addPoint(x1, y1);
        p.addPoint(x2, y1);
        p.addPoint(x2, y2);
        p.addPoint(x1, y2);
        p.addPoint(x1, y1);

        return p;
    }

}
