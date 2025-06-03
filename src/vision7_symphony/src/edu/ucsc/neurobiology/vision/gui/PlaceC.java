package edu.ucsc.neurobiology.vision.gui;


/**
 *
 * Dumitru Petrusca, University of California, Santa Cruz
 */
public class PlaceC {
    double x;
    double y;
    double ax;
    double ay;
    double wx;
    double wy;

    public PlaceC(double x, double y, double ax, double ay) {
        this(x, y, ax, ay, -1, -1);
    }


    public PlaceC(double x, double y, double ax, double ay, double wx, double wy) {
        this.x = x;
        this.y = y;
        this.ax = ax;
        this.ay = ay;
        this.wx = wx;
        this.wy = wy;
    }
}
