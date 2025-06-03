package edu.ucsc.neurobiology.vision.neuronviewer;


/**
 *
 * @author Matthew Grivich, The Salk Institute
 */
public class Arrow {
    public double x1, y1, x2, y2; //plot coordinates
    public float lineWidth = 1;
    public double headLength = 5; //pixels
    public double headWidth = 5; //pixels


    
    public Arrow(double x1, double y1, double x2, double y2) {
        this.x1 =x1;
        this.y1 = y1;
        this.x2 = x2;
        this.y2 = y2;
        
    }
    
    public double getAngle() {
        double angle = -Math.atan((y2-y1)/(x2-x1));
        if((x2-x1)<0) {
                angle = Math.PI - angle;
        }
        System.out.println(angle*180/Math.PI);
        return angle;
    }
    

}
