package edu.ucsc.neurobiology.vision.stimulus;

import java.awt.geom.*;


/**
 * Implementation of the moving bar stimulus.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class MovingBarMovie
    implements AdvancedMovie {

    private int nFrames;
    private int width;
    private int height;
    private double stixelWidth, stixelHeight;
    private double refreshTime;
    private double backgroundColor;
    private double barColor;
    private double length;
    private double thickness;
    private double velocity;
    private double x1;
    private double y1;
    private double x2;
    private double y2;
    private double angle;

    private int framesPerPeriod;
    private double[] xMiddle, yMiddle;




    public MovingBarMovie(
        int nFrames, int width, int height, double pixelSize, double refreshTime,
        double length, double thickness, double x1, double y1, double x2, double y2,
        double velocity, double backgroundColor, double barColor) {

//        super(NeuroTagSet.MOVING_BAR_MOVIE_TAG, 1);
        this.width = width;
        this.height = height;
//        this.stixelWidth = stixelHeight;
//        this.stixelHeight = stixelHeight;
        this.refreshTime = refreshTime;
        this.backgroundColor = backgroundColor;
        this.barColor = barColor;
        this.length = length;
        this.thickness = thickness;
        this.velocity = velocity;
        this.x1 = x1;
        this.y1 = y1;
        this.x2 = x2;
        this.y2 = y2;
        this.angle = Math.atan2(y2 - y1, x2 - x1);

        if (nFrames < 0) {
            this.nFrames = 20 + framesPerPeriod;
        } else {
            this.nFrames = nFrames;
        }

        // generate the frame params
        double distance = Math.sqrt( (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
//        System.out.println("distance " + distance);
        framesPerPeriod = (int) Math.ceil(distance / Math.abs(velocity));
        xMiddle = new double[framesPerPeriod];
        yMiddle = new double[framesPerPeriod];

//        System.out.println(
//            this.x1 + ", " + this.y1 + " --> " + this.x2 + ", " + this.y2 +
//            ": " + framesPerPeriod + "frames");

        double angle = Math.atan2(y2 - y1, x2 - x1);
//        System.out.println("angle " + angle);
        double vx = velocity * Math.cos(angle);
        double vy = velocity * Math.sin(angle);
        double x = x1, y = y1;

        for (int i = 0; i < framesPerPeriod; i++) {
            xMiddle[i] = x;
            yMiddle[i] = y;

            // move the bar
//            x += vx * refreshTime;
//            y += vy * refreshTime;
            x += vx * 1;
            y += vy * 1;
        }
    }


    public double getContrast() {
        return barColor;
    }


    private void fillBackground(STAFrame f) {
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                f.setPixel(x, y, 0, (float) backgroundColor);
                f.setPixel(x, y, 1, (float) backgroundColor);
                f.setPixel(x, y, 2, (float) backgroundColor);
            }
        }
    }


    private void drawBar(double x, double y, STAFrame f) {
        GeneralPath bar = new GeneralPath();
        bar.append(new Rectangle2D.Double(
            x - thickness / 2.0, y - length / 2.0, thickness, length), true);
        AffineTransform t = new AffineTransform();
        t.rotate(angle);
        bar.transform(t);

        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                if (bar.intersects(i, j, 1, 1)) {
                    f.setPixel(i, j, 0, (float) barColor);
                    f.setPixel(i, j, 1, (float) barColor);
                    f.setPixel(i, j, 2, (float) barColor);
                }
            }
        }
    }


    private void drawGaussian(double x0, double y0, STAFrame f) {
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                double dx =
                    (x + 0.5 - x0) * Math.cos(angle) - (y + 0.5 - y0) * Math.sin(angle);
                float c = 0.5f + (float) (
                    barColor * Math.exp( -0.5 * dx * dx / (thickness * thickness)));

                f.setPixel(x, y, 0, c);
                f.setPixel(x, y, 1, c);
                f.setPixel(x, y, 2, c);
            }
        }
    }


    public String getDescription() {
        return "Moving Bar Movie(" + framesPerPeriod + ")";
    }


    public double getRefreshTime() {
        return refreshTime;
    }


    public int getWidth() {
        return width;
    }


    public int getHeight() {
        return height;
    }


    public double getVelocity() {
        return velocity;
    }


    @Override
    public int size() {
        return nFrames;
    }


    public ImageFrame getFrame(int frameIndex) {
        STAFrame f = new STAFrame(width, height, stixelWidth, stixelHeight);
//        fillBackground(f);
        drawGaussian(xMiddle[frameIndex], yMiddle[frameIndex], f);
        
        return f;
    }


    public String toString() {
        String s = "Moving Bar Movie: ";

        s += "\n nFrames: " + nFrames;
        s += "\n width: " + width;
        s += "\n height: " + height;
        s += "\n stixelWidth: " + stixelWidth;
        s += "\n stixelHeight: " + stixelHeight;
        s += "\n refreshTime: " + refreshTime;
        s += "\n backgroundColor: " + backgroundColor;
        s += "\n barColor: " + barColor;
        s += "\n length: " + length;
        s += "\n thickness: " + thickness;
        s += "\n velocity: " + velocity;

        return s;
    }


    public void setFrameEncoding(FrameEncoding frameEncoding) {
        throw new Error("Method not implemented");
    }


    public double getY1() {
        return y1;
    }


    public double getX1() {
        return x1;
    }


    public double getThickness() {
        return thickness;
    }


    public double getAngle() {
        return angle;
    }


    public Object getEncodedFrame(int frameIndex, Object frame) {
        return getFrame(frameIndex).getBuffer();
    }
}
