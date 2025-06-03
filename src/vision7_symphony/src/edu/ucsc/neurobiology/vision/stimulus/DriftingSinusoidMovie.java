package edu.ucsc.neurobiology.vision.stimulus;


/**
 * Implementation of the drifting sinusoid stimulus.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class DriftingSinusoidMovie
    implements AdvancedMovie {
    
    private int nFrames;
    private int width;
    private int height;
    private double stixelWidth;
    private double stixelHeight;
    private double refreshTime;
    private double backgroundColor;
    private double barColor;
    private double spatialPeriod;
    private double temporalPeriod;
    private double velocity;
    private final double contrast;
    int n1;
    
    
    public DriftingSinusoidMovie(
        int nFrames, int width, int height, double pixelSize, double refreshTime,
        double spatialPeriod, double temporalPeriod, double contrast, int n1) {
        
        this.width = width;
        this.height = height;
//        this.stixelWidth = stixelWidth;
//        this.stixelHeight = stixelHeight;
        this.refreshTime = refreshTime;
//        this.backgroundColor = backgroundColor;
//        this.barColor = barColor;
        this.spatialPeriod = spatialPeriod;
        this.temporalPeriod = temporalPeriod;
        velocity = spatialPeriod / temporalPeriod;
        this.contrast = contrast;

        this.n1 = n1;
    }


    public double getContrast() {
        return contrast;
    }


    public double getVelocity() {
        return velocity;
    }


    public String getDescription() {
        return "Drifting Sinusoid";
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


    @Override
    public int size() {
        return nFrames;
    }


    public double getSpatialPeriod() {
        return spatialPeriod;
    }


    public double getTemporalPeriod() {
        return temporalPeriod;
    }


    public ImageFrame getFrame(int frameIndex) {
        STAFrame f = new STAFrame(width, height, stixelWidth, stixelHeight);
        if (frameIndex <= n1) {
            for (int i = 0; i < width; i++) {
                for (int j = 0; j < height; j++) {
                    double s = Math.sin(2 * Math.PI * (
                        (j + 0.5) / spatialPeriod - frameIndex / temporalPeriod));
                    float c = (float) (0.5 + contrast * s);

                    f.setPixel(i, j, 0, c);
                    f.setPixel(i, j, 1, c);
                    f.setPixel(i, j, 2, c);
                }
            }
        } else {
            for (int i = 0; i < width; i++) {
                for (int j = 0; j < height; j++) {
                    f.setPixel(i, j, 0, 0.5f);
                    f.setPixel(i, j, 1, 0.5f);
                    f.setPixel(i, j, 2, 0.5f);
                }
            }
        }
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

        return s;
    }


    public void setFrameEncoding(FrameEncoding frameEncoding) {
        throw new Error("Method not implemented");
    }


    public Object getEncodedFrame(int frameIndex, Object frame) {
        return getFrame(frameIndex).getBuffer();
    }
}