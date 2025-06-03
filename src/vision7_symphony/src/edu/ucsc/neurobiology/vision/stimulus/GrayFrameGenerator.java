package edu.ucsc.neurobiology.vision.stimulus;

/*
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class GrayFrameGenerator
    extends FrameGenerator {

    final int nFrames;

    public String toString() {
        return
            "GRAY: " +
            "nFrames " + nFrames + ", "
            ;
    }


    public GrayFrameGenerator(int nFrames) {
        super(640, 1, RandomNumberGenerator.JAVA_RANDOM, 0);

        this.nFrames = nFrames;
    }


    public int getNumberOfFrames() {
        return nFrames;
    }


    public void nextFrame(Object _frame) {
        if (_frame instanceof int[]) {
            throw new IllegalStateException(
                "Frames are not created as int[] with raw movies.");
        }
        float[] frame = (float[]) _frame;

        for (int x = 0; x < width; x++) {
            frame[3 * x + 0] = frame[3 * x + 1] = frame[3 * x + 2] = 0.5f;
        }
    }


    public Object createFrameBuffer() {
        return new int[width * height * 3];
    }

}
