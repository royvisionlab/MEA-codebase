package edu.ucsc.neurobiology.vision.stimulus;

import java.io.*;


/*
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class JoiningFrameGenerator
    extends FrameGenerator {

    FrameGenerator[] frameGenerators;
    int frameIndex = -1;
    int currentGeneratorIndex = 0;
    int n = 0;


    public JoiningFrameGenerator(FrameGenerator ...frameGenerators) {
        super(640, 1, RandomNumberGenerator.JAVA_RANDOM, 0);

        this.frameGenerators = frameGenerators;
    }


    public void nextFrame(Object _frame) throws IOException {
        n++;
        frameIndex++;
        if (frameIndex == frameGenerators[currentGeneratorIndex].getNumberOfFrames()) {
            currentGeneratorIndex++;
            frameIndex = 0;
        }

        if (n % 50000 == 0) {
            System.err.println(n / 1000 + "k frames");
        }

        frameGenerators[currentGeneratorIndex].nextFrame(_frame);
    }


    public int getNumberOfFrames() {
        int n = 0;
        for (FrameGenerator g : frameGenerators) {
            n += g.getNumberOfFrames();
        }
        return n;
    }


    public Object createFrameBuffer() {
        return new int[width * height * 3];
    }
}
