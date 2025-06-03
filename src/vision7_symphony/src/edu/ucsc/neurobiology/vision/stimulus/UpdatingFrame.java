package edu.ucsc.neurobiology.vision.stimulus;


/**
 * The superclass of all frames that can be updated with other frames. Such frames
 * are used to calculate spike triggered averages (STAs).
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public abstract class UpdatingFrame {
    protected int width, height;


    public UpdatingFrame(int width, int height) {
        this.width = width;
        this.height = height;
    }


    /**
     * Updates (averages) the current frame with a new PACKED_ENCODING frame.
     * @param newPackedFrame Object
     */
    public abstract void update(Object newPackedFrame);


    /**
     * Finishes the averaging process by cleaning and finilizing the calculation.
     */
    public abstract void finish();


    public abstract float[] getColorBuffer();


    public abstract float[] getErrorBuffer();


    public abstract float[] getCovarianceBuffer();
}
