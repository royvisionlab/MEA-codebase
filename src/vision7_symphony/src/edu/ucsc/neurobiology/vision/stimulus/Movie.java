package edu.ucsc.neurobiology.vision.stimulus;

import java.io.*;


/**
 * This class defines a movie as an array of frames.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface Movie extends Cloneable {
    /**
     * Returns the movie description.
     * @return String
     */
    public abstract String getDescription();


    /**
     * Returns the time between frames.
     * @return double
     */
    public abstract double getRefreshTime();


    /**
     * Returns the frame indexed by i. The call may be out of order
     *
     * @param index int
     * @return ImageFrame
     * @throws IOException
     */
    public abstract ImageFrame getFrame(int index) throws IOException;


    /**
     * Returns the number of frames in the movie.
     * @return int
     */
    public abstract int size();


    /**
     * Returns the number of stixes on the x-axis.
     * @return int
     */
    public abstract int getWidth();


    /**
     * Returns the number of stixes on the y-axis.
     * @return int
     */
    public abstract int getHeight();

}
