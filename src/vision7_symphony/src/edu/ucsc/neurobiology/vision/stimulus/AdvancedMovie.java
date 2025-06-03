package edu.ucsc.neurobiology.vision.stimulus;

import java.io.*;


/**
 * An advanced movie is able to encode frames in a special, movie-specific way.
 * With such an encoding it is possible to have custom frame adders that perform much
 * faster that a "naive" frame representation. The white, pink and natural power STAs
 * are calculated using PACKED_ENCODING movie frames.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface AdvancedMovie extends Movie {
    
    public enum FrameEncoding {
        FLOAT_ARRAY_ENCODING, PACKED_ENCODING, INT_ENCODING
    }
    
    public Object getEncodedFrame(int index, Object frame) throws IOException;
    
    public void setFrameEncoding(FrameEncoding frameEncoding);
    
    public int size();
}
