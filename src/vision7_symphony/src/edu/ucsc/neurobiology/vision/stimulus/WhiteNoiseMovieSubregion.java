package edu.ucsc.neurobiology.vision.stimulus;

import java.io.IOException;

import edu.ucsc.neurobiology.vision.io.chunk.ChunkFile;
import edu.ucsc.neurobiology.vision.io.chunk.WhiteNoiseMovieFile;
import edu.ucsc.neurobiology.vision.util.ParallelUtil;
import edu.ucsc.neurobiology.vision.util.ParallelUtil.StartEndIndices;


public class WhiteNoiseMovieSubregion extends WhiteNoiseMovie {

    public WhiteNoiseMovieSubregion(String movieFileName, String globalsFileName, int startPixel, int endPixel) {
        super(movieFileName, globalsFileName, false);
        frameGenerator.offsetSeed(startPixel, new float[0]);
        generateSeeds(); // Generate seeds based on full movie size
        
        int numPixels = endPixel - startPixel + 1;
        resize(1, numPixels); // Resize to subregion
    }
    
    // Use factory methods since arg profiles for different calls overlap
    public static WhiteNoiseMovieSubregion threadWhiteNoiseMovieSubregion(String movieFileName, 
            String globalsFileName, int thisThread, int numThreads) throws IOException {

        WhiteNoiseMovieFile movieFile = new WhiteNoiseMovieFile(movieFileName, ChunkFile.READ);
        WhiteNoiseMovieFile.WhiteMovieParams movieParams = movieFile.getWhiteMovieParams();
        int fullNumPixels = movieParams.width * movieParams.height;
        StartEndIndices myPixels = ParallelUtil.subIndex(thisThread, numThreads, 0, fullNumPixels);
        return new WhiteNoiseMovieSubregion(movieFileName, globalsFileName, myPixels.start, myPixels.end);
    }
}
