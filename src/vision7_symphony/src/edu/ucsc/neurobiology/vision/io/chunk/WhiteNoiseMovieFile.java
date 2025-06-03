package edu.ucsc.neurobiology.vision.io.chunk;

import java.io.*;

import edu.ucsc.neurobiology.vision.stimulus.FrameGenerator.*;
import edu.ucsc.neurobiology.vision.stimulus.MovieType;

public class WhiteNoiseMovieFile extends ChunkFile {
    static final boolean DEBUG = false;

    static final long fileID = 24L + ( 248L << 8L) + (3L << 16L) + (46L << 24L) + (171L << 32L) + (159L << 40L) + (13L << 48L) + (2L << 56L);
    
    public WhiteNoiseMovieFile(String fileName) throws IOException {
        super(fileName);
    }
    
    public WhiteNoiseMovieFile(String fileName, int mode) throws IOException {
        super(fileName, mode);
    }

    public long getFileID() { return fileID; }

    public WhiteMovieParams getWhiteMovieParams() throws IOException {
        return new WhiteMovieParams();
    }

    public void setWhiteMovieParams(int width, int height, int seed, ColorType colorType, MovieType movieType,
            RandomNumberGenerator rng, double contrastValue, boolean sparse, float probability) throws IOException {
        new WhiteMovieParams(width, height, seed, colorType, movieType, rng, contrastValue, sparse, probability);
    }
    
    @Override
    public String toString() {
        WhiteMovieParams wmp = null;
        try { wmp = getWhiteMovieParams(); } catch (IOException e) {}
        return wmp.toString();
    }
    
    public class WhiteMovieParams {
        private static final int TAG = 0;
        public int width, height, seed;
        public ColorType colorType;
        public MovieType movieType;
        public RandomNumberGenerator rng;
        public double contrastValue;
        public boolean sparse;
        public float probability;


        //For writing
        WhiteMovieParams(int width, int height, int seed, ColorType colorType, MovieType movieType,
                RandomNumberGenerator rng, double contrastValue, boolean sparse, float probability) throws IOException {
            ByteArrayOutputStream bos = new ByteArrayOutputStream(0);
            DataOutputStream dos = new DataOutputStream(bos);
            dos.writeInt(width);
            dos.writeInt(height);
            dos.writeInt(seed);
            dos.writeInt(colorType.ordinal());
            dos.writeInt(movieType.ordinal());
            dos.writeInt(rng.ordinal());
            dos.writeDouble(contrastValue);
            dos.writeBoolean(sparse);
            dos.writeFloat(probability);
            addChunk(TAG, bos.toByteArray());	

            System.out.println("White Movie Params:");
            System.out.println("width: "		 + width);
            System.out.println("height: "  		 + height);
            System.out.println("seed: "			 + seed);
            System.out.println("colorType: "	 + colorType);
            System.out.println("movieType: "	 + movieType);
            System.out.println("rng: "			 + rng);
            System.out.println("contrastValue: " + contrastValue);
            System.out.println("sparse: "        + sparse);
            System.out.println("probability: "   + probability);
        }

        //For reading
        WhiteMovieParams() throws IOException {
            byte[] buffer = getChunk(TAG);
            DataInputStream dis = new DataInputStream(new ByteArrayInputStream(buffer));
            width = dis.readInt();
            height = dis.readInt();
            seed = dis.readInt();
            colorType = ColorType.values()[dis.readInt()];
            movieType = MovieType.values()[dis.readInt()];
            rng = RandomNumberGenerator.values()[dis.readInt()];
            contrastValue = dis.readDouble();
            
            // If it's an old movie file, these won't be there.  Shouldn't really have to catch this, but nice for some debugging sets where easier not to have to reload movie files
            try {
                sparse = dis.readBoolean();
                probability = dis.readFloat();
            } catch (EOFException e) {};
            
            dis.close();
            
            if (DEBUG) System.out.println(this);
        }
        
        @Override
        public String toString() {
            String nl = System.getProperty("line.separator");
            String out = "White Movie Params:" + nl;
            out += "width: "	 	 + width + nl;
            out += "height: "  	 	 + height + nl;
            out += "seed: "		 	 + seed + nl;
            out += "colorType: "	 + colorType + nl;
            out += "movieType: "	 + movieType + nl;
            out += "rng: "			 + rng + nl;
            out += "contrastValue: " + contrastValue + nl;
            out += "sparse: "        + sparse + nl;
            out += "probability: "   + probability + nl;
            return out;
        }
    }
    
}