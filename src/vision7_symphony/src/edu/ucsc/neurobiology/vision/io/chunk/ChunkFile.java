package edu.ucsc.neurobiology.vision.io.chunk;

import java.io.*;

import java.util.*;
/**
 * Overall design does not give user access to the order of the chunks and does not allow duplicate tags.
 * 
 * Implementation note: if the file is opened for reading, this is implemented to read the whole buffer into 
 * memory and then close the file handle.  So ChunkFile#close() does nothing if opened for reading.
 */
public abstract class ChunkFile {
    
    private String fileName;
    private byte[] buffer;
    private DataOutputStream dos;
    public static final int READ = 0, WRITE = 1;
    private int mode;
    
    public ChunkFile(String fileName) throws IOException {
        this(fileName, new File(fileName).exists() ? READ : WRITE);
    }
    
    public ChunkFile(String fileName, int mode) throws IOException {
        this.fileName = fileName;
        this.mode = mode;
        if (mode == WRITE) {
            newFile();
        } else {
            ReadFullFile();
        }
    }
    
    void newFile() throws IOException {
        mode = WRITE;
        File file = new File(fileName);
        if (file.exists()) file.delete();
        
        dos = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(fileName)));
        dos.writeLong(getFileID());
    }
    
    public void addChunk(int tag, byte[] data) throws IOException {
        dos.writeInt(tag);
        dos.writeLong(data.length);
        dos.write(data);
    }
    
    public void close() throws IOException {
        if (mode == WRITE) {
            dos.flush();
            dos.close();
        }
    }
    
    void ReadFullFile() throws IOException {
        DataInputStream dis = new DataInputStream(new BufferedInputStream(new FileInputStream(fileName)));
        if (dis.readLong() != getFileID()) {
            dis.close();
            throw new IOException("File " + fileName + " is not of the proper type.");
        }

        buffer = new byte[dis.available()];
        dis.readFully(buffer);
        dis.close();
    }
    
    public byte[] getChunk(int tag) throws IOException {
        if (mode == WRITE) throw new IOException("Attempted to read from file.  Mode is write");

        DataInputStream dis = new DataInputStream(new ByteArrayInputStream(buffer));
        byte[] chunkData = null;  // if not found, chunk data is null
        while (dis.available() > 0) {
            int readTag = dis.readInt();
            int size = dis.readInt();
            if (size != 0) throw new IOException("Chunk sizes of larger than 2 gigabytes are not implemented.");

            size = dis.readInt();
            if (size < 0) throw new IOException("Chunk sizes of larger than 2 gigabytes are not implemented.");

            if (readTag == tag) {					
                chunkData = new byte[(int) size];
                dis.readFully(chunkData);
                break;
            } else {
                dis.skip(size);
            }
        }
    
        return chunkData;
    }
    
    abstract public long getFileID();
    
    public void copyTagsTo(ChunkFile toFile, Collection<Integer> tags) throws IOException {
        for (int tag : tags) {
            byte[] chunk = getChunk(tag);
            if (chunk == null) continue;
            toFile.addChunk(tag, chunk);
        }
    }
    
    public Set<Integer> getTags() throws IOException {
        if (mode == WRITE) throw new IOException("Attempted to read from file.  Mode is write");

        DataInputStream dis = new DataInputStream(new ByteArrayInputStream(buffer));
        Set<Integer> tags = new HashSet<Integer>();
        while (dis.available() > 0) {
            tags.add(dis.readInt());

            // Skip to next tag
            int size = dis.readInt();
            if (size != 0) throw new IOException("Chunk sizes of larger than 2 gigabytes are not implemented.");
            size = dis.readInt();
            if (size < 0) throw new IOException("Chunk sizes of larger than 2 gigabytes are not implemented.");
            dis.skip(size);
        }
        
        return tags;
    }
        
}