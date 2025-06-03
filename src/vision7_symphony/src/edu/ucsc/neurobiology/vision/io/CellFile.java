package edu.ucsc.neurobiology.vision.io;

import java.io.*;

import edu.ucsc.neurobiology.vision.util.*;


/**
 * This class is designed to act as a superclass to all cell-based files. Every cell
 * is accessed by an ID.
 *
 * The class id completely synchronized.
 * TODO make emptyCellIndex have a good value if a file is reopen for addition of extractions
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public abstract class CellFile {
    private final static int SLOT_SIZE = 16;

    public final int nSlots;

    private RandomAccessFile dis;
    private int[] cellID, cellSize;
    private long[] cellLocation;
    private final long directoryLocation;
    private int emptyCellIndex = 0;


    protected CellFile(String fileName, final int nSlots, byte[] userHeader) throws
        IOException {
        this.nSlots = nSlots;

        dis = new RandomAccessFile(fileName, "rw");

//        System.out.println("User Header Length: " + userHeader.length);
        dis.writeInt(userHeader.length);
        dis.write(userHeader);
        directoryLocation = dis.getFilePointer();

        cellID = new int[nSlots];
        cellSize = new int[nSlots];
        cellLocation = new long[nSlots];
        for (int i = 0; i < nSlots; i++) {
            cellID[i] = -1;
            cellSize[i] = -1;
            cellLocation[i] = -1;
        }
        flushDirectory();
    }


    protected CellFile(String fileName) throws IOException {
        dis = new RandomAccessFile(fileName, "rw");

        int headerSize = dis.readInt();
        this.directoryLocation = 4 + headerSize;

        dis.seek(directoryLocation);
        this.nSlots = dis.readInt();
        cellID = new int[nSlots];
        cellSize = new int[nSlots];
        cellLocation = new long[nSlots];
        for (int i = 0; i < nSlots; i++) {
            cellID[i] = dis.readInt();
            cellSize[i] = dis.readInt();
            cellLocation[i] = dis.readLong();
        }

//        System.err.println("headerSize " + headerSize);
//        System.err.println("nSlots: "+ nSlots);
//        for (int i = 0; i < cellID.length; i++) {
//            System.err.println(cellID[i]);
//        }
    }


    synchronized private void flushDirectory() throws IOException {
        dis.seek(directoryLocation);

        dis.writeInt(nSlots);
        for (int i = 0; i < nSlots; i++) {
            dis.writeInt(cellID[i]);
            dis.writeInt(cellSize[i]);
            dis.writeLong(cellLocation[i]);
        }
    }


    synchronized private void flushDirectoryEntry(int slotIndex) throws IOException {
        dis.seek(directoryLocation + SLOT_SIZE * slotIndex);
        dis.writeInt(cellID[slotIndex]);
        dis.writeInt(cellSize[slotIndex]);
        dis.writeLong(cellLocation[slotIndex]);
    }


    synchronized public int[] getIDList() {
        IntegerList idList = new IntegerList();
        for (int i = 0; i < cellID.length; i++) {
            if (cellID[i] != -1) {
                idList.add(cellID[i]);
            }
        }
        return idList.toArray();
    }


    synchronized protected byte[] readUserHeader() throws IOException {
        dis.seek(0);
        int userHeaderSIze = dis.readInt();
        byte[] header = new byte[userHeaderSIze];
        dis.readFully(header);
        return header;
    }


    synchronized private int getIndexOfCell(int id) {
        for (int i = 0; cellID[i] != -1; i++) {
            if (cellID[i] == id) {
                return i;
            }
        }

        return -1;
    }


    synchronized public int[] getExtractionIds() {
        IntegerList list = new IntegerList();
        for (int i = 0; i < nSlots; i++) {
            if (cellID[i] != -1) {
                list.add(cellID[i]);
            }
        }

        return list.toArray();
    }


    synchronized protected final void addCell(int id, byte[] cell) throws IOException {
        if (emptyCellIndex == nSlots) {
            throw new IOException("The capacity of the header is used up");
        }

        cellID[emptyCellIndex] = id;
        cellSize[emptyCellIndex] = cell.length;
        cellLocation[emptyCellIndex] = dis.length();
        flushDirectoryEntry(emptyCellIndex);

        dis.seek(cellLocation[emptyCellIndex]);
        dis.write(cell);

        // increment cell index
        emptyCellIndex++;
    }


    synchronized protected byte[] getCell(int id) throws IOException {
        int index = getIndexOfCell(id);
        if (index == -1) {
            return null;
        }
        byte[] cell = new byte[cellSize[index]];
        dis.seek(cellLocation[index]);
        dis.readFully(cell);
        return cell;
    }


    synchronized public void close() throws IOException {
        
        flushDirectory();
        dis.close();
    }

}
