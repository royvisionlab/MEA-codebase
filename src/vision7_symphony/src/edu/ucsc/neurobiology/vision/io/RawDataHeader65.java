package edu.ucsc.neurobiology.vision.io;

import java.io.*;


/**
 * This class implements the 65 style header. It is provided for backward compatibility only.
 *
 * @author Charles A. Loomis, University of California, Santa Cruz
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class RawDataHeader65 extends RawDataHeader {
    public String date;
    public String time;
    
    public final int frequency;
    public final int nElectrodes;
    public final int headerSize;
    
    
    public RawDataHeader65(int nElectrodes, int frequency) throws IOException {
        this.nElectrodes = nElectrodes;
        this.frequency = frequency;
        this.headerSize = getBinaryRepresentation().length;
    }
    
    
    public byte[] getBinaryRepresentation() throws IOException {
        // This is a fake old format header. When acquire reads old format files, it
        // checks to make sure that this header is there.
        // The second number - 202, gives the length of the header in bytes minus 4.
        // If there was a comment, it would go at the end. All the numbers in the middle
        // are not very useful but we have to keep them in case old code like acquire
        // depends on them.
        short header[] = {0, 202, 0, 5, 12332, 12588, 12800, 0, -19200, 0, 768, 0, 6448,
                         8224,
                         8224, 8224, 8224, 8224, 8224, 8224, 8224, 8224, 8224, 8224, 8224,
                         16384, 0, -16384,
                         0, 16800, 0, 1, 16544, 0, 1, 1, 14976, 0, 0, 0, 0, 25, 12576,
                         8224, 8224, 8224, 8224, 8224,
                         8224, 8224, 8224, 8224, 8224, 8224, 8256, 0, 192, 0, 65, -24576,
                         0, 320, -24576, 0, 256,
                         314, -32768, 0, 0, 0, 0, 6450, 8224, 8224, 8224, 8224, 8224,
                         8224, 8224, 8224, 8224, 8224,
                         8224, 8224, 16384, 0, -16384, 0, 16800, 0, 1, 16544, 0, 1, 1,
                         14976, 0, 0, 0, 17948, 16384, 0, 0};

        ByteArrayOutputStream bos = new ByteArrayOutputStream();
        DataOutputStream dos = new DataOutputStream(bos);

        for (int i = 0; i < header.length; i++) {
            dos.writeShort(header[i]);
        }

        return bos.toByteArray();
    }


    public RawDataHeader65(InputStream _stream) throws IOException {
        // First read the number of bytes in the header. We only support files in
        // which the 4-bytes of the header size can be interpreted as a positive int.
        DataInputStream stream = new DataInputStream(_stream);
        int _headerSize = stream.readInt();
        if (_headerSize < 0) {
            throw new IOException("Header too large");
        }

        byte[] buffer = new byte[_headerSize];
        int bytesRemaining = _headerSize;
        while (bytesRemaining > 0) {
            bytesRemaining -= stream.read(
                buffer, _headerSize - bytesRemaining, bytesRemaining);
        }
        stream = null;

        headerSize = _headerSize + 4;

        DataInputStream bufferStream =
            new DataInputStream(new ByteArrayInputStream(buffer));

        // Read the channel list string. The size is given as an unsigned 4-byte integer.
        // (We only support here a range in which the 4-bytes map to a positive int.)
        int stringSize = bufferStream.readInt();
        if (stringSize < 0) {
            throw new IOException("channel list string too large");
        }

        // Make a StringBuffer of the correct size.
        // Read the characters one at a time, and append them to the string.
        StringBuffer sb = new StringBuffer(stringSize);
        for (int i = 0; i < stringSize; i++) {
            char c = (char) bufferStream.readByte();
            sb.append(c);
        }

        // Read the length of the group channel setting array and the number of channels.
        int groupArraySize = bufferStream.readInt();
        int _nElectrodes = bufferStream.readInt();
        if (groupArraySize < 0 || _nElectrodes < 0) {
            throw new IOException();
        }

        frequency = 20000;

        // Look for the tab-delimited date.
        try {
            StringBuffer dtBuffer = new StringBuffer();
            char c = (char) bufferStream.readUnsignedByte();
            while (c != '\t') {
                dtBuffer.append(c);
                c = (char) bufferStream.readUnsignedByte();
            }
            date = bufferStream.toString();
        } catch (EOFException eof) {
            date = "";
        }

        // Look for the tab-delimited time.
        try {
            StringBuffer dtBuffer = new StringBuffer();
            char c = (char) bufferStream.readUnsignedByte();
            while (c != '\t') {
                dtBuffer.append(c);
                c = (char) bufferStream.readUnsignedByte();
            }
            time = bufferStream.toString();
        } catch (EOFException eof) {
            time = "";
        }

        nElectrodes = 65;
    }


    public int getFormat() {
        return FORMAT_16BIT_UNCOMPRESSED;
    }


//    public int getFileType() {
//        return 0x64;
//    }


    public int getSampleSize() {
        return 2 * nElectrodes;
    }


    public int getHeaderSize() {
        return headerSize;
    }


    /**
     * Get the number of channels.
     */
    public int getNumberOfElectrodes() {
        return nElectrodes;
    }


    public long getTime() {
        return 0;
    }


    public int getSamplingFrequency() {
        return frequency;
    }


    public int getNumberOfSamples() {
        return 16000000;
    }


    public String getDate() {
        return date;
    }


    @Override
    public String toString() {
        String s = "";

        s += "\n headerSize = " + headerSize;
        s += "\n nElectrodes = " + nElectrodes;
        s += "\n date = " + date;
        s += "\n time = " + time;

        return s;
    }


    public int getArrayID() {
        throw new Error("Not Supported");
    }


    public String getComment() {
        return "";
    }


    public int getElectrodeMapPart() {
        return 0;
    }


    public int getElectrodeMapPartsCount() {
        return 0;
    }


    public String getDatasetIdentifier() {
        return null;
    }


    public void setFormat(int format) {
        throw new Error("Method not implemented");
    }


    public void setNumberOfSamples(int nSamples) {
        throw new Error("Method not implemented");
    }
}
