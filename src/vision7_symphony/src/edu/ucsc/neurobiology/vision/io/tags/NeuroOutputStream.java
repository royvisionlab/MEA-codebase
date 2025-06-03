package edu.ucsc.neurobiology.vision.io.tags;

import java.io.*;


/**
 * Write out tagged blocks to the neurobiology file.
 *
 * @author Charles A. Loomis, University of California, Santa Cruz
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class NeuroOutputStream
    extends DataOutputStream {

    /**
     * Set of tags that can be used by this Stream
     */
    protected NeuroTagSet tagSet;


    public NeuroOutputStream(OutputStream os) throws IOException {
        super(os);
        this.tagSet = NeuroTagSet.getInstance();
    }


    /**
     * Write an unsigned byte. */
    public void writeUnsignedByte(int ub) throws IOException {

        write(ub);
    }


    /**
     * Write an unsigned short. */
    public void writeUnsignedShort(int s) throws IOException {
        write( (s >>> 8) & 0xFF);
        write(s & 0xFF);
    }


    /**
     * Write an unsigned integer. */
    public void writeUnsignedInt(long i) throws IOException {
        write( (int) ( (i >>> 24) & 0xFF));
        write( (int) ( (i >>> 16) & 0xFF));
        write( (int) ( (i >>> 8) & 0xFF));
        write( (int) (i & 0xFF));
    }


    /**
     * Writes the TagHeader to the output stream. <br>
     * The header is an extension of the SWF format:
     * <br>
     * The first 10 bits are the tagID (for identification).
     * <ul>
     *   <li> If the tag length is in te interval [0..0x3D] ([0..61]) then
     *   the remaining 6 bits contain the length of the tag.
     *
     *   <li> If the tag length it bigger that 0x3D (61) but less or equal than 0xFFFF
     *   (65535) then the 6 bits contain the value 0x3F and there are another 2 bytes of
     *   length information in "unsigned short" format.
     *
     *   <li> If the tag length it bigger that 0xFFFF (65535)
     *   then the 6 bits contain the value 0x3E and there are another 4 bytes of length
     *   information in "unsigned int" format.
     * </ul><br>
     * Therefore the maximim tag length is 4 Gb.
     */
    protected void writeTagHeader(int tagID, long tagLength) throws IOException {
        int th = (tagID << 6);

        if (tagLength < 0x3E) {
            th |= tagLength;
            writeUnsignedShort(th);
        } else if (tagLength <= 0xFFFF) {
            th |= 0x3F;
            writeUnsignedShort(th);
            writeUnsignedShort( (int) tagLength);
        } else {
            th |= 0x3E;
            writeUnsignedShort(th);
            writeUnsignedInt(tagLength);
        }

//    writeShort(tagID);
//    writeInt((int)tagLength);
    }


    public void writeString(String s) throws IOException {
        byte[] buffer = s.getBytes();
        this.writeInt(buffer.length);

        for (int i = 0; i < buffer.length; i++) {
            this.writeByte(buffer[i]);
        }
    }


    /**
     * Write a tag.
     */
    public void writeTag(Tag tag) throws IOException {
        int tagID = tag.getTag();

        if (!tagSet.exists(tagID)) {
            throw new IOException("Undefined tag: " + tagID);
        }

        ByteArrayOutputStream array = new ByteArrayOutputStream();
        NeuroOutputStream nos = new NeuroOutputStream(array);
        tag.write(tagID, nos);
        nos.close();

        int len = array.size();
        writeTagHeader(tagID, len);
        this.write(array.toByteArray());
    }

}
