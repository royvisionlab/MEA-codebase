package edu.ucsc.neurobiology.vision.io.tags;

import java.io.*;


/**
 * Read neurobiology tags from a stream.
 *
 * @author Charles A. Loomis, University of California, Santa Cruz
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class NeuroInputStream
    extends DataInputStream {

    public boolean debug = true;

    int tagID;
    int tagLength;


    /**
     * Set of tags that can be used by this Stream
     */
    protected NeuroTagSet tagSet;


    /**
     * Constructor for a new neurobiology input stream.  Initially the
     * tag and action sets are empty.
     */
    public NeuroInputStream(InputStream in) {
        super(in);
        this.tagSet = NeuroTagSet.getInstance();
    }


    /**
     * Read an unsigned integer.
     */
    public long readUnsignedInt() throws IOException {
        long i1 = readUnsignedByte();
        long i2 = readUnsignedByte();
        long i3 = readUnsignedByte();
        long i4 = readUnsignedByte();
        return (i1 << 24) + (i2 << 16) + (i3 << 8) + i4;
    }


    /**
     * Decodes and returns the TagHeader.
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
    protected void readTagHeader() throws IOException {
        int temp = readUnsignedShort();
        if (temp == -1) {
            this.tagID = -1;
            this.tagLength = -1;
            return;
        }

        int ID = temp >> 6;

        long length = temp & 0x3F;
        if (length == 0x3F) {
            length = readUnsignedShort();
        } else if (length == 0x3E) {
            length = readUnsignedInt();
        }

        this.tagID = ID;
        this.tagLength = (int) length;
    }


    public String readString() throws IOException {
        byte[] buffer = new byte[readInt()];

        for (int i = 0; i < buffer.length; i++) {
            buffer[i] = readByte();
        }

        return new String(buffer);
    }


    /**
     * Read a tag. */
    public Tag readTag() throws IOException {
        readTagHeader();
        if (tagID == -1) {
            return null;
        }

        // Look up the proper tag.
        Tag tag = tagSet.get(tagID);

        // set max tag length and read tag
        tag = tag.read(tagID, this, tagLength);

        return tag;
    }
}
