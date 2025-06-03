package edu.ucsc.neurobiology.vision.io.tags;

import java.util.*;


/**
 * A TagSet which contains all of the standard neurobiology tags.
 *
 * @author Charles A. Loomis, University of California, Santa Cruz<br>
 *         Dumitru Petrusca, University of California, Santa Cruz
 */
public class NeuroTagSet {
    public final static int WHITE_NOISE_MOVIE_TAG = 130;
    public final static int DOUBLE_TAG = 3;
    public final static int DOUBLE_ARRAY_TAG = 4;
    public final static int STRING_TAG = 5;


    /**
     * This holds the individual tags. */
    protected Map<Integer, Tag> tags;

    /** The shared instance for all neurobiology streams. */
    private static NeuroTagSet singleton = null;


    /**
     * Constructor for a new NeuroTagSet. Initally the TagSet contains
     * only the standard neurobiology tags.
     */
    private NeuroTagSet() {
        // Initialize the tag classes.
        tags = new HashMap<Integer, Tag>();

        addTag(new DoubleTag(0));
        addTag(new DoubleArrayTag(new double[1]));
        addTag(new StringTag(""));

//        try {
//            addTag(new WhiteNoiseMovie(1, 1, 0, 0, 1000, 0L, ColorType.DEPENDENT,
//                                       MovieType.BINARY_MOVIE,
//                                       RandomNumberGenerator.MAC_RANDOM,
//                                       .48f, new int[] {0, 0}
//                                       , new long[] {}));
//        } catch (IOException e) {
//            Vision.reportFatalException("Cannot create the NeuroTagSet", e);
//        }
    }


    /**
     * Return the shared instance of the neurobiology tag set. This is
     * a shared instance so adding new tags will be visible by all
     * neurobiology input and output streams.
     */
    public static NeuroTagSet getInstance() {
        if (singleton == null) {
            singleton = new NeuroTagSet();
        }
        return singleton;
    }


    /**
     * Add a new tag to this set.  If the
     * tagID returned is the DEFAULT_TAG,
     * then the default handler is
     * set to the given handler.
     */
    public void addTag(Tag tag) {
        int id = tag.getTag();
        if (id != Tag.DEFAULT_TAG) {
            tags.put(new Integer(id), tag);
        } else {
            System.out.println("id == Tag.DEFAULT_TAG");
        }
    }


    public Tag get(int tagID) {
        return (Tag) tags.get(new Integer(tagID));
    }


    public boolean exists(int tagID) {
        return (tags.get(new Integer(tagID)) != null);
    }


}
