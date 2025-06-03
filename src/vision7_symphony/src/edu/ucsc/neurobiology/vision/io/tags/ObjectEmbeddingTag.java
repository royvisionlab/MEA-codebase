package edu.ucsc.neurobiology.vision.io.tags;


/**
 * A tag that can embed a certain object and provide setter and getter methods to access it.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public abstract class ObjectEmbeddingTag
    extends Tag {

    public ObjectEmbeddingTag(int tagID) {
        super(tagID, 0);
    }


    /**
     * Returns the object contained in the tag.
     *
     * @return Object
     */
    public abstract Object getEmbeddedObject();


    /**
     * Stores the given object in the tag.
     *
     * @return Object
     */
    public void setEmbeddedObject(Object obj) {
        throw new IllegalAccessError("The method is not supported");
    }
}
