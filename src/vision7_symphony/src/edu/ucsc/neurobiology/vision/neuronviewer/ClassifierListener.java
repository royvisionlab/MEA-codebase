
package edu.ucsc.neurobiology.vision.neuronviewer;

import javax.swing.tree.*;


/**
 * The interface defines a callback method to tell the NeuronViewer that a classification happened.\
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface ClassifierListener {

    public void classificationHappened(int[] idList, TreePath classPath);

}
