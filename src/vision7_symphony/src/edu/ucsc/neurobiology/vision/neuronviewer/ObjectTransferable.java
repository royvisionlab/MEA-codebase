package edu.ucsc.neurobiology.vision.neuronviewer;

import java.io.*;

import java.awt.datatransfer.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ObjectTransferable
    implements Transferable {

    private DataFlavor supportedFlavor;
    private Object object;


    public ObjectTransferable(Object object) throws ClassNotFoundException {
        this.object = object;
        supportedFlavor = new DataFlavor(DataFlavor.javaJVMLocalObjectMimeType);
    }


    public DataFlavor[] getTransferDataFlavors() {
        return new DataFlavor[] {
            supportedFlavor};
    }


    public boolean isDataFlavorSupported(DataFlavor flavor) {
        return flavor.equals(supportedFlavor);
    }


    public Object getTransferData(DataFlavor flavor) throws UnsupportedFlavorException,
        IOException {
        return object;
    }


    public Object getTransferData() {
        return object;
    }
}
