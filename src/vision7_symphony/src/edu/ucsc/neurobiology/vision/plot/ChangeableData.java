package edu.ucsc.neurobiology.vision.plot;


/**
 * This interface can be implemented by sources of Histograms, ScatterPlots,
 * Functions etc. if the content of the corresponding plots will change with time.
 * If implemented, when the data gets added to the <tt>PlotPanel</tt> the <tt>PlotPanel
 * </tt> will call the <tt>setDrawingControl()</tt> method. The plot class has to save
 * the given <tt>DrawingControl</tt> and every time the data changes use it to repaint
 * the plot.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public interface ChangeableData {
    /**
     * The callback method used to give to the plot class the <tt>DrawingControl</tt>
     * associated with the component that actually displays the plot.
     */
    public void setDrawingControl(DrawingControl drawingControl);
}
