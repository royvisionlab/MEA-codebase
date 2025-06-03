package edu.ucsc.neurobiology.vision.dataview;

import java.beans.*;
import java.io.*;

import java.awt.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.stimulus.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class MovieView {
    private Movie movie;
    private int frameIndex;
    private PlotPanel movieView;
    private SimpleColorPlot colorFramePlot;
    private ColorPlotStyle style;


    public MovieView(Movie movie) throws IOException {
        this.movie = movie;

        movieView = new PlotPanel();
        movieView.setAxisVisible(false);
        movieView.addDataPlotter(new ColorPlotPlotter());
        movieView.setRange(0, movie.getFrame(0).getStixelWidth() * movie.getWidth(), 0,
                           movie.getFrame(0).getStixelHeight() * movie.getHeight());

        movieView.setLabels("x (\u03BCm)", "y (\u03BCm)");

        colorFramePlot = new SimpleColorPlot(SimpleColorPlot.UNCHANGED);
        colorFramePlot.setFrame(movie, 0);
        style = new ColorPlotStyle("Movie");
        movieView.addData(colorFramePlot, style);

        movieView.addToLegend("Size: " + movie.getWidth() + " x " + movie.getHeight());
        movieView.addToLegend("Interval: " + StringUtil.format(movie.getRefreshTime(), 2));
        movieView.addToLegend("Frames: " + movie.size());

        Vision app = Vision.getInstance();
        app.createFrame(movieView, getController(), null, movie.getDescription());
    }


    private JComponent getController() throws IOException {
        ParametersTable table = new ParametersTable();

        // the STA index parameter
        table.addParameter(
            new IntegerParameter("Frame Index", null, null, frameIndex, 0,
                                 Integer.MAX_VALUE),
            new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent event) {
                frameIndex = ( (IntegerParameter) event.getNewValue()).getValue();
                try {
                    colorFramePlot.setFrame(movie, frameIndex);
                } catch (IOException e) {
                    Vision.reportException(e);
                }
                movieView.replotAllData();
            }
        });

        JPanel p = new JPanel(new BorderLayout());
        p.add(new JScrollPane(table), BorderLayout.CENTER);
        return p;
    }


    public void saveFrames() throws IOException {
        for (int i = 0; i < 120; i++) {
            System.out.println("Saving frame " + i);

            colorFramePlot.setFrame(movie, i);
            movieView.replotAllData();

            try {
                Thread.sleep(250);
            } catch (InterruptedException e) {}

            movieView.saveAsPNG(new File(StringUtil.format(i, 0, 3) + ".png"));
        }
    }
}
