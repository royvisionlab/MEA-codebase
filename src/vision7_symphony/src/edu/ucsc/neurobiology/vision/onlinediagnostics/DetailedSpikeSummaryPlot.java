package edu.ucsc.neurobiology.vision.onlinediagnostics;

import java.util.*;

import java.awt.*;
import java.awt.geom.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.io.*;


/**
 * @author Matthew Grivich, University of California, Santa Cruz
 */
@SuppressWarnings("serial")
public class DetailedSpikeSummaryPlot
    extends JPanel {
    final int nElectrodes;

    private boolean[] showElectrodes;
    private int electrodePositions[][];
    ArrayList<Spike> spikes = new ArrayList<Spike> ();
    int electrodeLength = 100;
    int electrodeHeight = 30;
    public static int DEF_MAX_AMPLITUDE = 750;
    int maxAmplitude = DEF_MAX_AMPLITUDE; // maxAmplitude of spikes on screen, larger spikes get cropped
    int refreshRate = (int) SpikeSummary.DEF_REFRESH;
    int electrodeSpacing = 10;
    double baseTime = 0.0;
    int centerElectrode = 0; //center electrode in visible space
    JPanel drawingArea;
    Rectangle oldVisibleRectangle = new Rectangle(); //used to determine if the user moved the
    //scroll bars

    private Dimension sizeOfMap; // indicates range of scroll bars

    private int style; // determines whether electrodes are in rectangular or hexagonal format
    private Point2D point; //needed to extract coordinates from electrode map
    private Point2D point2; //needed to extract coordinates from electrode map
    private ElectrodeMap electrodeMap; //gives the coordinates of the electrodes in microns

    final int maxCharHeight = 15;
    final int minFontSize = 1;

    private boolean painting = false;
    private boolean writingSpikeArray = false;
    JScrollPane scroller;
    FontMetrics fontMetrics;


    public DetailedSpikeSummaryPlot(ElectrodeMap electrodeMap) {
        this.electrodeMap = electrodeMap;
        this.nElectrodes = electrodeMap.getNumberOfElectrodes();
        showElectrodes = new boolean[nElectrodes];
        electrodePositions = new int[nElectrodes][2];

        for (int i = 0; i < nElectrodes; i++) {
            showElectrodes[i] = false;
        }

        //Set up the drawing area and the painting.
        drawingArea = new JPanel() {
            protected void paintComponent(Graphics g) {
                super.paintComponent(g);
                Graphics2D g2 = (Graphics2D) g;

                // show electrode if it is within the visibleRectangle
                Rectangle visibleRectangle = this.getVisibleRect();
                for (int i = 0; i < nElectrodes; i++) {
                    showElectrodes[i] =
                        (electrodePositions[i][0] >
                         visibleRectangle.x - electrodeLength &&
                         electrodePositions[i][0] <
                         visibleRectangle.x + visibleRectangle.width &&
                         electrodePositions[i][1] >
                         visibleRectangle.y - electrodeHeight &&
                         electrodePositions[i][1] <
                         visibleRectangle.y + visibleRectangle.height);
                }

                fontMetrics = SpikeSummaryPlot.pickFont(
                    g2, "Filled and Stroked GeneralPath",
                    (electrodeLength + electrodeSpacing) * 4);
                fontMetrics = g2.getFontMetrics();
                g2.setPaint(Color.red);

                for (int i = 0; i < nElectrodes; i++) {
                    if (showElectrodes[i]) {
                        //draw baseline
                        g2.drawLine(
                            electrodePositions[i][0],
                            electrodePositions[i][1],
                            electrodePositions[i][0] + electrodeLength,
                            electrodePositions[i][1]);

                        // draw electrode number
                        g2.drawString(String.valueOf(i),
                                      electrodePositions[i][0],
                                      electrodePositions[i][1] +
                                      (electrodeHeight + electrodeSpacing) / 3);
                    }
                }

                //Wait until spike array is not changing
                //I need this to be able to use the sleep(float ms) method.  I use this
                // to avoid spinning while
                //loops when painting is waiting for the spike ArrayList to change, or
                //the spike ArrayList
                //is waiting for painting to change.
                while (writingSpikeArray == true) {
                    try {
                        Thread.sleep(10);
                    } catch (InterruptedException e) {
                    }
                }

                //Declare that we are painting to prevent spike array from changing
                painting = true;
                //draw all spikes
                for (int i = 0; i < spikes.size(); i++) {
                    int spikeElectrode = spikes.get(i).electrode;
                    if (spikeElectrode == 512) {
                        System.out.println(spikes.get(i));
                    }
                    if (showElectrodes[spikeElectrode]) {
                        double spikeAmplitude = spikes.get(i).amplitude;
                        double spikeTime = spikes.get(i).time;

                        g2.drawLine(
                            (int) (electrodePositions[spikeElectrode][0] +
                                   electrodeLength * (spikeTime - baseTime) /
                                   refreshRate),

                            (int) electrodePositions[spikeElectrode][1],

                            (int) (electrodePositions[spikeElectrode][0] +
                                   (electrodeLength * (spikeTime - baseTime)) /
                                   refreshRate),

                            (int) (electrodePositions[spikeElectrode][1] -
                                   Math.min( (spikeAmplitude * electrodeHeight) /
                                            maxAmplitude, electrodeHeight))
                            );
                    }
                }

                //Done painting
                painting = false;
            }
        }; // paintComponent

        drawingArea.setBackground(Color.black);

        //Put the drawing area in a scroll pane.
        scroller = new JScrollPane(drawingArea);

        //needed for scroller
        setLayout(new BorderLayout());

        add(scroller, BorderLayout.CENTER);
    }


    public boolean isElectrodeVisible(int electrode) {
        return showElectrodes[electrode];
    }


    //Electrode at the center of the visible area
    int getCenterElectrode() {
        return centerElectrode;
    }


    public void setMaxAmplitude(int newMaxAmplitude) {
        maxAmplitude = newMaxAmplitude;
    }


    public void setRefreshRate(int newRefreshRate) {
        refreshRate = newRefreshRate;
    }


    public void setElectrodeMap(int style) {
        this.style = style;
        int xSpacing = electrodeLength + electrodeSpacing;
        int ySpacing = electrodeHeight + electrodeSpacing;
        int columns = 9;
        int rows = nElectrodes / columns;
        while (columns * rows < nElectrodes) {
            rows++;
        }

        drawingArea.setPreferredSize(sizeOfMap);
        drawingArea.revalidate();

        if (style == SpikeSummary.RECTANGULAR) {
            int x = xSpacing / 6;
            int y = ySpacing;

            sizeOfMap = new Dimension(columns * (electrodeLength + electrodeSpacing) + 60, /*60 is for width of scrollbar*/
                                         rows * (electrodeHeight + electrodeSpacing) + 60);
            for (int j = 0; j < rows; j++) {
                for (int i = 0; i < columns; i++) {
                    if (j * columns + i < nElectrodes) {
                        electrodePositions[j * columns + i][0] = x;
                        electrodePositions[j * columns + i][1] = y;
                        x += xSpacing;
                    }
                }
                x = xSpacing / 6;
                y += ySpacing;
            }
        }

        //else use real electrode map (hexagonal)
        else {
            if (nElectrodes == 513) {
                double XMultiplier, YMultiplier;

                sizeOfMap = new Dimension( (electrodeLength + electrodeSpacing) * 32 + 60,
                                           (electrodeHeight + electrodeSpacing) * 17 + 60);

                point = electrodeMap.getPosition(1, point);
                point2 = electrodeMap.getPosition(9, point2);
                XMultiplier = (electrodeLength + electrodeSpacing) /
                              (point.getX() - point2.getX());
                point2 = electrodeMap.getPosition(5, point2);
                YMultiplier = (electrodeHeight + electrodeSpacing) /
                              (point.getY() - point2.getY());

                for (int i = 0; i < nElectrodes; i++) {
                    point = electrodeMap.getPosition(i, point);
                    electrodePositions[i][0] = (int) (point.getX() * XMultiplier
                        + (double) (electrodeLength + electrodeSpacing) * 15.75 +
                        electrodeSpacing);
                    electrodePositions[i][1] = (int) ( -point.getY() * YMultiplier
                        + (double) (electrodeHeight + electrodeSpacing) * 8.5 +
                        electrodeSpacing);

                }
            }

            else if (nElectrodes == 520) {
                double XMultiplier, YMultiplier;

                sizeOfMap = new Dimension( (electrodeLength + electrodeSpacing) * 26 + 60,
                                           (electrodeHeight + electrodeSpacing) * 55 + 60);

                point  = electrodeMap.getPosition(2, point);
                point2 = electrodeMap.getPosition(3, point2);
                XMultiplier = Math.abs( (electrodeLength + electrodeSpacing) / (point.getX() - point2.getX()));
                YMultiplier = Math.abs( (electrodeHeight + electrodeSpacing) / (point.getY() - point2.getY()));

                for (int i = 0; i < nElectrodes; i++) {
                    point = electrodeMap.getPosition(i, point);
                    electrodePositions[i][0] = (int) (  point.getX() * XMultiplier
                        + (double) (electrodeLength + electrodeSpacing) * 12.75 + electrodeSpacing);
                    electrodePositions[i][1] = (int) ( -point.getY() * YMultiplier
                        + (double) (electrodeHeight + electrodeSpacing) * 28    + electrodeSpacing);
                }
            }

            //if nElectrodes == 65
            else {
                double XMultiplier, YMultiplier;

                sizeOfMap = new Dimension( (electrodeLength + electrodeSpacing) *  9 + 60,
                                           (electrodeHeight + electrodeSpacing) * 17 + 60);

                point = electrodeMap.getPosition(1, point);
                point2 = electrodeMap.getPosition(60, point2);
                XMultiplier = (electrodeLength + electrodeSpacing) /
                              (2 * (point.getX() - point2.getX()));
                point2 = electrodeMap.getPosition(60, point2);
                YMultiplier = (electrodeHeight + electrodeSpacing) * 2 /
                              (point.getY() - point2.getY());

                for (int i = 0; i < nElectrodes; i++) {
                    point = electrodeMap.getPosition(i, point);
                    electrodePositions[i][0] = (int) (point.getX() * XMultiplier
                        + (double) (electrodeLength + electrodeSpacing) * 4 +
                        electrodeSpacing);
                    electrodePositions[i][1] = (int) ( -point.getY() * YMultiplier
                        + (double) (electrodeHeight + electrodeSpacing) * 8.5 +
                        electrodeSpacing);

                }
            }
        }
        drawingArea.setPreferredSize(sizeOfMap);
        drawingArea.revalidate();
    }


    //Called when there is new data available
    public void refresh(ArrayList<Spike> spikesIn, double baseTime,
                        int newCenterElectrode) {
        //wait until done painting
        while (painting == true) {
            try {
                Thread.sleep(10);
            } catch (InterruptedException e) {}
        }
        //Declare that we are writing spikes so that we don't start painting
        writingSpikeArray = true;
        spikes.clear();
        spikes.ensureCapacity(spikesIn.size());

        for (int i = 0; i < spikesIn.size(); i++) {
            spikes.add( (Spike) spikesIn.get(i).clone());
        }
        writingSpikeArray = false;

        this.baseTime = baseTime;

        Rectangle visibleRectangle = drawingArea.getVisibleRect();

        //if centerElectrode was changed in spikeSummaryPlot, use it.
        if (centerElectrode != newCenterElectrode) {
            this.centerElectrode = newCenterElectrode;
            //move visible area to match new centerElectrode.
            Rectangle newRectangle = new Rectangle(visibleRectangle);
            newRectangle.x = electrodePositions[centerElectrode][0] -
                             visibleRectangle.width / 2;
            newRectangle.y = electrodePositions[centerElectrode][1] -
                             visibleRectangle.height / 2;
            drawingArea.scrollRectToVisible(newRectangle);
            oldVisibleRectangle = newRectangle;
        }
        //else if visibleRectangle has changed, find new centerElectrode
        else if (oldVisibleRectangle != visibleRectangle) {
            //compute center electrode
            int minimumDistance = Integer.MAX_VALUE;
            int distance;
            int xPos = (visibleRectangle.x + visibleRectangle.width) / 2;
            int yPos = (visibleRectangle.y + visibleRectangle.height) / 2;

            for (int i = 0; i < nElectrodes; i++) {
                if (minimumDistance >
                    (distance = (xPos - electrodePositions[i][0]) *
                     (xPos - electrodePositions[i][0])
                     +
                     (yPos - electrodePositions[i][1]) * (yPos - electrodePositions[i][1]))) {
                    minimumDistance = distance;
                    centerElectrode = i;
                }
            }
            oldVisibleRectangle = visibleRectangle;
        }

        this.repaint();
    }


}
