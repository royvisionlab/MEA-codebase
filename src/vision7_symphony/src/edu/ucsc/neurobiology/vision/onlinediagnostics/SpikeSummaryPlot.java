package edu.ucsc.neurobiology.vision.onlinediagnostics;

import java.text.*;

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.electrodemap.*;


/**
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class SpikeSummaryPlot
    extends JPanel {

    int j = 0;
    int nElectrodes;
    int numberOfSpikes[];
    int maxMagnatudes[];
    public static int DEF_MAX_MAG_MAX = 750;
    private int maxMagnatudeMax = DEF_MAX_MAG_MAX; //gives largest magnatude that gets considered when giving cicle sizes

    private int spikesMax = (SpikeSummary.DEF_MAX_FREQ * SpikeSummary.DEF_REFRESH) /
                            1000;

    private int electrodePositions[][];
    private int rows = 15; // rows of electrodes - set default value so I don't get divide by zero
    private int columns = 30; // columns of electrodes - set default value so I don't get divide by zero
    private int centerElectrode = 0; //centerElectrode in detailedSpikeSummaryPlot
    private int mouseCenterElectrode = 0; //centerElectrode can be changed by clicking spikeSummaryPlot.
    //This tells you which electrode has been clicked.
    private int style; // determines whether electrodes are in rectangular or hexagonal format
    private Point2D point; //needed to extract coordinates from electrode map
    private Point2D point2; //needed to extract coordinates from electrode map
    private ElectrodeMap electrodeMap; //gives the coordinates of the electrodes in microns
    private int xSpacing, ySpacing; //distance between electrodes;
    FontMetrics fontMetrics;
    int diameters[];
    Color colors[];
    Ellipse2D.Double ellipse = new Ellipse2D.Double();

    static DecimalFormat formatter;
    static {
        formatter = new DecimalFormat();
        formatter.setMaximumFractionDigits(0);
        formatter.setMinimumIntegerDigits(5);
        formatter.setGroupingUsed(false);
    }


    public SpikeSummaryPlot(ElectrodeMap electrodeMap) {
        this.electrodeMap = electrodeMap;
        this.nElectrodes = electrodeMap.getNumberOfElectrodes();
        electrodePositions = new int[nElectrodes][2];
        numberOfSpikes = new int[nElectrodes];
        maxMagnatudes = new int[nElectrodes];
        diameters = new int[nElectrodes];
        colors = new Color[nElectrodes];

        addMouseListener(new MouseAdapter() {
            public void mousePressed(MouseEvent e) {
                int minimumDistance = Integer.MAX_VALUE;
                int distance;
                int xPos = e.getX();
                int yPos = e.getY();
                for (int i = 0; i < nElectrodes; i++) {
                    if (minimumDistance >
                        (distance = (xPos - electrodePositions[i][0]) *
                         (xPos - electrodePositions[i][0])
                         +
                         (yPos - electrodePositions[i][1]) *
                         (yPos - electrodePositions[i][1]))) {
                        minimumDistance = distance;
                        mouseCenterElectrode = i;
                    }
                }
            }
        });

        addComponentListener(new ComponentAdapter() {
            public void componentResized(ComponentEvent e) {
                setElectrodeMap(style);
            }
        });
    }


    //set electrode map positions

    public void setElectrodeMap(int style) {
        this.style = style;
        Dimension d = this.getSize();

        xSpacing = (int) (d.width  / (columns + .5));
        ySpacing = (int) (d.height / (rows    + .5));

        if (nElectrodes == 513) {
            if (style == SpikeSummary.RECTANGULAR) {
                rows = (int) Math.sqrt(26.0 * nElectrodes / 21.0);
                columns = (int) (21.0 * rows / 26.0);

                while (rows * columns < nElectrodes) {
                    rows++;
                }
            } else {

                rows = 17;
                columns = 33;

            }

            if (style == SpikeSummary.RECTANGULAR) {
                int x = xSpacing / 2;
                int y = ySpacing / 2;

                for (int j = 0; j < rows; j++) {
                    for (int i = 0; i < columns; i++) {
                        if (j * columns + i < nElectrodes) {
                            electrodePositions[j * columns + i][0] = x;
                            electrodePositions[j * columns + i][1] = y;
                            x += xSpacing;
                        }
                    }
                    x = xSpacing / 2;
                    y += ySpacing;
                }
            } else {
                double XMultiplier, YMultiplier;

                point = electrodeMap.getPosition(1, point);
                point2 = electrodeMap.getPosition(9, point2);
                XMultiplier = Math.abs(xSpacing / (point.getX() - point2.getX()));
                point2 = electrodeMap.getPosition(5, point2);
                YMultiplier = Math.abs(ySpacing / (point.getY() - point2.getY()));

                for (int i = 0; i < nElectrodes; i++) {
                    point = electrodeMap.getPosition(i, point);

                    electrodePositions[i][0] = (int) (point.getX() * XMultiplier +
                        (double) xSpacing * 15.75 + xSpacing / 2);
                    electrodePositions[i][1] = (int) ( -point.getY() * YMultiplier +
                        (double) ySpacing * 8.5);

                }
            }

        } // if nElectrodes == 513

        // if nElectrodes == 520
        else if (nElectrodes == 520) {
            if (style == SpikeSummary.RECTANGULAR) {
                rows = (int) Math.sqrt(26.0 * nElectrodes / 21.0);
                columns = (int) (21.0 * rows / 26.0);

                while (rows * columns < nElectrodes) {
                    rows++;
                }
            } else {
                rows    = 53;
                columns = 25;
            }

            if (style == SpikeSummary.RECTANGULAR) {
                int x = xSpacing / 2;
                int y = ySpacing / 2;

                for (int j = 0; j < rows; j++) {
                    for (int i = 0; i < columns; i++) {
                        if (j * columns + i < nElectrodes) {
                            electrodePositions[j * columns + i][0] = x;
                            electrodePositions[j * columns + i][1] = y;
                            x += xSpacing;
                        }
                    }
                    x = xSpacing / 2;
                    y += ySpacing;
                }
            } else {
                double XMultiplier, YMultiplier;
                point  = electrodeMap.getPosition(2, point);
                point2 = electrodeMap.getPosition(3, point2);
                XMultiplier = Math.abs(xSpacing / (point.getX() - point2.getX()));
                YMultiplier = Math.abs(ySpacing / (point.getY() - point2.getY()));

                for (int i = 0; i < nElectrodes; i++) {
                    point = electrodeMap.getPosition(i, point);
                    electrodePositions[i][0] = (int) ( point.getX() * XMultiplier +
                        (double) xSpacing * 13);
                    electrodePositions[i][1] = (int) (-point.getY() * YMultiplier +
                        (double) ySpacing * 29);
                }
            }
        }

        // else nElectrodes == 65
        else if (nElectrodes == 65) {
            if (style == SpikeSummary.RECTANGULAR) {
                rows = (int) Math.sqrt(26.0 * nElectrodes / 21.0);
                columns = (int) (21.0 * rows / 26.0);

                while (rows * columns < nElectrodes) {
                    rows++;
                }
            } else {
                rows = 9;
                columns = 18;
            }

            if (style == SpikeSummary.RECTANGULAR) {
                int x = xSpacing / 2;
                int y = ySpacing / 2;

                for (int j = 0; j < rows; j++) {
                    for (int i = 0; i < columns; i++) {
                        if (j * columns + i < nElectrodes) {
                            electrodePositions[j * columns + i][0] = x;
                            electrodePositions[j * columns + i][1] = y;
                            x += xSpacing;
                        }
                    }
                    x = xSpacing / 2;
                    y += ySpacing;
                }
            } else {
                double XMultiplier, YMultiplier;
                point = electrodeMap.getPosition(1, point);
                point2 = electrodeMap.getPosition(60, point2);
                XMultiplier = xSpacing / (point.getX() - point2.getX());
                point2 = electrodeMap.getPosition(60, point2);
                YMultiplier = ySpacing / (point.getY() - point2.getY());

                for (int i = 0; i < nElectrodes; i++) {
                    point = electrodeMap.getPosition(i, point);
                    electrodePositions[i][0] = (int) (point.getX() * XMultiplier +
                        (double) xSpacing * 9 + xSpacing / 2);
                    electrodePositions[i][1] = (int) ( -point.getY() * YMultiplier +
                        (double) ySpacing * 5);
                }
            }
        } else {
            // This should never happen.
            throw new Error("Unrecognized array size");
        }
    }


    public void paintComponent(Graphics g) {
        Graphics2D g2 = (Graphics2D) g;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                            RenderingHints.VALUE_ANTIALIAS_ON);
        Dimension d = getSize();
        xSpacing = (int) (d.width / (columns + .5));
        ySpacing = (int) (d.height / (rows + .5));
        int maxDiameter = Math.min(xSpacing, ySpacing);

        for (int i = 0; i < nElectrodes; i++) {
            diameters[i] = Math.min(maxDiameter,
                                    maxMagnatudes[i] * maxDiameter / maxMagnatudeMax);

            switch ( (int) Math.min(numberOfSpikes[i] * 6 / spikesMax, 5)) {
                case 0:
                    colors[i] = Color.red;
                    break;
                case 1:
                    colors[i] = Color.orange;
                    break;
                case 2:
                    colors[i] = Color.yellow;
                    break;
                case 3:
                    colors[i] = Color.green;
                    break;
                case 4:
                    colors[i] = Color.blue;
                    break;
                case 5:
                    colors[i] = Color.magenta;
                    break;
            }
        }

        fontMetrics = pickFont(g2, "Filled and Stroked GeneralPath", xSpacing * 4);
        fontMetrics = g2.getFontMetrics();

        // draw background
        g2.setPaint(Color.black);
        g2.fillRect(0, 0, d.width, d.height);

        // draw Ellipse2D.Double
        for (int i = 0; i < nElectrodes; i++) {
            g2.setPaint(colors[i]);
            ellipse.x = electrodePositions[i][0] - diameters[i] / 2;
            ellipse.y = electrodePositions[i][1] - diameters[i] / 2;
            ellipse.width = diameters[i];
            ellipse.height = diameters[i];
            g2.fill(ellipse);
        }

        // draw electrode number
        g2.setPaint(Color.red);
        for (int i = 0; i < nElectrodes; i++) {
            g2.drawString(String.valueOf(i),
                          electrodePositions[i][0] + xSpacing / 3,
                          electrodePositions[i][1] + ySpacing / 2);
        }
    }


    int frameIndex = 0;

    public void refresh(int numberOfSpikes[], int maxMagnatudes[], int newCenterElectrode) {
        for (int i = 0; i < nElectrodes; i++) {
            this.numberOfSpikes[i] = numberOfSpikes[i];
            this.maxMagnatudes[i] = maxMagnatudes[i];
        }

        //if centerElectrode was changed in detailedSpikeSummaryPlot, use it.
        if (centerElectrode != newCenterElectrode) {
            mouseCenterElectrode = centerElectrode = newCenterElectrode;
        } else {
            //else change it with this class.
            centerElectrode = mouseCenterElectrode;
        }

        this.repaint();

        /*
                try {
                    Thread.sleep(500);
                    File f = new File("frame" + formatter.format(frameIndex) + ".png");
//            new PNGExportFileType().exportToFile(file, this, null, null, "");
                    new PNGExportFileType().exportToFile(f, this, null, null, "");
                    frameIndex++;
                } catch (Exception ex) {
                    ex.printStackTrace();
                }
         */
    }


    public int getCenterElectrode() {
        return centerElectrode;
    }


    public void setMaxAmplitude(int newMaxAmplitude) {
        maxMagnatudeMax = newMaxAmplitude;
    }


    public void setSpikesMax(int newSpikesMax) {
        spikesMax = Math.max(newSpikesMax, 1);
    }


    public static FontMetrics pickFont(Graphics2D g2, String longString, int xSpace) {
        final int maxCharHeight = 15;
        final int minFontSize = 1;
        boolean fontFits = false;

        Font font = g2.getFont();
        FontMetrics fontMetrics = g2.getFontMetrics();

        int size = font.getSize();
        String name = font.getName();
        int style = font.getStyle();

        while (!fontFits) {
            if ( (fontMetrics.getHeight() <= maxCharHeight)
                && (fontMetrics.stringWidth(longString) <= xSpace)) {
                fontFits = true;
            } else {
                if (size <= minFontSize) {
                    fontFits = true;
                } else {
                    g2.setFont(font = new Font(name, style, --size));
                    fontMetrics = g2.getFontMetrics();
                }
            }
        }

        return fontMetrics;
    }

}
