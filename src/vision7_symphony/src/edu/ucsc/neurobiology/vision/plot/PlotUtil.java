package edu.ucsc.neurobiology.vision.plot;

import java.io.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.*;


/**
 *
 * @author Petrusca Dumitru, University of California, Santa Cruz
 */
public final class PlotUtil {
//    private final static ArrayList frameList = new ArrayList();

    public static Color[] rgb = {Color.red, Color.green, Color.blue};
    public static String[] rgbNames = {"Red", "Green", "Blue"};
    public static ScatterPlotStyle[] timeCoursesStyle = new ScatterPlotStyle[3];
    public static ScatterPlotStyle[] timeCoursesStyle2 = new ScatterPlotStyle[3];

    public static Window rootWindow = Vision.getInstance().getMainFrame();


    static int defaultWidth = 500;
    static int defaultHeight = 500;
    static int defaultX0 = 400;
    static int defaultY0 = 400;

    private static final Object[][]
        colorList = { {Color.black, "black"}
                    , {Color.green, "green"}
                    , {Color.orange, "orange"}
                    , {Color.blue, "blue"}
                    , {Color.cyan, "cyan"}
                    , {Color.gray, "gray"}
                    , {Color.darkGray, "darkGray"}
                    , {Color.magenta, "magenta"}
                    , {Color.pink, "pink"}
                    , {Color.lightGray, "lightGray"}
                    , {Color.yellow, "yellow"}
                    , {new Color(0.0f, 0.5f, 0.0f), "darkGreen"}
                    , {new Color(0.5f, 1.0f, 0.5f), "lightGreen"}
                    , {new Color(0.5f, 0.0f, 0.0f), "darkRed"},
                    {new Color(1.0f, 0.5f, 0.5f), "lightRed"}, {Color.red, "red"}
    };

    static {
        // create the styles for time courses
        for (int cIndex = 0; cIndex < 3; cIndex++) {
            timeCoursesStyle[cIndex] = new ScatterPlotStyle(
                rgbNames[cIndex] + " Time Course",
                SymbolType.SQUARE, 3, rgb[cIndex], true, rgb[cIndex], 1);
            timeCoursesStyle2[cIndex] = new ScatterPlotStyle(
                rgbNames[cIndex] + " Time Course",
                SymbolType.NONE, 3, rgb[cIndex], true, rgb[cIndex], .25f);
        }
    }


    private PlotUtil() {

    }


    public static Color getColor(int i) {
        return (Color) colorList[i % colorList.length][0];
    }


    public static String getColorName(int i) {
        return (String) colorList[i % colorList.length][1];
    }


    public static double[] getPlotBounds(ScatterPlotData scatterPlot) {
        double xMin = Integer.MAX_VALUE;
        double xMax = Integer.MIN_VALUE;
        double yMin = Integer.MAX_VALUE;
        double yMax = Integer.MIN_VALUE;

        int n = scatterPlot.getPointCount();
        if (n == 0) {
            return new double[] {0, 1, 0, 1};
        }
        double[] p = new double[4];

        for (int i = 0; i < n; i++) {
            scatterPlot.getDataPoint(i, p);

            double x = p[0], y = p[1];
            if (x < xMin) {
                xMin = x;
            }
            if (x > xMax) {
                xMax = x;
            }
            if (y < yMin) {
                yMin = y;
            }
            if (y > yMax) {
                yMax = y;
            }
        }

        return new double[] {xMin, xMax, yMin, yMax};
    }


    public static double[] getPlotBounds(HistogramData histogram) {
        double yMin = histogram.getBin(0);
        double yMax = histogram.getBin(0);

        int nBins = histogram.getBinCount();
        for (int i = 1; i < nBins; i++) {
            double bin = histogram.getBin(i);
            if (bin < yMin) {
                yMin = bin;
            }
            if (bin > yMax) {
                yMax = bin;
            }
        }

        return new double[] {
            histogram.getMin(), histogram.getMax(), yMin, yMax};
    }


    /*
        public static PlotPanel showDataInModalDialog(PlotData data, Object style) {
            PlotPanel p = new PlotPanel();
            p.addData(data, style);
            p.autoscale();
     JDialog f = new JDialog( (JFrame)null, "Plot of: " + data.getDescription(), true);
            f.addWindowListener(new WindowAdapter() {
                public void windowClosing(WindowEvent e) {
//                System.exit(0);
                }
            });
            f.add(p, BorderLayout.CENTER);
            f.setBounds(100, 100, defaultWidth, defaultHeight);
            f.setVisible(true);
            return p;
        }
     */


    public static void showData(String title, PlotPanel p, int w, int h) {
        Rectangle r = new Rectangle();
        if (rootWindow != null) {
            Point point = rootWindow.getLocation();
            Dimension d = rootWindow.getSize();
            r.setBounds(
                point.x + d.width / 2 - w / 2,
                point.y + d.height / 2 - h / 2,
                w, h);
        } else {
            r.setBounds(defaultX0, defaultY0, w, h);
        }

        showData(title, p, r);
    }


    public static void showData(String title, PlotPanel p) {
        showData(title, p, defaultWidth, defaultHeight);
    }


    public static void showData(String title, PlotPanel p, Rectangle r) {
        final JFrame f = new JFrame(title);
        f.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        f.addKeyListener(new KeyAdapter() {
            public void keyTyped(KeyEvent e) {
            }


            public void keyPressed(KeyEvent e) {
                if (e.getKeyCode() == KeyEvent.VK_ESCAPE) {
                    f.setVisible(false);
                    f.dispose();
                }
            }

            public void keyReleased(KeyEvent e) {
            }
        });
//        frameList.add(f);
        f.getContentPane().add(p, BorderLayout.CENTER);
        f.setBounds(r);
        f.setVisible(true);
//        f.setAlwaysOnTop(true);
    }


    /*
        public static void showData(String title, PlotPanel p[], Rectangle r) {
            final JFrame f = new JFrame(title);
            f.addKeyListener(new KeyAdapter() {
                public void keyTyped(KeyEvent e) {
                }


                public void keyPressed(KeyEvent e) {
                    if (e.getKeyCode() == KeyEvent.VK_ESCAPE) {
                        f.setVisible(false);
                        f.dispose();
                    }
                }


                public void keyReleased(KeyEvent e) {
                }
            });
            frameList.add(f);
            f.getContentPane().setLayout(new GridLayout(p.length, 1));
            for (int i = 0; i < p.length; i++) {
                f.add(p[i]);
            }
            f.setBounds(r);
            f.setVisible(true);
        }
     */


    public static PlotPanel showData(
        String title, PlotData data, PlotStyle style, int width, int height, int x0,
        int y0) {

        PlotPanel p = new PlotPanel();
        p.addData(data, style);
        p.autoscale();

        final JFrame f = new JFrame(title);
        f.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
        f.addKeyListener(new KeyAdapter() {
            public void keyTyped(KeyEvent e) {
            }


            public void keyPressed(KeyEvent e) {
                if (e.getKeyCode() == KeyEvent.VK_ESCAPE) {
                    f.setVisible(false);
                    f.dispose();
                }
            }


            public void keyReleased(KeyEvent e) {
            }
        });
        f.getContentPane().add(p, BorderLayout.CENTER);

        f.setBounds(x0, y0, width, height);

        f.setVisible(true);
//        f.setAlwaysOnTop(true);

        return p;
    }


    public static PlotPanel showData(String title, PlotData data) {
        if (data instanceof ScatterPlotData) {
            return showData(title, data, new ScatterPlotStyle());
        } else if (data instanceof HistogramData) {
            return showData(title, data, new HistogramStyle());
        } else if (data instanceof FunctionData) {
            return showData(title, data, new FunctionStyle(title));
        } else {
            throw new IllegalArgumentException("Unknown data " + data.getClass());
        }
    }


    public static PlotPanel showData(String title, PlotData data, PlotStyle style) {
        return showData(title, data, style, defaultWidth, defaultHeight);
    }


    public static PlotPanel showData(PlotData data, PlotStyle style) {
        return showData(data.getDescription(), data, style, defaultWidth, defaultHeight);
    }


    public static PlotPanel showData(
        String title, PlotData data, PlotStyle style, int width, int height) {
        Rectangle r = new Rectangle();
        if (rootWindow != null) {
            Point point = rootWindow.getLocation();
            Dimension d = rootWindow.getSize();

            return showData(title, data, style, width, height,
                            point.x + d.width / 2 - width / 2,
                            point.y + d.height / 2 - height / 2);
        } else {
            return showData(title, data, style, width, height, defaultX0, defaultY0);

        }

    }


    public static PlotPanel showArray(String title, double[] array) {
        ScatterPlot scatter = new ScatterPlot(array);
        return showData(title, scatter, new ScatterPlotStyle(SymbolType.FILLED_SQUARE,
            2, Color.black, true, Color.black, 1));
    }


    public static PlotPanel showArray(String title, float[] array) {
        ScatterPlot scatter = new ScatterPlot(null, array, null, null);
        return showData(title, scatter, new ScatterPlotStyle(SymbolType.FILLED_SQUARE,
            2, Color.black, true, Color.black, 1));
    }
    
    public static PlotPanel showArray(String title, short[] array) {
        ScatterPlot scatter = new ScatterPlot(null, array, null, null);
        return showData(title, scatter, new ScatterPlotStyle(SymbolType.FILLED_SQUARE,
            2, Color.black, true, Color.black, 1));
    }


    public static void saveToPNG(PlotData data, PlotStyle style, int w, int h,
                                 String fileName) throws IOException {

        PlotPanel p = new PlotPanel();
        p.addData(data, style);
        p.autoscale();
        p.setSize(w, h);

        p.saveAsPNG(fileName);
    }


    /**
     * This method returns a color with a high contrast compared to the given one
     * (i.e. easy to see).
     */
    public static Color getContrastColor(Color color) {
        // First calculate the intensity of this color.
        float[] rgb = new float[3];
        rgb = color.getColorComponents(rgb);
        float intensity = 0.3f * rgb[0] + 0.59f * rgb[1] + 0.11f * rgb[2];

        return (intensity < 0.5) ? Color.white : Color.black;
    }
}
