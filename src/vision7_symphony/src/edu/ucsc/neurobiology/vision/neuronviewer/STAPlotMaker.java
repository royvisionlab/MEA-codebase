package edu.ucsc.neurobiology.vision.neuronviewer;

import java.io.*;
import java.util.*;
import java.util.List;

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.tree.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.gui.*;
import edu.ucsc.neurobiology.vision.io.chunk.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.stimulus.*;
import edu.ucsc.neurobiology.vision.util.*;
import edu.ucsc.neurobiology.vision.neuronviewer.PlotMaker.KeybindingBuilder;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 * @owner Matthew Grivich, University of California, Santa Cruz
 */
public class STAPlotMaker extends CompoundPlotMaker implements KeybindingBuilder {
    public static int drawThreads;
    public static boolean autoZoom;
    public static double azPadFactor;
    public static double azAspectRatio;
    public static String azHotkey;
    public static boolean profile; // Allow rudimentary profiling without Eclipse
    
    public boolean showSTAGaussian = true;
    public boolean showSTAContour = false;
    public double significance = 3;
    public boolean showSigPixels = false;
    public boolean showFrameSpinner = true;
    public boolean showRFwithEI = false;
    public boolean isSTV = false;
    public boolean showPixelMask = false;
    public double rig = -1.0;
    
    public int colorToUseForContour = 2;
    public double contourCutValue = 3;
    
    int width, height;
    double pixelWidth;
    double pixelHeight;
    double xOffset, yOffset;
    
    GlobalsFile.RunTimeMovieParams runParams;
    
    static JDialog timeCoursesDialog;
    static JPanel timeCoursesFrame;
    static JPanel buttonsFrame;
    static JRadioButton superimposeRadioButton;
    static JRadioButton showAllRadioButton;
    static JRadioButton addAllRadioButton;
    static ButtonGroup group;
    static {
        timeCoursesDialog = new JDialog( (JFrame)null, "Time Courses", false);
        timeCoursesDialog.setSize(400, 200);
        group = new ButtonGroup();
        timeCoursesDialog.setDefaultCloseOperation(JFrame.HIDE_ON_CLOSE);
        timeCoursesFrame = new JPanel(new GridLayout(0, 1));
        buttonsFrame = new JPanel(new FlowLayout());
        showAllRadioButton = new JRadioButton("Show Individually", false);
        addAllRadioButton = new JRadioButton("Sum", false);
        superimposeRadioButton = new JRadioButton("Superimpose", true);
        group.add(showAllRadioButton);
        group.add(addAllRadioButton);
        group.add(superimposeRadioButton);
        
        buttonsFrame.add(showAllRadioButton);
        buttonsFrame.add(addAllRadioButton);
        buttonsFrame.add(superimposeRadioButton);
        timeCoursesDialog.add(buttonsFrame,BorderLayout.NORTH);
        timeCoursesDialog.add(timeCoursesFrame,BorderLayout.CENTER);
    }

    private List<Keybinding> keybindings = new LinkedList<Keybinding>();
    
    public STAPlotMaker() {
        super("STA", CLASS_PLOT);
    }
    
    
    public void initialize(NeuronViewer viewer) {
        super.initialize(viewer);

    }


    public Component makePlot(IntegerList list, int plotType, final TreePath classPath) {
        long startTime = System.currentTimeMillis(); // For profiling if activated; see end of method
        
        if (staFile == null) return new JLabel("No STA available");

        try {
            runParams = viewer.globalsFile.getRunTimeMovieParams();
            pixelWidth = runParams.micronsPerStixelX;
            pixelHeight = runParams.micronsPerStixelY;
        } catch (Exception e){ 	
            return new JLabel("Could not make STA plot");
        }

        width = viewer.staCollection.getWidth();
        height = viewer.staCollection.getHeight();
        xOffset = (640 - (int) (width * pixelWidth)) / 2;
        yOffset = (480 - (int) (height * pixelHeight)) / 2;
        final int neuronID = list.get(0);
        
        ParametricEllipse ellipse = null;
        if (paramsFile.hasParameter("x0")) {
            double x0 = paramsFile.getDoubleCell(neuronID, "x0");
            double y0 = paramsFile.getDoubleCell(neuronID, "y0");
            double sigmaX = paramsFile.getDoubleCell(neuronID, "SigmaX");
            double sigmaY = paramsFile.getDoubleCell(neuronID, "SigmaY");
            double theta = paramsFile.getDoubleCell(neuronID, "Theta");

            //rescaled by the makeSTAPanel
            ellipse = new ParametricEllipse(x0, y0, sigmaX, sigmaY, theta, 1, 1);
        }
        
        
        if (isSTV && viewer.stvCollection == null) {
            System.err.println("Neuron STV not available.  Showing STA.");
            isSTV = false;
        }
        
        PlotPanel plPanel = null;
        try {
            plPanel = makeSTAPanel(
                    isSTV ? viewer.stvCollection.getSTA(neuronID) : staFile.getSTA(neuronID),
                            showSTAContour, colorToUseForContour - 1, contourCutValue,
                            significance,
                            showSigPixels,
                            showFrameSpinner,
                            isSTV,
                            viewer.globalsFile,
                            showSTAGaussian ? ellipse : null);
        } catch (IOException ex) {
            ex.printStackTrace();
        }

        
        // Add autozoom hotkey
        if (ellipse != null) {
            ParametricEllipse scaledEllipse = scaleEllipse(ellipse, runParams);
            AutozoomToggle autozoomToggle = new AutozoomToggle(plPanel, scaledEllipse, azPadFactor, azAspectRatio);
            
            Keybinding autozoomKey = new Keybinding();
            autozoomKey.action = autozoomToggle;
            autozoomKey.handle = "Autozoom Toggle";
            autozoomKey.key = KeyStroke.getKeyStroke(azHotkey);
            keybindings.add(autozoomKey);
            
            // Start panel with autoZoom activated?
            if (autoZoom) autozoomToggle.toggle();
        }
        
        
        // Add overlaid EI (SS)
        if (showRFwithEI && imgFile != null) {
            // This feature is broken
            System.err.println("Overlaying STAs on RFs is not supported. Uncheck that option to hide this warning.");
            
            /*float[][][] img = null;
            try {
                img = imgFile.getImage(neuronID);
            } catch (IOException ex1) {
                ex1.printStackTrace();
            }
            if (img != null) {
                double xPanelSize = (plPanel.getRange()[1] - plPanel.getRange()[0]);
                double yPanelSize = (plPanel.getRange()[3] - plPanel.getRange()[2]);
                double xArraySize = viewer.electrodeMap.getNumberOfElectrodes() == 513 ?
                        1890 : 480;
                double yArraySize = viewer.electrodeMap.getNumberOfElectrodes() == 513 ?
                        900 : 480;
                double deltaX, deltaY, xPixelSize, yPixelSize, xCenter, yCenter;
                deltaX = 0;
                deltaY = 0;
                xCenter = 0;
                yCenter = 0;

                if (viewer.arrayCorners != null) {
                    //  Calculate Array size from the Map
                    for (int np = 1; np < viewer.arrayCorners.length; np++) {
                        deltaX = deltaX + Math.abs(viewer.arrayCorners[np].getX() -
                                viewer.arrayCorners[np - 1].getX()) / 2.;
                        deltaY = deltaY + Math.abs(viewer.arrayCorners[np].getY() -
                                viewer.arrayCorners[np - 1].getY()) / 2.;
                    }

                    //  Calculate Pixel Size
                    xPixelSize = xArraySize / deltaX;
                    yPixelSize = yArraySize / deltaY;

                    //  Calculate Array Center Location
                    for (int np = 0; np < viewer.arrayCorners.length - 1; np++) {
                        xCenter = xCenter + viewer.arrayCorners[np].getX() /
                        (viewer.arrayCorners.length - 1);
                        yCenter = yCenter + viewer.arrayCorners[np].getY() /
                        (viewer.arrayCorners.length - 1);
                    }
                    double xOffset = (640 - (int) (width * pixelWidth)) / 2;
                    double yOffset = (480 - (int) (height * pixelHeight)) / 2;
                    xCenter = xCenter - xOffset;
                    yCenter = yCenter - yOffset;
                } else {
                    xPixelSize = 6;
                    yPixelSize = 6;
                    xCenter = 0.5 * xPanelSize;
                    yCenter = 0.5 * yPanelSize;
                }

                //  Calculate compression factors and ref. point coordinates
                double xCompress = (electrodeMap.getBounds()[1] -
                        electrodeMap.getBounds()[0]) /
                        xPixelSize / xPanelSize;
                double yCompress = (electrodeMap.getBounds()[3] -
                        electrodeMap.getBounds()[2]) /
                        yPixelSize / yPanelSize;

                double xStart;
                if (electrodeMap.getNumberOfElectrodes() > 500) {
                    xStart = (1. - xCompress) / 2. + rig * (0.5 - xCenter / xPanelSize);
                } else {
                    xStart = (1. - xCompress) / 2. - rig * (0.5 - xCenter / xPanelSize);
                }
                double yStart = (1. - yCompress) / 2. + rig * (0.5 - yCenter / yPanelSize);

                System.out.println(xStart + " " + yStart + " " + xCenter + " " + yCenter);

                plPanel.add(new PhysiologicalImagePanel(
                        img, null, 2, electrodeMap, neuronFile.getElectrode(neuronID),
                        new File(viewer.filePathRoot).getParent(), neuronID + ".gif",
                        showRFwithEI, true, rig),
                        new PlaceC(xStart, yStart, 0, 0, xCompress, yCompress));
            } else {
                JLabel label = new JLabel("No EI for this neuron", JLabel.CENTER);
                label.setBackground(Color.white);
                plPanel.add(label, new PlaceC(0, 0, 0, 0, 1, 1));
            }*/
        }
        
        for (int i = 0; i < additionalPlots.size(); i++) {
            plPanel.addData(additionalPlots.get(i), additionalStyles.get(i));
        }

        if (profile) System.out.println("STAPlotMaker#makePlot: " + (System.currentTimeMillis() - startTime) + " ms.");
        return plPanel;
    }

    public static PlotPanel makeSTAPanel(
            final STA sta, boolean showSTAContour, int color,
            double significance, boolean showSignificantPixels, boolean showFrameSpinner,
            boolean isSTV, GlobalsFile globals, ParametricEllipse ...ellipse)throws IOException {
        return makeSTAPanel(
                sta, showSTAContour, color, 3,
                significance, showSignificantPixels, showFrameSpinner,
                isSTV, globals, ellipse);
    }


    public static PlotPanel makeSTAPanel(
            final STA sta, boolean showSTAContour, int color,  double contourCut,
            double significance, boolean showSignificantPixels, boolean showFrameSpinner,
            boolean isSTV, GlobalsFile globals, ParametricEllipse ...ellipse) throws IOException {
        final GlobalsFile.RunTimeMovieParams runParams = globals.getRunTimeMovieParams();

        final PlotPanel staPlotPanel = new PlotPanel("sta", true, false, false, false);
    //	final double stixelWidth = sta.getFrame(0).getStixelWidth();
    //	final double stixelHeight = sta.getFrame(0).getStixelHeight();
        final int width = sta.getWidth();
        final int height = sta.getHeight();

        ScatterPlotStyle significantPixelsStyle = new ScatterPlotStyle(
                "Significant Pixels", SymbolType.SQUARE, 1, Color.red, false, Color.red, 1);
        FunctionStyle gaussianFitStyle = new FunctionStyle("Gaussian Fit",
                Color.black, .5f);
        FunctionStyle simpleContour = new FunctionStyle("Simple Contour", Color.black, .5f);
        FunctionStyle fullContour = new FunctionStyle("Contour", Color.red, .5f);
        
        // add the STA
        staPlotPanel.setAxisVisible(false);
        staPlotPanel.addDataPlotter(new ColorPlotPlotter());
        final SimpleColorPlot scp;
        if (!isSTV) {
            scp = new SimpleColorPlot(SimpleColorPlot.NORMALIZE);
        } else {
            scp = new SimpleColorPlot(SimpleColorPlot.NORMALIZE_ALTERNATE);
        }
        int frameToShow = sta.getMainFrameIndex();
        if (frameToShow == -1) frameToShow = sta.size() - 1;
        scp.setFrame(sta, frameToShow);
        scp.setBounds(runParams.xOffset/1000, 
                (runParams.micronsPerStixelX*runParams.width + runParams.xOffset)/1000,
                runParams.yOffset/1000,
                (runParams.micronsPerStixelY*runParams.height+runParams.yOffset)/1000);

        ColorPlotStyle staStyle = new ColorPlotStyle("STA");
        staPlotPanel.addData(scp, staStyle);
        
        // add the click listener
        staPlotPanel.addClickListener(new ScaledMouseListener() {
            public void clickPerformed(Component source, MouseEvent event, double x,
                    double y) {
                // On Mac w/1-button mouse, command-click sends a right mouse event, but ALSO sends a left mouse event :(
                if (SwingUtilities.isRightMouseButton(event)) {	return; }
                
                if (!timeCoursesDialog.isVisible()) {
                    timeCoursesFrame.removeAll();
                }
                
                double[] r = staPlotPanel.getRange();
                
                int i = (int) ( (x - r[0]) / (runParams.micronsPerStixelX/1000));
                int j = sta.getHeight() - (int) ( (y - r[2]) / (runParams.micronsPerStixelY/1000)) - 1;
                double[][] tc = sta.getTimeCourse(i, j);

                PlotPanel timeCoursesPanel;

                int nPanels = timeCoursesFrame.getComponentCount();
                //superimpose time courses
                if (superimposeRadioButton.isSelected() && nPanels > 0) {
                    timeCoursesPanel = (PlotPanel) timeCoursesFrame.getComponent(
                            nPanels - 1);
                } else {
                    // make a new plot for each time course
                    if (showAllRadioButton.isSelected() || nPanels == 0) {
                        timeCoursesPanel = new PlotPanel();
                        timeCoursesFrame.add(timeCoursesPanel);
                        timeCoursesDialog.setSize(400, 200 * (nPanels + 1));
                    } else {
                        // add time courses to each other
                        timeCoursesPanel = (PlotPanel) timeCoursesFrame.getComponent(
                                nPanels - 1);
                        ArrayList<PlotData> data = timeCoursesPanel.getPlotData();
                        
                        // for each time course in the current plot
                        int vectorCount = 0;
                        for (Iterator<PlotData> iter = data.iterator(); iter.hasNext();) {
                            ScatterPlot lastTc = (ScatterPlot)iter.next();
                            double[] pointAdding = new double[lastTc.getPointSize()];
                            int index = 0;
                            int color = vectorCount%3;
                            for(int currPoint = 0; currPoint < tc[color].length; currPoint++) {    
                                // get a point from the time course we are adding
                                lastTc.getDataPoint(currPoint, pointAdding);
                                // add it to the newly selected pixel's time course
                                tc[color][index] = tc[color][index] + pointAdding[1];
                                index++;
                            }
                            vectorCount = vectorCount + 1;
                        }
                        timeCoursesPanel.removeAllData();
                    }
                }

                for (int c = 0; c < 3; c++) {
                    timeCoursesPanel.addData(
                            new ScatterPlot(tc[c]), PlotUtil.timeCoursesStyle2[c]);
                }
                timeCoursesPanel.autoscale();
                
                timeCoursesDialog.validate();
                timeCoursesDialog.setVisible(true);
            }

            public void enteredPerformed(Component source, MouseEvent event,
                    double x, double y) {		
            }

            public void exitedPerformed(Component source, MouseEvent event,
                    double x, double y) {
            }

            public void pressPerformed(Component source, MouseEvent event,
                    double x, double y) {
            }

            public void releasePerformed(Component source, MouseEvent event,
                    double x, double y) {	
            }
            

        });

        // add the significant pixels
        if (showSignificantPixels) {
            int mainIndex = sta.getMainFrameIndex();
            if (mainIndex != -1) {
                ScatterPlot sp = new ScatterPlot();
                ImageFrame mainFrame = sta.getFrame(mainIndex);
                for (int y = 0; y < height; y++) {
                    for (int x = 0; x < width; x++) {
                        if (mainFrame.getPixelSignificance(x, y) > significance) {
                            sp.add( runParams.xOffset/1000 + (x + 0.5) * (runParams.micronsPerStixelX/1000)
                                    , runParams.yOffset/1000 + (height - y - 0.5) * (runParams.micronsPerStixelY/1000));
                        }
                    }
                }
                staPlotPanel.addData(sp, significantPixelsStyle);
            }
        }

        // add the ellipse
        // Must use built in scaling functions for parametric ellipse.
        // Otherwise, asymmetry between X and Y will bite you.
        // major and minor axes do not point in the direction of x and y.
        ParametricEllipse scaledEllipse = null;
        if (ellipse != null) {
            for (int i = 0; i < ellipse.length; i++) {
                if (ellipse[i] != null) {
                    scaledEllipse = scaleEllipse(ellipse[i], runParams);
                    staPlotPanel.addData(scaledEllipse, gaussianFitStyle);
                }
            }
        }
        
        staPlotPanel.setRange(runParams.getSTARange());
        
        if ((color < 0 || color > 2) && showSTAContour == true) {
            System.err.println("colorToUseForContour is an illegal value.");
            System.err.println("You should use 1 for Red, 2 for Green, or 3 for Blue.");
            showSTAContour = false;
        }
        
        // add the contour
        if (showSTAContour) {
//			int polygonPoints = paramsFile.getDoubleCell(neuronID, "x0");
            staPlotPanel.addDataPlotter(new PolygonPlotter());
            Polygon2DAdapter p2d;
            Polygon2D polygon = sta.getContour(contourCut, color, true);
            Polygon2D[] polygons = polygon.divideIntoSeparatePolygons();

            p2d = new Polygon2DAdapter();
            Polygon2D toshow = polygon.simplifyPolygon();
            toshow.transform(runParams.micronsPerStixelX/1000, 
                    runParams.micronsPerStixelY/1000,
                    runParams.xOffset/1000, runParams.yOffset/1000);
            p2d.setPolygon(toshow);
            staPlotPanel.addData(p2d, simpleContour);

//			double xValsSimple[] = polygon.xPoints;
//			double yValsSimple[] = polygon.yPoints;
//			for(int j=0; j<polygon.nPoints; j++) {
//			System.out.println(xValsSimple[j] + " " +  " " + yValsSimple[j]);
//			}

            for (int i = 0; i < polygons.length; i++) {
                p2d = new Polygon2DAdapter();
                polygons[i].transform(runParams.micronsPerStixelX/1000, 
                        runParams.micronsPerStixelY/1000,
                        runParams.xOffset/1000, runParams.yOffset/1000);
                p2d.setPolygon(polygons[i]);
                staPlotPanel.addData(p2d, fullContour);

            }

        }

        if (showFrameSpinner) {
            final JSpinner frameBox = new JSpinner(new SpinnerNumberModel(frameToShow, 0,
                    sta.size() - 1, 1));
            frameBox.addChangeListener(new ChangeListener() {
                public void stateChanged(ChangeEvent event) {
                    try {
                        scp.setFrame(sta, (Integer) frameBox.getValue());
                    } catch (IOException e) {
                        Vision.reportException(e);
                    }
                    staPlotPanel.replotAllData();
                }
            });
            staPlotPanel.add(frameBox, new PlaceC(1, 0, 1, 0));
        }

        staPlotPanel.addPopupAction(new AbstractAction("Print Data") {
            public void actionPerformed(ActionEvent e) {
                ImageFrame f = sta.getMainFrame();

//				for (int i = 0; i < f.getWidth(); i++) {
//				for (int j = 0; j < f.getHeight(); j++) {
//				System.err.print(f.getPixel(i, j, 1) + "\t");
//				}
//				System.err.println();
//				}

                int[] p = sta.getMainFrameParams();
                ScatterPlot sp = new ScatterPlot();
                for (int i = 0; i < f.getWidth(); i++) {
                    sp.add(i, f.getPixel(i, p[2], 1));
                }
                PlotUtil.showData(sp, new ScatterPlotStyle());
            }
        });

        staPlotPanel.addSelectionAction(new SelectionAction("Print") {
            public void selectionPerformed(JComponent source, Selection selection) {
                SelectionRegion r = selection.getSelection();
                Rectangle2D rect = r.getBounds2D();

                int i1 = (int) (rect.getX() / sta.getStixelWidth());
                int j1 = (int) (rect.getY() / sta.getStixelHeight());
                int i2 = (int) ( (rect.getX() + rect.getWidth()) / sta.getStixelWidth());
                int j2 = (int) ( (rect.getY() + rect.getHeight()) / sta.getStixelHeight());
                System.err.println(i1 + ", " + i2);
                System.err.println(j1 + ", " + j2);

                ImageFrame f = sta.getMainFrame();
                for (int i = i1; i < i2; i++) {
                    for (int j = j1; j < j2; j++) {
                        System.err.print(f.getPixel(i, j, 1) + "\t");
                    }
                    System.err.println();
                }
            }
        });
        
        return staPlotPanel;
    }

    private static ParametricEllipse scaleEllipse(ParametricEllipse ellipse, GlobalsFile.RunTimeMovieParams runParams) {
        return new ParametricEllipse(
            ellipse.getX0()+runParams.xOffset/runParams.micronsPerStixelX,
            ellipse.getY0()+runParams.yOffset/runParams.micronsPerStixelY,
            ellipse.getA(),
            ellipse.getB(), 
            ellipse.getTheta() * -1,
            runParams.micronsPerStixelX/1000,
            runParams.micronsPerStixelY/1000);
    }
    
    // @Override only works on Java >= 6 
    public List<Keybinding> getKeybindings()   { return keybindings;  }
    public void             clearKeybindings() { keybindings.clear(); }

    
    class AutozoomToggle extends AbstractAction {
        private static final long serialVersionUID = 1L;
        private final PlotPanel staPlotPanel;
        private final ParametricEllipse ellipse;
        private final double padFactor;
        private final double aspectRatio;
        
        public AutozoomToggle(PlotPanel staPlotPanel, ParametricEllipse ellipse, Double padFactor, Double aspectRatio) {
            this.staPlotPanel = staPlotPanel;
            this.ellipse = ellipse;			
            this.padFactor = padFactor;
            this.aspectRatio = aspectRatio;
        }
        
        public void toggle() {
            actionPerformed(null);
        }
        
        // @Override only works on Java >= 6
        public void actionPerformed(ActionEvent e) {
            Rectangle panelBounds = staPlotPanel.getBounds();
            double panelAspectRatio = panelBounds.getWidth() / panelBounds.getHeight();
            double targetAspectRatio = Double.isNaN(panelAspectRatio) ? aspectRatio : panelAspectRatio * aspectRatio;
            
            Rectangle2D b = ellipse.getBounds2D(padFactor + 1.0, targetAspectRatio);
            double[] bounds = new double[]{b.getMinX(), b.getMaxX(), b.getMinY(), b.getMaxY()};
            
            double[] curRange = staPlotPanel.getRange();
            boolean alreadyZoomed = Arrays.equals(curRange, bounds);
            if (alreadyZoomed) {
                staPlotPanel.zoomOut();
            } else {
                staPlotPanel.setRange(bounds);
                staPlotPanel.addZoomHistory(curRange);
            }
        }
    }
}