package edu.ucsc.neurobiology.vision.neuronviewer;

import java.io.*;
import java.util.*;

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import java.awt.image.*;
import javax.swing.*;
import javax.swing.event.*;

import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.gui.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.plot.gif.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class PhysiologicalImagePanel
    extends JPanel implements MouseListener {

    private ElectrodeMap map;
    private int nElectrodes;
    public int markerElectrode;
    private int timepoint;
    private final double[] amplitude;
    private double minAmplitude, maxAmplitude;
    private float[][] average, error;
    private float[] sigma;
    private double threshold;
    private double rig;
    private boolean isOverlaying;
    private boolean skipDisconnected;
    private PlotPanel inset = new PlotPanel();
    public ScatterPlotStyle style = new ScatterPlotStyle(
        SymbolType.NONE, 0, Color.black, true, Color.black, 1 / 4f);
    private boolean animate = false;
    private JPopupMenu popup;
    private JFileChooser fd;
    private String defaultName;
    private double amplitudeFactor = 1;
    private int[] markElectrodes = null;
    
    private int currentElectrode;
    Ellipse2D.Double ellipse = new Ellipse2D.Double();
    ArrayList<String> backgroundText = new ArrayList<String>();
    ArrayList<Double> backgroundTextX = new ArrayList<Double>();
    ArrayList<Double> backgroundTextY = new ArrayList<Double>();
    ArrayList<Font> backgroundTextFont = new ArrayList<Font>();
    ArrayList<Color> backgroundTextColor = new ArrayList<Color>();
    int t0;
    double scaleBarSize = -1;
    float scaleBarThickness = 1;
    double scaleBarX = -1;
    double scaleBarY = -1;
    ArrayList<Arrow> arrows = new ArrayList<Arrow>();

    float[] axonTemplate, dendriteTemplate, bodyTemplate;
    double ampThreshold = 0;
    double maxChi2 = 0.5;
    
    //config.xml editable parameters
    boolean showAmplitudeSpinner;      // hide amplitude spinner?
    boolean showAnimationDelaySpinner; // hide animation framerate spinner?
    boolean showFrameSlider;           // hide animation frame slider?
    private final JSpinner animationDelaySpinner = new JSpinner(); // SpinnerModel set during Panel construction
    private final JSlider frameSlider = new JSlider();
    boolean negativeAmplitudes; 	// Show negative amplitudes as different color in EI movie playback
    double defaultAmplitude;  		// value to start frame spinner at?
    double minDisplayAmplitude;		// minimum ei amplitude to display

    public void addArrow(Arrow a) {
        arrows.add(a);
    }


//    public void addArrow(int electrode, Arrow a) {
//        a.x = map.getXPosition(electrode);
//        a.y = -map.getYPosition(electrode);
//        arrows.add(a);
//    }


    public synchronized void setScaleBar(double scaleBarSize, float scaleBarThickness,
                                         double scaleBarX, double scaleBarY) {
        if (scaleBarSize < 0) {
            throw new IllegalArgumentException("scaleBar cannot be negative");
        }
        this.scaleBarSize = scaleBarSize;
        this.scaleBarThickness = scaleBarThickness;
        this.scaleBarX = scaleBarX;
        this.scaleBarY = scaleBarY;
    }


    public void setAmplitudeFactor(double amplitudeFactor) {
        this.amplitudeFactor = amplitudeFactor;

        for (int electrode = 0; electrode < nElectrodes; electrode++) {
            amplitude[electrode] = MathUtil.maxAbs(average[electrode]);
        }

        for (int electrode = 0; electrode < nElectrodes; electrode++) {
            if (amplitude[electrode] < minDisplayAmplitude) {
                amplitude[electrode] = 0;
            }

            amplitude[electrode] *= amplitudeFactor;
            if (amplitude[electrode] > maxAmplitude) {
                amplitude[electrode] = maxAmplitude;
            }
        }
    }

    public PhysiologicalImagePanel(float[][][] image, float[] sigma, double threshold,
                                   ElectrodeMap map, int markerElectrode) {
        this(image, sigma, threshold, map, markerElectrode, null, null, false, false, 1., true, false, true, 40, true, 1, 0);
    }


    public PhysiologicalImagePanel(float[][][] image, float[] sigma, double threshold,
                                   ElectrodeMap map, int markerElectrode,
                                   String directory, String defaultName) {
        this(image, sigma, threshold, map, markerElectrode, directory, defaultName, false, false,
             1., true, false, true, 40, true, 1, 0);
    }
    
    public PhysiologicalImagePanel(float[][][] image, float[] sigma, double threshold,
            ElectrodeMap map, int markerElectrode,
            String directory, String defaultName,
            boolean isOverlaying, boolean skipDisconnected) {
        
        this(image, sigma, threshold, map, markerElectrode, directory, defaultName,
                isOverlaying, skipDisconnected, 1., true, false, true, 40, true, 1, 0);
    }


    public PhysiologicalImagePanel(float[][][] image, float[] sigma, double threshold,
                                   ElectrodeMap map, int markerElectrode,
                                   String directory, String defaultName,
                                   boolean isOverlaying, boolean skipDisconnected, double rig) {
        this(image, sigma, threshold, map, markerElectrode, directory, defaultName,
             isOverlaying, skipDisconnected, rig, true, false, true, 40, true, 1, 0);
    }


    public PhysiologicalImagePanel(float[][][] image, float[] sigma, double threshold,
                                   ElectrodeMap map, int markerElectrode,
                                   String directory, String defaultName,
                                   boolean isOverlaying, boolean skipDisconnected,
                                   double rig, boolean showAmplitudeSpinner, boolean showAnimationDelaySpinner, boolean showFrameSlider, int defaultAnimationDelay, boolean negativeAmplitudes, double defaultAmplitude, double minDisplayAmplitude) {

        super(new PlaceLayout());
        this.setName("ei");
        this.setBackground(Color.white);

        this.minDisplayAmplitude = minDisplayAmplitude;
        this.showAmplitudeSpinner      = showAmplitudeSpinner;
        this.showAnimationDelaySpinner = showAnimationDelaySpinner;
        this.showFrameSlider           = showFrameSlider;
        this.negativeAmplitudes = negativeAmplitudes;
        this.defaultAmplitude = defaultAmplitude;
        this.map = map;
        this.nElectrodes = map.getNumberOfElectrodes();
        this.markerElectrode = markerElectrode;
        average = image[0];
        error = image[1];
        this.sigma = sigma;
        this.threshold = threshold;
        this.isOverlaying = isOverlaying;
        this.skipDisconnected = skipDisconnected;
        this.rig = rig;
        addMouseListener(this);

//      Define if the inset is transparent (SS)
        if (isOverlaying) {
            this.setOpaque(false);
        }

        //////////////////////////////////////////
        amplitude = new double[nElectrodes];

        for (int electrode = 0; electrode < nElectrodes; electrode++) {
            amplitude[electrode] = MathUtil.maxAbs(average[electrode]);
        }
        minAmplitude = MathUtil.min(amplitude);
        maxAmplitude = MathUtil.max(amplitude);
 
        int maxAmplitudeElectrode = MathUtil.maxIndex(amplitude);
        t0 = MathUtil.maxAbsIndex(average[maxAmplitudeElectrode]);

        setAmplitudeFactor(defaultAmplitude);

        this.add(inset, new PlaceC(1, 0, 1, 0, 0.5, 0.5));
        inset.setVisible(false);
        style.setConnectingPoints(true);
        style.setSymbolType(SymbolType.NONE);
        inset.setLabels("Time (ms)", "Amplitude (ADC)");
        inset.addClickListener(new ScaledMouseListener() {
            public void clickPerformed(Component source, MouseEvent event, double x, double y) { hideInset(); }
            public void pressPerformed(Component source, MouseEvent event, double x, double y)   {};
            public void releasePerformed(Component source, MouseEvent event, double x, double y) {};
            public void enteredPerformed(Component source, MouseEvent event, double x, double y) {};
            public void exitedPerformed(Component source, MouseEvent event, double x, double y)  {};
        });
        inset.getYAxis().setFixedTickSpacing(0, 1000, 2); // 1000 is replaced after autoscaling
        
        initPopup();

        fd = new JFileChooser();
        if (directory != null) {
            fd.setCurrentDirectory(new File(directory));
        }
        this.defaultName = defaultName;

        animationDelaySpinner.setModel(new SpinnerNumberModel(defaultAnimationDelay, 1, 1000, 5));
        animationDelaySpinner.setToolTipText("Animation frame delay");
        frameSlider.setModel(new DefaultBoundedRangeModel(0, 0, 0, average[0].length-1));
        frameSlider.setToolTipText("Frame selector");
        frameSlider.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                if (frameSlider.getValueIsAdjusting()) {
                    pauseAnimation();
                }
                setFrame((Integer) frameSlider.getValue());
                repaint();
            }
        });
        addControls();
    }


    private void addControls() {
        if (showAmplitudeSpinner) {
            // add the amplitude spinner
            final JSpinner spinner = new JSpinner(new SpinnerNumberModel(
                    (int) amplitudeFactor, 1, 1000000, 1));
            PhysiologicalImagePanel.this.add(spinner, new PlaceC(0, 0, 0, 0, -1, -1));
            spinner.addChangeListener(new ChangeListener() {
                public void stateChanged(ChangeEvent e) {
                    //                System.err.println("ev " +(Integer) spinner.getValue());
                    setAmplitudeFactor( (Integer) spinner.getValue());
                    PhysiologicalImagePanel.this.repaint();
                }
            });
            spinner.setToolTipText("Amplitude display multiplier");
        }

        if (showAnimationDelaySpinner) {
            // add the animation delay spinner at bottom left
            PhysiologicalImagePanel.this.add(animationDelaySpinner, new PlaceC(0, 1, 0, 1));
        }
        
        if (showFrameSlider) {
            // add the animation frame slider at bottom right
            PhysiologicalImagePanel.this.add(frameSlider, new PlaceC(1, 1, 1, 1));
        }
            
        this.validate();
    }


    private void removeControls() {
        int componentCount = getComponentCount();
        for (int i = componentCount - 1; i > -1; i--) {
            Component component = getComponent(i);
            if (component instanceof JSpinner || component instanceof JSlider) {
                remove(component);
            }
        }

        this.validate();
    }

    // Fields used by popup
    PlotPanel lastPanel = null;
    AbstractAction addTrace = new AbstractAction("Add Trace") {
        public void actionPerformed(ActionEvent event) {
            lastPanel.addData(getPlotFor(currentElectrode), style);
            lastPanel.autoscale();
            lastPanel.padY();
        }
    };

    
    private void initPopup() {
        popup = new JPopupMenu();

        JMenuItem item = new JMenuItem("Play/Stop");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                synchronized (PhysiologicalImagePanel.this) {
                    if (animate) {
                        endAnimation();
                    } else {
                        startAnimation();
                    }
                }
            }
        });
        popup.add(item);

        item = new JMenuItem("Save Animated GIF");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (defaultName != null) {
                    fd.setSelectedFile(new File(defaultName + ".gif"));
                }
                if (fd.showSaveDialog(PhysiologicalImagePanel.this) ==
                    JFileChooser.APPROVE_OPTION) {

                    makeGIF(fd.getSelectedFile().getAbsolutePath());
                }
            }
        });
        popup.add(item);

        item = new JMenuItem("Save Image");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                removeControls();
                GraphicsIO.saveImage(PhysiologicalImagePanel.this,
                                     PhysiologicalImagePanel.this, "Save Image");
                addControls();
            }
        });
        popup.add(item);

        item = new JMenuItem("Save QuickTime Movie");
        item.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                // save all frames as JPG
                ArrayList<String> list = new ArrayList<String>();

                int w = 500;
                int h = 250;
                removeControls();
                for (int t = 0; t < average[0].length; t++) {
                    setFrame(t);
                    try {
                        String name =
                            defaultName + "-" + StringUtil.format(t, 0, 3) + ".jpg";
                        GraphicsIO.saveComponentToFile(
                            PhysiologicalImagePanel.this, name, "JPG", w, h);
                        list.add(name);

//                        saveAsEPS(defaultName + "-" + StringUtil.format(t, 0, 3), 7, 3.5);
                        try {
                            GraphicsIO.saveComponentToFile(
                                PhysiologicalImagePanel.this,
                                defaultName + "-" + StringUtil.format(t, 0, 3) + ".png",
                                "PNG", 800, 400);
                        } catch (IOException ex1) {
                            ex1.printStackTrace();
                        }
                    } catch (IOException ex) {
                        ex.printStackTrace();
                    }
                }

                JpegImagesToMovie.save("movie.mov", w, h, 20, list);

                addControls();

                // restore
                for (int electrode = 0; electrode < nElectrodes; electrode++) {
                    amplitude[electrode] = MathUtil.maxAbs(average[electrode]);
                }
                repaint();
            }
        });
        popup.add(item);

        popup.add(new AbstractAction("Copy to Clipboard") {
            public void actionPerformed(ActionEvent event) {
                Clipboard.saveComponentToClipboard(PhysiologicalImagePanel.this);
            }
        });

        popup.add(new AbstractAction("Show Trace") {
            public void actionPerformed(ActionEvent event) {
                lastPanel = new PlotPanel();
                lastPanel.setLabels("Time (ms)", "Amplitude (ADC)");
                lastPanel.addData(getPlotFor(currentElectrode), style);
                lastPanel.autoscale();
                lastPanel.padY();
                PlotUtil.showData("" + currentElectrode, lastPanel);
                
                addTrace.setEnabled(true);
            }
        });
        
        addTrace.setEnabled(false);
        popup.add(addTrace);
    }


    public synchronized void setMarkElectrodes(int[] markElectrodes) {
        Arrays.sort(markElectrodes);
        this.markElectrodes = markElectrodes;
    }


    public void drawColored(double x, double y, double rx, double ry, Color color,
                            Graphics2D g) {
        g.setColor(color);
        ellipse.setFrame(x - rx, y - ry, 2 * rx, 2 * ry);
        g.fill(ellipse);
    }


    public void addBackgroundText(String text, double x, double y, Font font) {
        addBackgroundText(text, x, y, font, Color.black);
    }


    public void addBackgroundText(String text, double x, double y, Font font, Color color) {
        backgroundText.add(text);
        backgroundTextX.add(x);
        backgroundTextY.add(y);
        backgroundTextFont.add(font);
        backgroundTextColor.add(color);
    }


    public void removeAllBackgroundText() {
        backgroundText.clear();
        backgroundTextX.clear();
        backgroundTextY.clear();
        backgroundTextFont.clear();
        backgroundTextColor.clear();
    }


    public void mouseClicked(MouseEvent e) {

    }


    public void mouseEntered(MouseEvent e) {

    }


    public void mouseExited(MouseEvent e) {

    }


    public int getClosestElectrode(int x0, int y0) {
        double[] bounds = map.getBounds();

        double x2 = getWidth();
        double y2 = getHeight();

        double x1p = bounds[0];
        double x2p = bounds[1];
        double y1p = bounds[2];
        double y2p = bounds[3];
      
        double ax = x2 *.9 / (x2p - x1p); //pixels/ electrode map unit
        double ay = y2 *.9 / (y2p - y1p);  //pixels/ electrode map unit
        double bx = -x2*x1p/(x2p-x1p); //pixels offset (50%)
        double by = -y2*y1p/(y2p-y1p); //pixels offset (50%)

        double minDistance = Double.POSITIVE_INFINITY;
        int electrode = -1;

        for (int e = 0; e < nElectrodes; e++) {

              
            double x = (ax * map.getXPosition(e) + bx);        
            double y = (-ay * map.getYPosition(e) + by);
            
            double d = Math.sqrt( (x - x0) * (x - x0) + (y - y0) * (y - y0));
            if (d < minDistance) {
                minDistance = d;
                electrode = e;
            }
        }

        return electrode;
    }


    boolean showing = false;
    public void hideInset() {
        showing = false;
        inset.setVisible(false);
    }
    
    public void mousePressed(MouseEvent event) {
        currentElectrode = getClosestElectrode(event.getX(), event.getY());

        if (SwingUtilities.isRightMouseButton(event)) {
            if (lastPanel == null || !lastPanel.isValid()) {
                addTrace.setEnabled(false);
            }
            popup.show(PhysiologicalImagePanel.this, event.getX(), event.getY());
        } else if (SwingUtilities.isLeftMouseButton(event)) {
            if (showing) {
                hideInset();
                return;
            } else {
                // Didn't seem like this was used, given the inset.removeAllData() directly below?
//                if (sigma != null) {
//                    ScatterPlot sp = new ScatterPlot();
//                    for (int i = 0; i < average[currentElectode].length; i++) {
//                        sp.add(i, -sigma[currentElectode] * threshold);
//                    }
//                    inset.addData(sp, style);
//                }

                inset.removeAllData();
                inset.addData(getPlotFor(currentElectrode), style);
                inset.autoscale();
                inset.getYAxis().setNumRoundedFixedLinearTicks(4);
                inset.removeAllBackgroundText();
                inset.addBackgroundText("" + currentElectrode, PlotPanel.RIGHT,
                                        PlotPanel.BOTTOM, Font.decode("Arial PLAIN 20"));

                inset.setVisible(true);
                showing = true;
            }
        }
    }


    public ScatterPlot getPlotFor(int electrode) {
        ScatterPlot sp = new ScatterPlot();
        for (int i = 0; i < average[electrode].length; i++) {
            sp.add( (i - t0) / 20.0, average[electrode][i]);
        }
        return sp;
    }


    public void mouseReleased(MouseEvent e) {
//        inset.setVisible(false);
    }


    public void showInAWindow(String title) {
        JFrame fr = new JFrame(title);
        fr.getContentPane().add(this);
        fr.setBounds(100, 100, 800, 400);
        fr.setVisible(true);
    }


    public void setFrame(int t) {
        timepoint = t;
        for (int electrode = 0; electrode < nElectrodes; electrode++) {
            amplitude[electrode] = MathUtil.maxAbs(average[electrode]);
            if (amplitude[electrode] < minDisplayAmplitude) {
                amplitude[electrode] = 0;
            } else {
                amplitude[electrode] = Math.abs(average[electrode][t] * amplitudeFactor);
                if (amplitude[electrode] > maxAmplitude) {
                    amplitude[electrode] = maxAmplitude;
                }
            }
        }
    }


    synchronized public void startAnimation() {
        animate = true;
        new Thread() {
            public void run() {
                boolean _animate;
                for (int t = timepoint; t < average[0].length; t++) {
                    synchronized (PhysiologicalImagePanel.this) {
                        _animate = animate;
                    }
                    
                    if (! (_animate && isShowing())) {
                        // finish the animation
                        return;
                    } else {
                        frameSlider.setValue(t);
                        try {
                            Thread.sleep( (Integer) animationDelaySpinner.getValue());
                        } catch (InterruptedException ex) {}
                    }
                }
                while (true) {
                    for (int t = 0; t < average[0].length; t++) {
                        synchronized (PhysiologicalImagePanel.this) {
                            _animate = animate;
                        }
                        
                        if (! (_animate && isShowing())) {
                            // finish the animation
                            return;
                        } else {
                            frameSlider.setValue(t);
                            try {
                                Thread.sleep( (Integer) animationDelaySpinner.getValue());
                            } catch (InterruptedException ex) {}
                        }
                    }
                }
            }
        }.start();
    }
    
    synchronized public void pauseAnimation() {
        animate = false;
    }
    
    synchronized public void endAnimation() {
        setAmplitudeFactor(amplitudeFactor);
        repaint();
        animate = false;
    }


    public void makeGIF(String name) {
        AnimatedGifEncoder e = new AnimatedGifEncoder();
        e.start(name);
        e.setDelay(20); // 1 frame per sec

        for (int t = 0; t < average[0].length; t++) {
            for (int electrode = 0; electrode < nElectrodes; electrode++) {
                amplitude[electrode] = Math.abs(average[electrode][t] * amplitudeFactor);
            }

            BufferedImage image = new BufferedImage(
                getWidth(), getHeight(), BufferedImage.TYPE_BYTE_INDEXED);
            this.paintComponent(image.getGraphics());

            e.addFrame(image);
        }

        // restore
        for (int electrode = 0; electrode < nElectrodes; electrode++) {
            amplitude[electrode] = MathUtil.maxAbs(average[electrode]);
        }
        repaint();

        e.finish();
    }


    final synchronized public void paintComponent(Graphics _g) {
        super.paintComponent(_g);

        Graphics2D g = (Graphics2D) _g;
        g.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                           RenderingHints.VALUE_ANTIALIAS_ON);

        double pitch = map.getPitch();
        double[] bounds = map.getBounds();

        double x2 = getWidth();
        double y2 = getHeight();

        double x1p = bounds[0];
        double x2p = bounds[1];
        double y1p = bounds[2];
        double y2p = bounds[3];

        
   
        double ax = x2 *.9 / (x2p - x1p); //pixels/ electrode map unit
        double ay = y2 *.9 / (y2p - y1p);  //pixels/ electrode map unit
        double bx = -x2*x1p/(x2p-x1p); //pixels offset (50%)
        double by = -y2*y1p/(y2p-y1p); //pixels offset (50%)
       
        double rx = map.getPitch() / 2 * ax;  //sizing scale factor
        double ry = map.getPitch() / 2 * ay;  //sizing scale factor


//        g.setColor(Color.white);
//        g.fillRect(0, 0, getWidth(), getHeight());

//      set alpha=0.5 and draw white frame if inset is overlaying (SS)
        if (isOverlaying) {
            g.setComposite(AlphaComposite.getInstance(3, 0.5f));
            g.setColor(Color.white);
            g.setStroke(new BasicStroke(3f));
            g.drawRect(0, 0, getWidth(), getHeight());
        }

        g.setColor(Color.black);
        g.setStroke(new BasicStroke(1f / 2f));
        g.drawRect(0, 0, getWidth(), getHeight());
        g.setFont(new Font("Arial", Font.PLAIN, 11));

        double r = 0, r1 = 0;

        // main drawing loop
        for (int electrode = 1; electrode < nElectrodes; electrode++) {
            if (skipDisconnected && map.isDisconnected(electrode)) {
                continue;
            }
          
            
            double x = (ax * map.getXPosition(electrode) + bx);
            double y = (-ay * map.getYPosition(electrode) + by);
 
            r = (amplitude[electrode] - minAmplitude) / (maxAmplitude - minAmplitude);
            
            if (sigma != null) {
                r1 = (sigma[electrode] * threshold - minAmplitude) /
                     (maxAmplitude - minAmplitude);
                r1 = Math.min(r1, 1);
            }
            
            if (markElectrodes != null &&
                Arrays.binarySearch(markElectrodes, electrode) > 0) {
                drawColored(x, y, rx /* * r*/, ry /* * r*/, Color.pink, g);
            } else {
                if (sigma != null) {
                    Color thresholdColor = (electrode != markerElectrode) ? new Color(.9f, .9f, .9f) : new Color(.9f, .0f, .0f);
                    if (amplitude[electrode] > sigma[electrode] * threshold) {
                        // draw amplitude
                        drawColored(x, y, rx * r, ry * r, Color.black, g);
                        // draw threshold
                        drawColored(x, y, rx * r1, ry * r1, thresholdColor, g);
                    } else {
                        // draw threshold
                        drawColored(x, y, rx * r1, ry * r1, thresholdColor, g);
                        // draw amplitude
                        drawColored(x, y, rx * r, ry * r, Color.black, g);
                    }
                } else {
                    // draw amplitude
                    if (electrode != markerElectrode) {

//                        switch (getElectrodeType(electrode)) {
//                            case AXON:
//                                drawColored(x, y, rx * r, ry * r, Color.green, g);
//                                break;
//                            case BODY:
//                                drawColored(x, y, rx * r, ry * r, Color.red, g);
//                                break;
//                            case DENDRITE:
//                                drawColored(x, y, rx * r, ry * r, Color.black, g);
//                                break;
//                        }

                        // Show negative electrodes in blue during movie playback
                        if (negativeAmplitudes && average[electrode][timepoint] < 0) {
                            drawColored(x, y, rx * r, ry * r, Color.blue, g);
                        } else {
                            drawColored(x, y, rx * r, ry * r, Color.black, g);
                        }
                    } else {
                        drawColored(x, y, rx * r, ry * r, new Color(.9f, .0f, .0f), g);
                    }
                }
            }
        }

        Insets ins = this.getInsets();
        int w = getWidth();
        int h = getHeight();

        // draw the scale bar
        if (scaleBarSize > 0) {
            double dx = scaleBarX *
                        (w - ins.left - ins.right) /*- c.ax * comp.getWidth()*/;
            double dy = scaleBarY *
                        (h - ins.top - ins.bottom) /* - c.ay * comp.getHeight()*/;

            double x = 0 + ins.left + dx;
            double y = h - ins.bottom - dy;

            g.setStroke(new BasicStroke(scaleBarThickness));
            g.setColor(Color.black);
            g.draw(new Line2D.Double(x, y, x + ax * scaleBarSize, y));
        }

        // draw the arrows
//        for (Arrow a : arrows) {
//            double x = ax * a.x + bx;
//            double y = ay * a.y + by;
//
//            // screen coords
//            double xb = x + a.r * cos(a.angle * PI / 180);
//            double yb = y + a.r * sin(a.angle * PI / 180);
//            double xc = x + (a.r + a.length) * cos(a.angle * PI / 180);
//            double yc = y + (a.r + a.length) * sin(a.angle * PI / 180);
//
//            g.setStroke(new BasicStroke(a.lineWidth));
//            g.setColor(Color.gray);
//            g.draw(new Line2D.Double(xb, h - yb, xc, h - yc));
//
//            g.draw(new Line2D.Double(
//                xb, h - yb,
//                xb + a.arrowSize * cos( (a.angle - a.arrowAngle) * PI / 180),
//                h - (yb + a.arrowSize * sin( (a.angle - a.arrowAngle) * PI / 180))));
//            g.draw(new Line2D.Double(
//                xb, h - yb,
//                xb + a.arrowSize * cos( (a.angle + a.arrowAngle) * PI / 180),
//                h - (yb + a.arrowSize * sin( (a.angle + a.arrowAngle) * PI / 180))));
//        }

        PlotPanel.drawBackgroundText(this, g, backgroundText,
                                     backgroundTextX, backgroundTextY,
                                     backgroundTextFont, backgroundTextColor);

        g.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                           RenderingHints.VALUE_ANTIALIAS_OFF);
    }


    public static double compare(float[] template, float[] data) {
        int t0 = MathUtil.minIndex(template);
        int t = MathUtil.minIndex(data);
        double s1 = MathUtil.sumAbs(template);
        double s2 = MathUtil.sumAbs(data);

        // i goes over the data
        double sum = 0;
        int n = 0;
        for (int i = 0; i < data.length; i++) {
            int j = i + (t0 - t);
            if (j >= 0 && j < data.length) {
                sum += Math.pow(template[j] / s1 - data[i] / s2, 2);
                n++;
            }
        }
        return 1000 * sum / n;
    }


    public void setTemplates(int axonTemplateElectrode,
                             int dendriteTemplateElectrode) throws IOException {

        axonTemplate = new float[average[0].length];
        System.arraycopy(average[axonTemplateElectrode], 0, axonTemplate, 0,
                         axonTemplate.length);

        dendriteTemplate = new float[average[0].length];
        System.arraycopy(average[dendriteTemplateElectrode], 0, dendriteTemplate,
                         0, dendriteTemplate.length);

        bodyTemplate = new float[average[0].length];
        System.arraycopy(average[dendriteTemplateElectrode], 0, bodyTemplate,
                         0, bodyTemplate.length);
    }


    public boolean isAxon(int electrode) {
        int t = MathUtil.minIndex(average[electrode]);
        double chi2 = compare(axonTemplate, average[electrode]);
        return (t > t0 && Math.abs(average[electrode][t]) > ampThreshold &&
                chi2 < maxChi2);
    }


    public boolean isDendrite(int electrode) {
        int t = MathUtil.minIndex(average[electrode]);
        double chi2 = compare(dendriteTemplate, average[electrode]);
        return (Math.abs(average[electrode][t]) > ampThreshold && chi2 < maxChi2);
    }


    public static final int AXON = 0;
    public static final int DENDRITE = 1;
    public static final int BODY = 2;
    public static final int NOISE = 3;

    public int getElectrodeType(int electrode) {
        int t = MathUtil.minIndex(average[electrode]);
        double aChi2 = compare(axonTemplate, average[electrode]);
        double dChi2 = compare(dendriteTemplate, average[electrode]);
        double bChi2 = compare(bodyTemplate, average[electrode]);
        
        if (Math.abs(average[electrode][t]) < 5) {
            return NOISE;
        }
        
        if (t > t0 && aChi2 < dChi2 && aChi2 < bChi2 && aChi2 < maxChi2) {
            return AXON;
        } else if (bChi2 <= aChi2 && bChi2 <= dChi2 && bChi2 < maxChi2) {
            return BODY;
        } else {
            return DENDRITE;
        }
    }

}
