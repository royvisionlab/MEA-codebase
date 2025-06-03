package edu.ucsc.neurobiology.vision.plot;

import java.io.*;
import java.util.*;

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import java.awt.image.*;
import javax.swing.*;
import javax.swing.event.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.gui.*;
import edu.ucsc.neurobiology.vision.neuronviewer.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class PlotPanel
extends JPanel implements MouseListener, MouseMotionListener, ChangeListener,
Cloneable, ComponentListener {

    public static String saveFolder = System.getProperty("user.dir");

    public enum Location { UPRIGHT, DOWNRIGHT, UPLEFT, DOWNLEFT; }


    public static final double LEFT = -Integer.MAX_VALUE + 0;
    public static final double CENTER = -Integer.MAX_VALUE + 1;
    public static final double RIGHT = -Integer.MAX_VALUE + 2;
    public static final double TOP = LEFT;
    public static final double BOTTOM = RIGHT;


    //    public static final int DRAW_TO_SCREEN = 0;
    //    public static final int DRAW_TO_FILE = 1;

    final Config configuration = Vision.getInstance().getConfig();


    ArrayList<PlotData> plotData;
    ArrayList<PlotStyle> styles;
    ArrayList<DataPlotter> dataPlotters;
    public AxesBorder axesBorder;
    Axis hAxis, vAxis;
    private double xMin, xMax, yMin, yMax;
    public DrawingControl drawingControl;

    private final int hPadding = 5;
    private boolean legendVisible;
    private boolean drawGrid = false;
    private ArrayList<String> additionalLegend = new ArrayList<String>();
    private Location legendLocation = Location.UPRIGHT;
    private LinkedList<double[]> zoomHistory = new LinkedList<double[]>();
    public ArrayList<ScaledMouseListener> scaledMouseListeners = new ArrayList<ScaledMouseListener>();
    private final JPopupMenu menu = new JPopupMenu();
    private double overscaleFactor = 1;
    protected Selection[] selections;

    private Selection mainSelection, currentSelection;
    protected ArrayList<Selection> passiveSelections;
    private BufferedImage backScreen;
    public Raster raster;
    private boolean doubleBuffered;
    boolean generalRepaintNeeded = false;
    ChangeListener passiveSlectionListener;

    private boolean allowPasiveSelectionAdding = true;
    private boolean allowPassiveSelectionRemoval = true;
    private boolean allowPassiveSelectionManipulation = true;

    JCheckBoxMenuItem showXAxisLabelItem, showXAxisTicksItem, showYAxisLabelItem,
    showYAxisTicksItem;
    JCheckBoxMenuItem showGridItem;
    boolean ignoreNextClick = false;
    ArrayList<String> backgroundText = new ArrayList<String>();
    ArrayList<Double> backgroundTextX = new ArrayList<Double>();
    ArrayList<Double> backgroundTextY = new ArrayList<Double>();
    ArrayList<Font> backgroundTextFont = new ArrayList<Font>();
    ArrayList<Color> backgroundTextColor = new ArrayList<Color>();

    ArrayList<String> options = new ArrayList<String>();
    int selectedOption = 0;
    Rectangle optionsBounds, legendBounds;
    boolean legendHighlighted = false;
    ArrayList<Arrow> arrows = new ArrayList<Arrow>();
    boolean antiallias = false;

    //controls limits for plot, when the user resizes it.
    private double[] limits = new double[]{Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, 
            Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY};


    public PlotPanel() {
        this(null, true, false, false, false);
    }


    public PlotPanel(String name) {
        this(name, true, false, false, false);
    }


    public PlotPanel(String name,
            boolean doubleBuffered,
            boolean allowPasiveSelectionAdding,
            boolean allowPassiveSelectionRemoval,
            boolean allowPassiveSelectionManipulation) {

        super(new PlaceLayout(), false);
        this.setName(name);

        this.doubleBuffered = doubleBuffered;
        this.allowPasiveSelectionAdding = allowPasiveSelectionAdding;
        this.allowPassiveSelectionRemoval = allowPassiveSelectionRemoval;
        this.allowPassiveSelectionManipulation = allowPassiveSelectionManipulation;

        axesBorder = new AxesBorder();
        hAxis = axesBorder.getHorizontalAxis();
        vAxis = axesBorder.getVerticalAxis();
        hAxis.setScreenRange(0, 100);
        vAxis.setScreenRange(0, 100);
        axesBorder.getHorizontalAxis().setChangeListener(this);
        axesBorder.getVerticalAxis().setChangeListener(this);
        setBorder(axesBorder);

        passiveSelections = new ArrayList<Selection>();
        this.addMouseListener(this);
        this.addMouseMotionListener(this);
        this.addComponentListener(this);

        setBackground(Color.white);
        plotData = new ArrayList<PlotData>();
        styles = new ArrayList<PlotStyle>();
        dataPlotters = new ArrayList<DataPlotter>();

        drawingControl = new DrawingControlImpl();
        addDataPlotter(new HistogramPlotter());
        addDataPlotter(new Histogram2DPlotter());
        addDataPlotter(new ScatterPlotPlotter());
        addDataPlotter(new FunctionPlotter());
        addDataPlotter(new ColorPlotPlotter());

        this.legendVisible = true;

        selections = new Selection[] {
                new RectangularSelection(this, hAxis, vAxis),
                new ParallelogramSelection(this, hAxis, vAxis),
                //                     new RotatedRectangularSelection(this, hAxis, vAxis),
                new EllipticSelection(this, hAxis, vAxis),
        };
        mainSelection = selections[0];

        // add the "Zoom In" option
        addSelectionAction(new SelectionAction("Zoom In") {
            public void selectionPerformed(JComponent source, Selection selection) {
                addZoomHistory(getRange());
                Rectangle2D box = selection.getSelection().getBounds2D();
                setRange(box.getX(), box.getX() + box.getWidth(),
                        box.getY(), box.getY() + box.getHeight());

            }
        });

        menu.add(new AbstractAction("Zoom Out") {
            public void actionPerformed(ActionEvent e) {
                zoomOut();
            }
        });
        menu.add(new AbstractAction("Autoscale") {
            public void actionPerformed(ActionEvent e) {
                autoscale();
            }
        });

        menu.add(new JSeparator());
        showGridItem = new JCheckBoxMenuItem("Show Grid Lines", drawGrid);
        showGridItem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                setGridVisible(showGridItem.isSelected());
            }
        });
        //        menu.add(showXAxisItem);
        //        menu.add(showYAxisItem);
        menu.add(showGridItem);

        menu.add(new JSeparator());
        menu.add(new AbstractAction("Rectangle Selection") {
            public void actionPerformed(ActionEvent e) {
                setCurrentSelection(0);
            }
        });
        menu.add(new AbstractAction("Rotated Rectangle Selection") {
            public void actionPerformed(ActionEvent e) {
                setCurrentSelection(1);
            }
        });
        menu.add(new AbstractAction("Ellipse Selection") {
            public void actionPerformed(ActionEvent e) {
                setCurrentSelection(2);
            }
        });
        menu.add(new JSeparator());

        menu.add(new AbstractAction("Save Image...") {
            public void actionPerformed(ActionEvent event) {
                GraphicsIO.saveImage(PlotPanel.this, PlotPanel.this, "Save Image");
            }
        });

        menu.add(new AbstractAction("Copy to Clipboard") {
            public void actionPerformed(ActionEvent event) {
                Clipboard.saveComponentToClipboard(PlotPanel.this);
            }
        });



        menu.add(new JSeparator());

        menu.add(new AbstractAction("Properties") {
            public void actionPerformed(ActionEvent e) {
                new PlotCustomizerDialog(drawingControl, PlotPanel.this);
            }
        });

        if (allowPasiveSelectionAdding) {
            addSelectionAction(new SelectionAction("Add Box") {
                public void selectionPerformed(JComponent source, Selection selection) {
                    addPassiveSelection(selection.getCopy());
                    replotAllData();
                }
            });
        }

        this.setRange(1, 2, 1, 2);
    }
    
    
    public void addZoomHistory(double[] range) {
        zoomHistory.addLast(range);
    }
    
    public void addZoomHistory(Rectangle2D bounds) {
        addZoomHistory(new double[]{bounds.getMinX(), bounds.getMaxX(), bounds.getMinY(), bounds.getMaxY()});
    }

    public void zoomOut() {
        if (zoomHistory.isEmpty()) return;
        setRange(zoomHistory.removeLast());
    }

    
    synchronized public void setDoubleBuffered(boolean doubleBuffered) {
        this.doubleBuffered = doubleBuffered;
    }


    public void setAntiallias(boolean antiallias) {
        this.antiallias = antiallias;
    }


    public void addArrow(Arrow a) {
        arrows.add(a);
    }


    public Selection getPassiveSelection(int i) {
        return (Selection) passiveSelections.get(i);
    }


    public void resetPassiveSelections() {
        passiveSelections.clear();
    }


    public void addPopupAction(Action action) {
        menu.add(action);
    }


    public void mouseDragged(MouseEvent e) {
        if (currentSelection != null) {
            currentSelection.mouseDragged(e);
        }
    }


    public void mouseMoved(MouseEvent e) {
        if (currentSelection != null) {
            currentSelection.mouseMoved(e);
        }
    }


    public void mouseClicked(MouseEvent e) {
        if (ignoreNextClick) {
            ignoreNextClick = false;
            return;
        }

        if (SwingUtilities.isLeftMouseButton(e)) {
            // check whether the options were clicked
            if (optionsBounds.contains(e.getX(), e.getY())) {
                int i = (int) Math.round( (e.getY() - optionsBounds.y) /
                        (optionsBounds.height / (double) options.size()));
                setSelectedOption(i);
                return;
            }

            // deliver click events
            for (int i = 0; i < scaledMouseListeners.size(); i++) {
                ScaledMouseListener l = (ScaledMouseListener) scaledMouseListeners.get(i);
                double x = axesBorder.getHorizontalAxis().getPlotCoord(e.getX());
                double y = axesBorder.getVerticalAxis().getPlotCoord(e.getY());
                l.clickPerformed(PlotPanel.this, e, x, y);
            }
        }

        if (allowPassiveSelectionRemoval &&
                (passiveSelections.size() != 0) &&   
                (SwingUtilities.isLeftMouseButton(e) &&
                        (e.getModifiersEx() & MouseEvent.ALT_DOWN_MASK) == MouseEvent.ALT_DOWN_MASK) 
        ) {
            
            // remove the closest box
            double minD = Double.POSITIVE_INFINITY;
            int jMin = -1;
            double x = hAxis.getPlotCoord(e.getX());
            double y = vAxis.getPlotCoord(e.getY());
            for (int j = 0; j < passiveSelections.size(); j++) {
                Selection s = (Selection) passiveSelections.get(j);
                double d = s.getSelection().getMinimumDistance(x, y);
                if (d < minD) {
                    minD = d;
                    jMin = j;
                }
            }
            passiveSelections.remove(jMin);
            mainSelection.resetSelection();

            if (passiveSlectionListener != null) {
                passiveSlectionListener.stateChanged(new ChangeEvent(this));
            }
        }
    }


    public void mouseEntered(MouseEvent e) {
    }


    public void mouseExited(MouseEvent e) {
    }


    public void mousePressed(MouseEvent e) {

        // deliver pressed events
        for (int i = 0; i < scaledMouseListeners.size(); i++) {

            ScaledMouseListener l = (ScaledMouseListener) scaledMouseListeners.get(i);
            double x = axesBorder.getHorizontalAxis().getPlotCoord(e.getX());
            double y = axesBorder.getVerticalAxis().getPlotCoord(e.getY());
            l.pressPerformed(PlotPanel.this, e, x, y);
        }

        // one and two button mouse friendly.
        if (SwingUtilities.isRightMouseButton(e) || 
                (e.getModifiersEx() & MouseEvent.CTRL_DOWN_MASK) == MouseEvent.CTRL_DOWN_MASK) {
            menu.show(PlotPanel.this, e.getX(), e.getY());
            return;
        }

        // check whether the legend was clicked
        if (SwingUtilities.isLeftMouseButton(e)) {
            if (legendBounds != null && legendBounds.contains(e.getX(), e.getY())) {
                legendHighlighted = true;
                resetBackScreen();
                repaint();
                ignoreNextClick = true;
                return;
            }
        }

        if (!SwingUtilities.isLeftMouseButton(e)) {
            mainSelection.resetSelection();
            return;
        }

        //if beginning drag without modifiers, make selection area
        if (!e.isControlDown() && !e.isShiftDown() && !e.isAltDown()) {  //pass on to listener (cluster mover)
            for (int i = 0; i < passiveSelections.size(); i++) {
                Selection s = (Selection) passiveSelections.get(i);
                if (s.nearWhichControlPoint(e) != -1) {
                    currentSelection = s;
                    currentSelection.mousePressed(e);
                    return;
                }
            }

            currentSelection = mainSelection;
            currentSelection.mousePressed(e);
        }
    }


    public void mouseReleased(MouseEvent e) {

        // deliver pressed events
        for (int i = 0; i < scaledMouseListeners.size(); i++) {
            ScaledMouseListener l = (ScaledMouseListener) scaledMouseListeners.get(i);
            double x = axesBorder.getHorizontalAxis().getPlotCoord(e.getX());
            double y = axesBorder.getVerticalAxis().getPlotCoord(e.getY());
            l.releasePerformed(PlotPanel.this, e, x, y);
        }

        if (legendBounds != null && legendBounds.contains(e.getX(), e.getY())) {
            legendHighlighted = false;
            resetBackScreen();
            repaint();
            return;
        }

        if (!SwingUtilities.isLeftMouseButton(e)) {
            return;
        }

        if (currentSelection != null) {
            currentSelection.mouseReleased(e);
            if (!currentSelection.active && passiveSlectionListener != null) {
                passiveSlectionListener.stateChanged(new ChangeEvent(this));
            }
            currentSelection = null;
        }
    }


    public int getPassiveSelectionsCount() {
        return passiveSelections.size();
    }


    public void addPassiveSelectionChangeListener(ChangeListener list) {
        this.passiveSlectionListener = list;
    }


    protected void setCurrentSelection(Selection selection) {
        mainSelection = selection;
    }


    public void addPassiveSelection(Selection s) {
        s.active = false;
        s.visible = true;
        passiveSelections.add(s);

        if (passiveSlectionListener != null) {
            passiveSlectionListener.stateChanged(new ChangeEvent(this));
        }

        this.repaint();
    }


    private Selection getCurrentSelection() {
        return mainSelection;
    }


    synchronized public void resetBackScreen() {
        generalRepaintNeeded = true;
    }


    final synchronized public void paintComponent(Graphics g) {
        //        System.err.println("paintComponent in " + this.getName() + " " + N);
        if (doubleBuffered) {
            // recreate the backscreen if needed
            if (backScreen == null || getWidth() != backScreen.getWidth() ||
                    getHeight() != backScreen.getHeight()) {
                
                backScreen = new BufferedImage(
                        getWidth(), getHeight(), BufferedImage.TYPE_INT_ARGB);
                raster = backScreen.getWritableTile(0, 0);
                generalRepaintNeeded = true;
            }

            // repaint the whole data if needed
            if (generalRepaintNeeded) {
                Graphics backGraphics = backScreen.createGraphics();
                super.paintComponent(backGraphics);
                spiPaintComponent(backGraphics, raster);
            }

            // flush the backscreen to the screen
            ( (Graphics2D) g).drawImage(backScreen, 0, 0, this);

            // draw the selections
            if (mainSelection != null && mainSelection.isValidSelection()) {
                mainSelection.paintSelection(g);
            }
            for (int i = 0; i < passiveSelections.size(); i++) {
                Selection s = (Selection) passiveSelections.get(i);
                if (s != null) {
                    s.paintSelection(g);
                }
            }
        } else {
            // repaint the whole data directly to the screen
            super.paintComponent(g);
            spiPaintComponent(g, null);

            // draw the selections
            if (mainSelection != null && mainSelection.isValidSelection()) {
                mainSelection.paintSelection(g);
            }
            for (int i = 0; i < passiveSelections.size(); i++) {
                Selection s = (Selection) passiveSelections.get(i);
                if (s != null) {
                    s.paintSelection(g);
                }
            }
        }

        generalRepaintNeeded = false;
        //        System.err.println("paintComponent out " + this.getName() + " " + N);
    }

    public ArrayList<PlotData> getPlotData() {
        return plotData;
    }


    /**
     * As the panel is being set up, it is possible for this to be called before 
     * mainSelection and passiveSelections are initialized.  This seems to happen 
     * on Linux (Ubuntu Gnome) but not Windows or Mac!  Therefore, if mainSelection
     * is not initialized yet, we abort and wait until later.
     */
    public void resizeHappened() {
        if (mainSelection == null) return;
        
        mainSelection.componentResized();
        for (int i = 0; i < passiveSelections.size(); i++) {
            Selection s = (Selection) passiveSelections.get(i);
            s.componentResized();
            s.parentComponent.repaint();
        }
    }


    public void stateChanged(ChangeEvent e) {
        resizeHappened();
    }


    public Axis getXAxis() {
        return axesBorder.getHorizontalAxis();
    }


    public Axis getYAxis() {
        return axesBorder.getVerticalAxis();
    }


    public void setCurrentSelection(int i) {
        this.setCurrentSelection(selections[i]);
    }


    public void addClickListener(ScaledMouseListener l) {
        scaledMouseListeners.add(l);
    }


    public void saveAsEPS(File file, double width, double height) throws IOException {
        GraphicsIO.saveAsEPS(this, file.getAbsolutePath(), width, height, false);
    }


    public void saveAsEPS(File file, double width, double height, boolean insetsSize) throws
    IOException {
        GraphicsIO.saveAsEPS(this, file.getAbsolutePath(), width, height, insetsSize);
    }


    //    public void saveAsEPS(File file, double width, double height) throws IOException {
    //        saveAsEPS(file, width, height, true);
    //    }

    /*
     public void saveAsEPS(File file, double width, double height, boolean boxSize) throws
            IOException {

            Rectangle oldBounds = getBounds();
            this.replotAllData();

            // 72 is default dpi.
            if (boxSize) {
                Insets i = this.getInsets();
                setSize( (int) Math.round(width * 72.0 + i.left + i.right),
                        (int) Math.round(height * 72.0 + i.top + i.bottom));
            } else {
     setSize( (int) Math.round(width * 72.0), (int) Math.round(height * 72.0));
            }

            EpsGraphics2D g = new EpsGraphics2D(
                "image", getXAxis().getLabelFont(), 0, 0, getWidth(), getHeight());
            g.setAccurateTextMode(true);
            g.setColorDepth(g.CYMK);

            paintBorder(g);
            spiPaintComponent(g, null);

            FileWriter writer = new FileWriter(file);
            writer.write(g.toString());
            writer.flush();
            writer.close();

            setBounds(oldBounds);
        }
     */

    public void saveAsPNG(File file) throws IOException {
        GraphicsIO.saveAsPNG(this, file, -1, -1);
        //        BufferedImage bi = new BufferedImage(
        //            getWidth(), getHeight(), BufferedImage.TYPE_INT_RGB);
        //        Graphics2D g2 = bi.createGraphics();
        //        g2.setColor(Color.white);
        //        g2.fillRect(0, 0, getWidth(), getHeight());
        //        paintBorder(g2);
        //        spiPaintComponent(g2, null);
        ////            paint(g2);  paint leads to errors when writing report in serial
        ////            neuron finder, due to multithreading problems
        //        ImageIO.write(bi, "PNG", file);
        //
        //        //This does not always work correctly
        //        //  new PNGExportFileType().exportToFile(file, this, null, null, "");
        ////        new PNGExportFileType().exportToFile(file, this, null, null, "");
    }


    public void saveAsPNG(String fileName) throws IOException {
        saveAsPNG(new File(fileName));
    }


    public void saveAsPNG(File file, int w, int h) throws IOException {
        GraphicsIO.saveAsPNG(this, file, w, h);
    }


    public void saveAsPNG(String fileName, int w, int h) throws IOException {
        saveAsPNG(new File(fileName), w, h);
    }


    public void saveSerial(File file) throws IOException {
        ObjectOutputStream os = new ObjectOutputStream(new FileOutputStream(file));

        os.writeDouble(xMin);
        os.writeDouble(xMax);
        os.writeDouble(yMin);
        os.writeDouble(yMax);

        os.writeUTF(axesBorder.getHorizontalAxis().getLabel());
        os.writeUTF(axesBorder.getVerticalAxis().getLabel());

        os.writeObject(plotData);
        os.writeObject(styles);

        os.writeObject(additionalLegend);
        os.writeBoolean(legendVisible);
        os.writeInt(legendLocation.ordinal());

        os.close();
    }


    public void setOverscaleFactor(double overscaleFactor) {
        this.overscaleFactor = overscaleFactor;
    }


    public int getPlotsCount() {
        return plotData.size();
    }


    /*
        public double[] getSelection() {
            Rectangle box = currentSelection.getBoundingBox();
            Axis h = axesBorder.getHorizontalAxis();
            Axis v = axesBorder.getVerticalAxis();
            Insets i = PlotPanel.this.getInsets();
            return new double[] {
                h.getPlotCoord(box.x),
                v.getPlotCoord(box.y),
                h.getPlotCoord(box.x + box.width),
                v.getPlotCoord(box.y + box.height)
            };
        }
     */

    public PlotPanel autoscale() {
        //FIXME
        //        this.setXRange(0, 1);
        //        this.setYRange(0, 1);

        double[] union = null;

        for (int i = 0; i < plotData.size(); i++) {
            PlotData data = (PlotData) plotData.get(i);

            double[] bounds;
            if (data instanceof ScatterPlotData) {
                bounds = PlotUtil.getPlotBounds( (ScatterPlotData) data);
            } else if (data instanceof HistogramData) {
                bounds = PlotUtil.getPlotBounds( (HistogramData) data);
            } else if (data instanceof Histogram2DData) {
                Histogram2DData h2d = (Histogram2DData) data;
                bounds = new double[4];
                bounds[0] = h2d.getMinX();
                bounds[1] = h2d.getMaxX();
                bounds[2] = h2d.getMinY();
                bounds[3] = h2d.getMaxY();
            } else {
                continue;
            }

            if (union == null) {
                union = new double[4];
                System.arraycopy(bounds, 0, union, 0, 4);
            } else if (bounds != null) {
                if (bounds[0] < union[0]) {
                    union[0] = bounds[0];
                }
                if (bounds[1] > union[1]) {
                    union[1] = bounds[1];
                }
                if (bounds[2] < union[2]) {
                    union[2] = bounds[2];
                }
                if (bounds[3] > union[3]) {
                    union[3] = bounds[3];
                }
            }
        }

        // All this messy ifs are needed to solve the problem of the PlotPanel
        // freezing on the screen in some JDKs. They just do not allow the ranges
        // to be equal. In Sun "1.4.2_03" it does not freeze. It freezes with
        // JBuilder9's 1.4.1_02-b06.
        if (union != null) {
            if (union[0] == union[1]) {
                union[0]--;
                union[1]++;
            }
            if (union[2] == union[3]) {
                union[2]--;
                union[3]++;
            }

            // added to provide overscale
            double x0 = (union[0] + union[1]) / 2;
            double y0 = (union[2] + union[3]) / 2;
            union[0] = x0 + (union[0] - x0) * overscaleFactor;
            union[1] = x0 + (union[1] - x0) * overscaleFactor;
            union[2] = y0 + (union[2] - y0) * overscaleFactor;
            union[3] = y0 + (union[3] - y0) * overscaleFactor;

            setRange(union);
        } else {
            setRange(0, 1, 0, 1);
        }

        //        IOUtil.printArray(this.getRange());

        return this;
    }


    public void pad() {
        double[] r = getRange();
        setRange(r[0] + ( -0.1) * (r[1] - r[0]),
                r[1] + ( +0.1) * (r[1] - r[0]),
                r[2] + ( -0.1) * (r[3] - r[2]),
                r[3] + ( +0.1) * (r[3] - r[2]));
    }


    public void padX() {
        double[] r = getRange();
        setRange(r[0] + ( -0.1) * (r[1] - r[0]),
                r[1] + ( +0.1) * (r[1] - r[0]),
                r[2],
                r[3]);
    }


    public void padY() {
        double[] r = getRange();
        setRange(r[0],
                r[1],
                r[2] + ( -0.1) * (r[3] - r[2]),
                r[3] + ( +0.1) * (r[3] - r[2]));
    }


    public void saveAsASCII(int index) {
        PlotData data = (PlotData) plotData.get(index);

        if (data instanceof ScatterPlotData) {
            double[] p = new double[3];
            int nPoints = ( (ScatterPlotData) data).getPointCount();
            System.out.println("==> Start: " + data.getDescription());
            for (int i = 0; i < nPoints; i++) {
                ( (ScatterPlotData) data).getDataPoint(i, p);
                System.out.println(p[0] + "\t" + p[1]);
            }
            System.out.println("==> End: " + data.getDescription());
        } else if (data instanceof HistogramData) {
            HistogramData h = (HistogramData) data;
            double x1 = h.getMin();
            double dx = h.getBinInterval();
            for (int i = 0; i < h.getBinCount(); i++) {
                System.out.println( (x1 + i * dx) + "\t" + h.getBin(i));
            }
        } else {
            System.out.println("Unable to save, sorry!");
            System.out.println(data);
        }
    }


    public void addSelectionAction(final SelectionAction selectionAction) {
        for (int sIndex = 0; sIndex < selections.length; sIndex++) {
            final Selection s = selections[sIndex];
            String name = selectionAction.getDescription();
            s.addPossibleAction(new AbstractAction(name) {
                public void actionPerformed(ActionEvent e) {
                    selectionAction.selectionPerformed(PlotPanel.this, s);
                    s.resetSelection();
                }
            });
        }
    }


    public void setXAxisVisible(boolean visible) {
        axesBorder.setXAxisVisible(visible);
        //        showXAxisItem.setSelected(visible);

        replotAllData();
    }


    public void setYAxisVisible(boolean visible) {
        axesBorder.setYAxisVisible(visible);
        //        showYAxisItem.setSelected(visible);

        replotAllData();
    }


    public void setGridVisible(boolean visible) {
        this.drawGrid = visible;
        showGridItem.setSelected(visible);

        replotAllData();
    }


    public void setAxisVisible(boolean visible) {
        axesBorder.setXAxisVisible(visible);
        axesBorder.setYAxisVisible(visible);
        //        showXAxisItem.setSelected(visible);
        //        showYAxisItem.setSelected(visible);

        replotAllData();
    }


    /*
        public boolean isXAxisVisible() {
            return axesBorder.isXAxisVisible();
        }


        public boolean isYAxisVisible() {
            return axesBorder.isYAxisVisible();
        }
     */

    private boolean insidePlottingArea(MouseEvent e) {
        int x = e.getX();
        int y = e.getY();
        Insets i = axesBorder.getBorderInsets(this);

        if ( (x > i.left) && (x < getWidth() - i.right) &&
                (y > i.top) && (y < getHeight() - i.bottom)) {
            return true;
        } else {
            return false;
        }
    }


    private boolean insideHorizontalAxis(MouseEvent e) {
        int x = e.getX();
        int y = e.getY();
        Insets i = axesBorder.getBorderInsets(this);

        if ( (x > i.left) && (y > getHeight() - i.bottom)) {
            return true;
        } else {
            return false;
        }
    }


    private boolean insideVerticalAxis(MouseEvent e) {
        int x = e.getX();
        int y = e.getY();
        Insets i = axesBorder.getBorderInsets(this);

        if ( (x < i.left) && (y < getHeight() - i.bottom)) {
            return true;
        } else {
            return false;
        }
    }


    public void setXRange(final double xMin, final double xMax) {
        if ( (this.xMin == xMin) && (this.xMax == xMax)) {
            return;
        } else {
            this.xMin = xMin;
            this.xMax = xMax;
            for (int i = 0; i < plotData.size(); i++) {
                PlotData data = (PlotData) plotData.get(i);
                if (data instanceof ScaleChangeListener) {
                    ( (ScaleChangeListener) data).xScaleChanged(xMin, xMax);
                }
            }
            axesBorder.getHorizontalAxis().setPlotRange(xMin, xMax);
        }

        replotAllData();
    }


    public void setYRange(final double yMin, final double yMax) {
        if ( (this.yMin == yMin) && (this.yMax == yMax)) {
            return;
        } else {
            this.yMin = yMin;
            this.yMax = yMax;
            for (int i = 0; i < plotData.size(); i++) {
                PlotData data = (PlotData) plotData.get(i);
                if (data instanceof ScaleChangeListener) {
                    ( (ScaleChangeListener) data).yScaleChanged(yMin, yMax);
                }
            }
            axesBorder.getVerticalAxis().setPlotRange(yMin, yMax);
        }

        replotAllData();
    }


    public void setLimits(double[] limits) {
        System.arraycopy(limits, 0, this.limits, 0, limits.length);
    }
    
    
    public void setRange(double xMin, double xMax, double yMin, double yMax) {
        setXRange(xMin, xMax);
        setYRange(yMin, yMax);
    }

    public void setRange(double[] bounds) {
        setRange(bounds[0], bounds[1], bounds[2], bounds[3]);
    }
    
    public void setRange(Rectangle2D bounds) {
        setRange(bounds.getMinX(), bounds.getMaxX(), bounds.getMinY(), bounds.getMaxY());
    }


    public double[] getRange() {
        return new double[] {
                xMin, xMax, yMin, yMax};
    }

    public double[] getLimits() {
        return limits;
    }


    public void setLabels(String xLabel, String yLabel) {
        axesBorder.getHorizontalAxis().setLabel(xLabel);
        axesBorder.getVerticalAxis().setLabel(yLabel);
    }


    public String[] getLabels() {
        return new String[] {
                axesBorder.getHorizontalAxis().getLabel(),
                axesBorder.getVerticalAxis().getLabel()
        };
    }


    private class DrawingControlImpl implements DrawingControl {
        public void updateNeeded(Object data) {
            replotAllData();
        }
    }


    synchronized public void addDataPlotter(DataPlotter plotter) {
        dataPlotters.add(plotter);
    }


    int N = 0;
    public void addData(PlotData data, PlotStyle style) {
        boolean accepted = false;
        N++;
        
        for (int i = 0; i < dataPlotters.size(); i++) {
            DataPlotter plotter = (DataPlotter) dataPlotters.get(i);
            if (plotter.accept(data, style)) {
                accepted = true;
                break;
            }
        }
        
        if (accepted) {
            plotData.add(data);
            if (data instanceof ChangeableData) {
                ( (ChangeableData) data).setDrawingControl(drawingControl);
            }
            styles.add(style);

            replotAllData();
        } else {
            throw new IllegalArgumentException(
                    "Illegal Data or Style. Data: " + data + ", Style: " + style);
        }
        
        //        System.err.println("addData " + this.getName() + " " + N);
    }


    public PlotStyle getStyleWithName(String name) {
        for (PlotStyle s : styles) {
            if (s.getDescription().equals(name)) {
                return s;
            }
        }
        return null;
    }


    /*
        public void removeStyleWithName(String name) {
            for (PlotStyle s : styles) {
                if (s.getDescription().equals(name)) {
                    return s;
                }
            }
            return null;
        }
     */

    public void addData(PlotData data, String styleString) {
        PlotStyle style;
        if (data instanceof ScatterPlot) {
            style = new ScatterPlotStyle(styleString);
        } else {
            throw new Error("This method cannot be called for " + data.getClass().getName());
        }
        addData(data, style);
    }


    synchronized public void removeData(PlotData data) {
        if (data == null) {
            throw new IllegalArgumentException("Null Data!");
        }

        int index = plotData.indexOf(data);
        if (index < 0) {
            throw new IllegalArgumentException("No Such Data!");
        } else {
            plotData.remove(index);
            styles.remove(index);

            replotAllData();
        }
    }


    synchronized public void removeDataOfType(Class c) {
        for (int i = 0; i < plotData.size(); ) {
            PlotData data = (PlotData) plotData.get(i);
            if (data.getClass().equals(c)) {
                plotData.remove(i);
                styles.remove(i);
            } else {
                i++;
            }
        }

        replotAllData();
    }


    synchronized public void removeData(int index) {
        if (index < 0) {
            throw new IllegalArgumentException("No Such Data!");
        } else {
            plotData.remove(index);
            styles.remove(index);

            replotAllData();
        }
    }


    synchronized public void removeDataWithStyle(PlotStyle style) {
        for (int i = 0; i < plotData.size(); ) {
            PlotStyle s = (PlotStyle) styles.get(i);
            if (s == style) {
                plotData.remove(i);
                styles.remove(i);
            } else {
                i++;
            }
        }

        replotAllData();
    }


    synchronized public void addData(DataIterator dataIterator, PlotStyle style) {
        plotData.add(dataIterator);
        styles.add(style);

        replotAllData();
    }


    synchronized public void removeAllData() {
        plotData.clear();
        styles.clear();

        replotAllData();
    }


    public void moveToBack(int dataIndex) {
        PlotData data = plotData.remove(dataIndex);
        PlotStyle style = styles.remove(dataIndex);
        plotData.add(0, data);
        styles.add(0, style);
    }


    private DataPlotter getPlotterFor(Object data, Object style) {
        for (int i = 0; i < dataPlotters.size(); i++) {
            DataPlotter plotter = (DataPlotter) dataPlotters.get(i);
            if (plotter.accept(data, style)) {
                return plotter;
            }
        }
        return null;
    }


    boolean left = true, right = true, top = true, bottom = true;

    public void setDrawBoxLines(boolean left, boolean right, boolean top, boolean bottom) {
        this.left = left;
        this.right = right;
        this.top = top;
        this.bottom = bottom;
    }


    public void spiPaintComponent(Graphics _g, Raster raster) {
        Graphics2D g = (Graphics2D) _g;

        if (antiallias) {
            g.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                    RenderingHints.VALUE_ANTIALIAS_ON);
        }

        // check the bounds
        if ( (xMin == 0 && xMax == 0) || (yMin == 0 && yMax == 0)) {
            xMin = 0;
            xMax = 1;
            yMin = 0;
            yMax = 1;
        }

        // draw
        Insets ins = this.getInsets();
        int w = getWidth();
        int h = getHeight();

        g.setColor(Color.black);
        float strokeW = 1.f / 2.f;
        ( (Graphics2D) g).setStroke(new BasicStroke(strokeW));

        //        g.drawRect(ins.left, ins.top, w - ins.left - ins.right, h - ins.top - ins.bottom);

        int x1 = ins.left;
        int y1 = ins.top;
        int x2 = w - ins.right;
        int y2 = h - ins.bottom;

        if (left) {
            g.drawLine(x1, y1, x1, y2);
        }
        if (right) {
            g.drawLine(x2, y1, x2, y2);
        }
        if (top) {
            g.drawLine(x1, y1, x2, y1);
        }
        if (bottom) {
            g.drawLine(x1, y2, x2, y2);
        }

        g.setClip(new Rectangle2D.Double(
                ins.left + 0.5 * strokeW, ins.top + 0.5 * strokeW,
                w - ins.left - ins.right - 1 * strokeW,
                h - ins.top - ins.bottom - 1 * strokeW));

        // draw grids
        if (drawGrid) {
            g.setColor(Color.darkGray);
            int x0 = (int) Math.round(hAxis.getScreenCoord(0));
            int y0 = (int) Math.round(vAxis.getScreenCoord(0));
            g.drawLine( (int) hAxis.getScreenCoord(xMin), y0,
                    (int) hAxis.getScreenCoord(xMax), y0);
            g.drawLine(x0, (int) vAxis.getScreenCoord(yMin),
                    x0, (int) vAxis.getScreenCoord(yMax));
        }

        // draw the plots
        Axis horAxis = axesBorder.getHorizontalAxis();
        Axis verAxis = axesBorder.getVerticalAxis();
        for (int i = 0; i < plotData.size(); i++) {
            PlotData data = (PlotData) plotData.get(i);
            Object style = styles.get(i);
            
            if (data instanceof DataIterator) {
                DataIterator iterator = (DataIterator) data;
                iterator.reset();
                while (iterator.hasNextData()) {
                    Object nextData = iterator.nextData();
                    getPlotterFor(nextData, style).draw(horAxis, verAxis, g, raster, nextData, style);
                }
            } else {
                getPlotterFor(data, style).draw(horAxis, verAxis, g, raster, data, style);
            }
        }

        // draw the scale bar
        if (scaleBarSize > 0) {
            double dx = scaleBarX * (w - ins.left - ins.right);
            double dy = scaleBarY * (h - ins.top - ins.bottom);

            double x = 0 + ins.left + dx;
            double y = h - ins.bottom - dy;
            double scaleBarW = hAxis.getScreenCoord(scaleBarSize) -
            hAxis.getScreenCoord(0);

            g.setStroke(new BasicStroke(scaleBarThickness));
            g.setColor(Color.black);
            g.draw(new Line2D.Double(x, y, x + scaleBarW, y));
        }

        // draw the arrows
        for (Arrow a : arrows) {
            g.setColor(Color.black);
            g.setStroke(new BasicStroke(a.lineWidth));
            double x1a = hAxis.getScreenCoord(a.x1);
            double x2a = hAxis.getScreenCoord(a.x2);
            double y1a = vAxis.getScreenCoord(a.y1);
            double y2a = vAxis.getScreenCoord(a.y2);

            //draw body
            g.draw(new Line2D.Double(x1a, y1a,
                    x2a, y2a));

            //draw head
            double[] seg1 = {-a.headLength, +a.headWidth/2};
            double[] seg2 = {-a.headLength, -a.headWidth/2};

            rotate2D(seg1, getAngle(x1a, y1a, x2a, y2a));
            rotate2D(seg2, getAngle(x1a, y1a, x2a, y2a));
            //       	double[] xPoints = {x2a, x2a-a.headLength, x2a-a.headLength, x2a};
            //      	double[] yPoints = {y2a, y2a+a.headWidth/2, y2a-a.headWidth/2, y2a};
            Polygon arrowHead = new Polygon();
            arrowHead.addPoint((int) x2a, (int) y2a);
            arrowHead.addPoint((int) (x2a+seg1[0]), (int) (y2a+seg1[1]));
            arrowHead.addPoint((int) (x2a+seg2[0]), (int) (y2a+seg2[1]));
            arrowHead.addPoint((int) x2a, (int) y2a);

            g.draw(arrowHead);
            g.fill(arrowHead);

            //    	hAxis.getScreenCoord(a.x);
            //            double x = ax * a.x + bx;
            //            double y = ay * a.y + by;

            // screen coords
            //            double xb = a.x + a.r * cos(a.angle * PI / 180);
            //            double yb = a.y + a.r * sin(a.angle * PI / 180);
            //            double xc = a.x + (a.r + a.length) * cos(a.angle * PI / 180);
            //            double yc = a.y + (a.r + a.length) * sin(a.angle * PI / 180);
            //
            //            g.setStroke(new BasicStroke(a.lineWidth));
            //            g.setColor(Color.black);
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
        }

        updateLegend(g);

        g.setClip(0, 0, w, h);
        drawBackgroundText(this, g, backgroundText, backgroundTextX, backgroundTextY,
                backgroundTextFont, backgroundTextColor);

        g.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                RenderingHints.VALUE_ANTIALIAS_OFF);
    }


    public void printComponent(Graphics g) {
        this.paintComponent(g);
        this.paintBorder(g);
    }


    public void setLegendVisible(boolean legendVisible) {
        this.legendVisible = legendVisible;

        replotAllData();
    }


    /*
     * Utility for arrow drawing 
     * 
     */

    public void rotate2D(double[] point, double angle) {
        double newX = point[0]*Math.cos(angle)-point[1]*Math.sin(angle);
        double newY = point[0]*Math.sin(angle)+point[1]*Math.cos(angle);
        point[0] = newX;
        point[1] = newY;
    }


    /*
     * Utility for arrow drawing 
     * 
     */  
    public double getAngle(double x1, double y1, double x2, double y2) {
        double angle = Math.atan((y2-y1)/(x2-x1));
        if((x2-x1)<0) {  		
            angle += Math.PI ;  			
        }
        return angle;
    }


    public boolean isLegendVisible() {
        return legendVisible;
    }


    public void setLegendLocation(Location legendLocation) {
        this.legendLocation = legendLocation;

        replotAllData();
    }


    public Location getLegendLocation() {
        return legendLocation;
    }


    public void setAdditionalLegend(ArrayList<String> additionalLegend) {
        this.additionalLegend = additionalLegend;
        replotAllData();
    }


    public void addToLegend(String s) {
        this.additionalLegend.add(s);
        replotAllData();
    }


    public void addBackgroundText(String text) {
        addBackgroundText(text, 0, 0, null, Color.black);
    }


    public void removeAllBackgroundText() {
        backgroundText.clear();
        backgroundTextX.clear();
        backgroundTextY.clear();
        backgroundTextFont.clear();
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


    public void clearOptions() {
        this.options.clear();
    }


    public void clearLegend() {
        this.additionalLegend.clear();
    }


    public void addOption(String option) {
        this.options.add(option);
    }


    public int getSelectedOption() {
        return selectedOption;
    }


    public void setSelectedOption(int selectedOption) {
        this.selectedOption = selectedOption;
        resetBackScreen();
        repaint();
    }


    int vPadding = 5;

    private void updateLegend(Graphics _g) {
        Graphics2D g = (Graphics2D) _g;
        //        g.setFont(new Font("sans", Font.PLAIN, 12));
        int fontH = g.getFontMetrics().getAscent() + 2;

        // draw the legend
        if (legendVisible) {
            int n = additionalLegend.size();
            int legendW = 0;
            for (int i = 0; i < n; i++) {
                String line = (String) additionalLegend.get(i);
                int lineW = g.getFontMetrics().stringWidth(line);
                if (lineW > legendW) {
                    legendW = lineW;
                }
            }
            legendW += 2 * hPadding;

            Insets ins = this.getInsets();

            switch (legendLocation) {
            case UPLEFT:
                legendBounds = new Rectangle(
                        ins.left + vPadding,
                        ins.top + vPadding,
                        legendW, n * fontH);
                break;

            case UPRIGHT:
                legendBounds = new Rectangle(
                        getWidth() - ins.right - legendW - vPadding,
                        ins.top + vPadding,
                        legendW, n * fontH);
                break;

            case DOWNLEFT:
                legendBounds = new Rectangle(
                        ins.left + vPadding,
                        getHeight() - ins.bottom - n * fontH - vPadding,
                        legendW, n * fontH);
                break;

            case DOWNRIGHT:
                legendBounds = new Rectangle(
                        getWidth() - ins.right - legendW - vPadding,
                        getHeight() - ins.bottom - n * fontH - vPadding,
                        legendW, n * fontH);
                break;
            }

            // draw the box over the legend
            if (legendHighlighted) {
                g.setColor(Color.white);
                g.fill(legendBounds);
            }

            // draw the texts
            g.setColor(Color.black);
            for (int i = 0; i < n; i++) {
                String line = (String) additionalLegend.get(i);
                //                int lineW = g.getFontMetrics().stringWidth(line);
                g.drawString(line, legendBounds.x + hPadding,
                        legendBounds.y + (i + 1) * fontH - 2);
            }
        }

        // draw options  //FIXME replace by simple components in the plotPanel
        int n = options.size();
        int legendW = 0;
        for (int i = 0; i < n; i++) {
            String line = (String) options.get(i);
            int lineW = g.getFontMetrics().stringWidth(line);
            if (lineW > legendW) {
                legendW = lineW;
            }
        }
        legendW += 2 * hPadding;

        Insets ins = this.getInsets();
        optionsBounds = new Rectangle(
                ins.left + 4, ins.top + 4, legendW, n * fontH);

        // draw the box over the legend
        if (options.size() > 0) {
            g.setColor(new Color(0.9f, .9f, .9f, 0.5f));
            g.fill(optionsBounds);
        }

        // draw the texts
        for (int i = 0; i < n; i++) {
            String line = (String) options.get(i);
            //            int lineW = g.getFontMetrics().stringWidth(line);
            if (i == selectedOption) {
                g.setColor(new Color(0.5f, .5f, .5f, 0.5f));
                g.drawRect(optionsBounds.x, optionsBounds.y + i * fontH,
                        optionsBounds.width - 1, fontH);
            }
            g.setColor(Color.black);
            g.drawString(line, optionsBounds.x + hPadding,
                    optionsBounds.y + (i + 1) * fontH - 2);
        }

    }


    public void replotAllData() {
        //        System.err.println("replotAllData in " + this.getName() + " " + N);
        resetBackScreen();
        repaint();
        //        System.err.println("replotAllData out " + this.getName() + " " + N);
    }


    public void componentResized(ComponentEvent e) {
        resizeHappened();
    }


    public void componentMoved(ComponentEvent e) {
    }


    public void componentShown(ComponentEvent e) {
        resizeHappened();
    }


    public void componentHidden(ComponentEvent e) {
    }


    public void setAxesType(AxisType hType, AxisType vType) {
        hAxis.setAxisType(hType);
        vAxis.setAxisType(vType);
        this.replotAllData();
    }


    public void setLabelFont(Font font) {
        hAxis.setLabelFont(font);
        vAxis.setLabelFont(font);
        replotAllData();
    }


    public Font getLabelFont() {
        return hAxis.getLabelFont();
    }


    private static int getInt(HashMap<String, String> m, String name) {
        return Integer.parseInt(m.get(name));
    }


    private static double getDouble(HashMap<String, String> m, String name) {
        return Double.parseDouble(m.get(name));
    }


    private static boolean getBool(HashMap<String, String> m, String name) {
        return Boolean.valueOf(m.get(name)).booleanValue();
    }


    public void loadStyle(String style, Config config) {
        HashMap<String, String> m = config.getParameterList(style + " PlotPanel Style");
        this.setRange(getDouble(m, "x1"), getDouble(m, "x2"),
                getDouble(m, "y1"), getDouble(m, "y2"));
        this.setLabelFont(Font.decode(m.get("font")));
        this.setAxesType(AxisType.valueOf(m.get("xAxis Type")),
                AxisType.valueOf(m.get("yAxis Type")));
        this.setLabels(m.get("xAxis Label"), m.get("yAxis Label"));
        this.getXAxis().setShowlabel(getBool(m, "show xAxis Label"));
        this.getYAxis().setShowlabel(getBool(m, "show yAxis Label"));
        this.getXAxis().setShowTicks(getBool(m, "show xAxis Ticks"));
        this.getYAxis().setShowTicks(getBool(m, "show yAxis Ticks"));

        this.getXAxis().setTickSize(getInt(m, "x Tick Size"));
        this.getYAxis().setTickSize(getInt(m, "y Tick Size"));

        this.getXAxis().setTickToMarkSpacing(getInt(m, "x Tick to Mark Spacing"));
        this.getYAxis().setTickToMarkSpacing(getInt(m, "y Tick to Mark Spacing"));

        this.getXAxis().setMarkToLabelSpacing(getInt(m, "x Mark to Label Spacing"));
        this.getYAxis().setMarkToLabelSpacing(getInt(m, "y Mark to Label Spacing"));

        this.axesBorder.setTopPadding(getInt(m, "top Padding"));
        this.axesBorder.setBottomPadding(getInt(m, "bottom Padding"));
        this.axesBorder.setLeftPadding(getInt(m, "left Padding"));
        this.axesBorder.setRightPadding(getInt(m, "right Padding"));
    }


    double scaleBarSize = -1;
    float scaleBarThickness = 1;
    double scaleBarX = -1;
    double scaleBarY = -1;
    public void setScaleBar(double scaleBarSize, float scaleBarThickness,
            double scaleBarX, double scaleBarY) {
        if (scaleBarSize < 0) {
            throw new IllegalArgumentException("scaleBar cannot be negative");
        }
        this.scaleBarSize = scaleBarSize;
        this.scaleBarThickness = scaleBarThickness;
        this.scaleBarX = scaleBarX;
        this.scaleBarY = scaleBarY;
    }


    public static void drawBackgroundText(JComponent c, Graphics2D g,
            ArrayList<String> backgroundText,
            ArrayList<Double> backgroundTextX,
            ArrayList<Double> backgroundTextY,
            ArrayList<Font> backgroundTextFont,
            ArrayList<Color> backgroundTextColor) {
        
        int dx = 4, dy = 2;
        Insets ins = c.getInsets();
        JLabel text = new JLabel();
        
        for (int i = 0; i < backgroundText.size(); i++) {
            text.setText(backgroundText.get(i));
            text.setFont(backgroundTextFont.get(i));
            text.setForeground(backgroundTextColor.get(i));
            text.setSize(text.getPreferredSize());
            
            int stringW = text.getWidth();
            double x = backgroundTextX.get(i);
            double y = backgroundTextY.get(i);
            
            if (x == LEFT) {
                x = ins.left + dx;
            } else if (x == CENTER) {
                x = ins.left + (c.getWidth() - ins.left - ins.right) / 2 - stringW / 2;
            } else if (x == RIGHT) {
                x = c.getWidth() - ins.right - stringW - dx;
            }
            
            if (y == TOP) {
                y = ins.top + dy;
            } else if (y == CENTER) {
                y = ins.top + (c.getHeight() - ins.top - ins.bottom) / 2 -
                text.getHeight() / 2;
            } else if (y == RIGHT) {
                y = c.getHeight() - ins.bottom - text.getHeight() - dy;
            }
            
            g.translate(x, y);
            text.paint(g);
            g.translate( -x, -y);
        }
    }

}
