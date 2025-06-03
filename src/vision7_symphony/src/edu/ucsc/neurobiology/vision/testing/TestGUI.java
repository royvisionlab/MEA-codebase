package edu.ucsc.neurobiology.vision.testing;

import java.io.*;
import java.util.*;

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.tree.*;

import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.gui.*;
import edu.ucsc.neurobiology.vision.gui.Desktop;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.neuronviewer.*;
import edu.ucsc.neurobiology.vision.plot.*;


/**
 * @author nobody, anyone can change
 */
public class TestGUI {

    public static void testPlotPanel() throws Exception {
        PlotPanel plot = new PlotPanel();

        double[] x = {
                     0.16837284482758622, 0.33674568965517243, 0.6734913793103449,
                     1.3469827586206897, 2.6939655172413794, 5.387931034482759,
                     10.775862068965518, 21.551724137931036, 43.10344827586207};
        double[] y = {
                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00338705573048595, 0.0};

        plot.addData(new ScatterPlot(x, y, null), new ScatterPlotStyle());
//        plot.loadStyle("f1f2");
//        plot.setAxesType(AxisType.LINEAR, AxisType.LINEAR);

 //       plot.addArrow(new Arrow(50, 50, 10, 0, 20, 15, 10, 1));
        plot.autoscale();

        plot.addToLegend("bbb");

        JFrame f = new JFrame("Axis Test");
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.add(plot);
        f.setBounds(100, 100, 400, 400);
        f.setVisible(true);
    }


    public static void testInteractiveTree() throws Exception {
        ParametersFile params = new ParametersFile("2003-08-06-2-data000.params");
        LinkedHashMap<Integer,? extends Object> c = params.getClassIDs();

        DefaultTreeModel m = InteractiveTree.makeModel(c);
        InteractiveTree t1 = new InteractiveTree(m);
        InteractiveTree t2 = new InteractiveTree(m);
        JFrame f = new JFrame();
        f.add(new JScrollPane(t1), BorderLayout.WEST);
        f.add(new JScrollPane(t2), BorderLayout.EAST);
        f.setBounds(100, 100, 500, 500);
        f.setVisible(true);

        t1.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
            }
        });
    }


    public static void testParametersTable() {
//        ParametersTable t = new ParametersTable();
//        final IntegerParameter p = new IntegerParameter("name", 5, 0, 50);
//        t.addParameter(p);
//
//        JButton b =  new JButton("get");
//        b.addActionListener(new ActionListener() {
//            public void actionPerformed(ActionEvent e) {
//                System.out.println(p.getValue());
//            }
//        });
//
//        JButton b1 =  new JButton("set");
//        b1.addActionListener(new ActionListener() {
//            public void actionPerformed(ActionEvent e) {
//                p.setValue(35);
//            }
//        });
//
//        // display a Frame with the Function
//        JFrame f = new JFrame("Function Example");
//        f.setDefaultCloseOperation(f.EXIT_ON_CLOSE);
//        f.setBounds(50, 50, 400, 400);
//        f.add(t, BorderLayout.CENTER);
//        JPanel panel = new JPanel(new GridLayout(2, 1));
//        panel.add(b);
//        panel.add(b1);
//        f.add(panel, BorderLayout.SOUTH);
//        f.setVisible(true);
    }


    public static void testImagingPanel() throws Exception {
        int id = 1670;
        String n = "f:\\data\\2005-04-26-0\\data009\\data009.ei";
        PhysiologicalImagingFile f = new PhysiologicalImagingFile(n);
        ElectrodeMap m = ElectrodeMapFactory.getElectrodeMap(f.arrayID);

        final PhysiologicalImagePanel p = new PhysiologicalImagePanel(
            f.getImage(id), null, 2, m, 0, null, "a");
//        p.addArrow(53, new PhysiologicalImagePanel.Arrow(0, 0, 5, -45, 50, 1f));

        final JFrame fr = new JFrame();
        fr.add(p);
        fr.setBounds(100, 100, 500, 300);
        fr.setVisible(true);
    }


    public static void testWindowBar() {
        try {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (Exception e) {
            e.printStackTrace();
        }

        JFrame f = new JFrame();
        f.setBounds(100, 100, 600, 60);

        WindowBar w = new WindowBar();
        w.addButton("Window 1");
        w.addButton("Window 2");
        w.addButton("Window 3");

        w.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                System.out.println("Selected " + e.getID());
            }
        });

        f.add(w);
        f.setVisible(true);
    }


    public static void testDataManager() throws IOException {
        JFrame f = new JFrame();
        f.setBounds(200, 200, 400, 800);
        DataManager manager = new DataManager();
        manager.setFolder("Y:\\2005-04-06-0");
        f.add(new JScrollPane(manager), BorderLayout.CENTER);
        f.setVisible(true);
    }


    public static void testDesktop() {
        try {
            UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
        } catch (Exception e) {
            e.printStackTrace();
        }

        JFrame f = new JFrame();
        f.setBounds(100, 100, 400, 400);

        Desktop w = new Desktop();
        JInternalFrame frame1 = new JInternalFrame();
        JInternalFrame frame2 = new JInternalFrame();
        JInternalFrame frame3 = new JInternalFrame();

        w.addFrame(frame1, "Window 1");
        w.addFrame(frame2, "Window 2");
        w.addFrame(frame3, "Window 3");

        f.add(w);
        f.setVisible(true);

        System.err.println(w.getAllFrames().length);
        w.removeFrame(frame1);
        System.err.println(w.getAllFrames().length);
    }


    public static void testHistogramPlotter() throws Exception {
        final DoubleHistogram h = new DoubleHistogram("", 0, 14400, 1);
        Random r = new Random(1234);
        for (int i = 0; i < h.getBinCount(); i++) {
            h.setBin(i, r.nextDouble());
        }
        final HistogramStyle s = new HistogramStyle();
        final PlotPanel p = PlotUtil.showData("", h, s);

        Thread.sleep(10);

        final HistogramPlotter plotter = new HistogramPlotter();
        final int n = 200;

        Thread t = new Thread() {
            public void run() {
                Graphics g = p.getGraphics();
                Insets ins = p.getInsets();
                int strokeW = 0;
                g.setClip(new Rectangle2D.Double(
                    ins.left + strokeW, ins.top + strokeW,
                    p.getWidth() - ins.left - ins.right - 2 * strokeW,
                    p.getHeight() - ins.top - ins.bottom - 2 * strokeW));

                long t1 = System.nanoTime();
                for (int i = 0; i < n; i++) {
                    plotter.draw(p.getXAxis(), p.getYAxis(), g, p.raster, h, s);
                }
                long t2 = System.nanoTime();
                System.err.println( (t2 - t1) / 1e9 / n);
            }
        };
        t.setName("Test Thread");

        t.start();
    }


    public static void main(String[] args) throws Exception {
        testPlotPanel();
    }
}
