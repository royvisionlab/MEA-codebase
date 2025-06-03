package edu.ucsc.neurobiology.vision.test;

import java.io.*;
import java.text.*;

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import javax.swing.*;
import javax.swing.event.*;

import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.plot.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ImagingAnimation
    extends JPanel implements ScaledMouseListener {

    final int nElectrodes;
    final double[] amplitude;
    final ElectrodeMap map;
    final double minAmplitude;
    final double maxAmplitude;
    final float[][][] image;
    final PlotPanel p;
    final ScatterPlotStyle style;


    public ImagingAnimation(final float[][][] image) {
        super(new BorderLayout());
        this.image = image;

        nElectrodes = image[0].length;
        amplitude = new double[nElectrodes];

        for (int electrode = 0; electrode < nElectrodes; electrode++) {
            amplitude[electrode] = -MathUtil.min(image[0][electrode]);
        }
        minAmplitude = MathUtil.min(amplitude);
        maxAmplitude = MathUtil.max(amplitude);

        map = ElectrodeMapFactory.getElectrodeMap(501);

        style = new ScatterPlotStyle();
        style.setSymbolType(SymbolType.NONE);
        style.setErrorSymbolType(SymbolType.CIRCLE);
        p = new PlotPanel();

        final JSpinner spin = new JSpinner(new SpinnerNumberModel(0, -1,
            image[0][0].length - 1, 1));
        spin.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                drawImage( (Integer) spin.getValue());
            }
        });

        spin.setValue(new Integer( -1));

        add(p, BorderLayout.CENTER);
        add(spin, BorderLayout.NORTH);
    }


    public void clickPerformed(Component source, MouseEvent e, double x, double y) {
        Point2D.Double p = new Point2D.Double(0, 0);
        double d2Min = Double.POSITIVE_INFINITY;
        int electrode = -1;
        for (int i = 0; i < map.getNumberOfElectrodes(); i++) {
            map.getPosition(i, p);
            double d2 = (x - p.x) * (x - p.x) + (y - p.y) * (y - p.y);
            if (d2 < d2Min) {
                d2Min = d2;
                electrode = i;
            }
        }
        System.out.println("Electrode " + electrode);
        ScatterPlotStyle style = new ScatterPlotStyle();
        style.setConnectingPoints(true);
        PlotUtil.showData("" + electrode, new ScatterPlot(
            null, image[0][electrode], image[1][electrode], null), style);
    }


    private void drawImage(int i) {
        ScatterPlot sp = new ScatterPlot("Immaging");
        if (i == -1) {
            for (int electrode = 0; electrode < nElectrodes; electrode++) {
                Point2D x = map.getPosition(electrode, null);
                double r = 30.0 * (amplitude[electrode] - minAmplitude) /
                           (maxAmplitude - minAmplitude);
                sp.add(x.getX(), x.getY(), r);
            }
        } else {
            for (int electrode = 0; electrode < nElectrodes; electrode++) {
                Point2D x = map.getPosition(electrode, null);
                double r = 30.0 * ( -image[0][electrode][i] - minAmplitude) /
                           (maxAmplitude - minAmplitude);
                sp.add(x.getX(), x.getY(), r);
            }
        }

        p.removeAllData();
        p.addData(sp, style);
        p.setAxisVisible(false);
        p.setRange( -1000, 1000, -500, 500);
        p.repaint();
    }


    public void saveFrames(String dir) throws IOException {
        DecimalFormat formatter = new DecimalFormat();
        formatter.setMaximumFractionDigits(0);
        formatter.setMinimumIntegerDigits(5);
        formatter.setGroupingUsed(false);

        File d = new File(dir);
        d.mkdir();
        String fileSeparator = System.getProperty("file.separator");

        p.setXAxisVisible(false);
        p.setYAxisVisible(false);
        p.setSize(600, 300);

        for (int i = 0; i < image[0][0].length - 1; i++) {
            System.out.println("Saving frame " + i);
            drawImage(i);
            File f = new File(
                dir + fileSeparator + "frame" + formatter.format(i) +
                ".png");
//            new PNGExportFileType().exportToFile(f, p, null, null, "dumitru");
        }
    }


    public static void main(String[] args) throws Exception {
        PhysiologicalImagingFile imgFile = new PhysiologicalImagingFile(
            "d:\\data\\2003-09-19-3\\data001\\data001.ei");
        ImagingAnimation p = new ImagingAnimation(imgFile.getImage(1984));

//        p.saveFrames("177-img");

        JFrame f = new JFrame();
        f.add(p);
        f.setBounds(100, 100, 800, 400);
        f.setVisible(true);
    }


    public void enteredPerformed(Component source, MouseEvent event, double x,
            double y) {
    }


    public void exitedPerformed(Component source, MouseEvent event, double x,
            double y) {	
    }


    public void pressPerformed(Component source, MouseEvent event, double x,
            double y) {

    }


    public void releasePerformed(Component source, MouseEvent event, double x,
            double y) {

    }
}
