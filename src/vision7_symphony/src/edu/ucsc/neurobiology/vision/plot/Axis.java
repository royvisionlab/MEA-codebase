package edu.ucsc.neurobiology.vision.plot;

import static java.lang.Math.abs;
import static java.lang.Math.ceil;
import static java.lang.Math.floor;
import static java.lang.Math.log10;
import static java.lang.Math.pow;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.text.DecimalFormat;
import java.util.LinkedHashMap;

import javax.swing.JLabel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;


/**
 * Abstract implementation of a plot Axis. Defines most of the axis functionality but
 * delegates the following methods to concrete classes.
 <pre>
    public abstract void paintAxis(Graphics2D g, int x0, int y0, int width, int height);
    public abstract double getScreenCoord(double x);
    public abstract double getPlotCoord(double s);
    public abstract int updateTicks();
    public abstract double getK();
 </pre>
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public abstract class Axis {

    public int maxTickMarkSize;
    boolean DEBUG = false;

    protected AxisType type = AxisType.LINEAR;

    protected double x1, x2;
    protected int s2, s1;
    protected LinkedHashMap<Double, String> ticksAndLabels;
    protected int minTickSeparation = 20;

    protected int tickSize = 3;
    protected int tickToMarkSpacing = 3;
    protected int markToLabelSpacing = 3;

    protected int gridPadding = 10;
    protected DecimalFormat df;
    protected Font currentFont;
    protected FontMetrics metrics;
    protected JLabel labelComponent;
    protected boolean showGridLines = false;
    protected ChangeListener listener;
    protected double secondaryK, secondaryB;
    protected String secondaryLabel = null;

    protected boolean showTicks = true;
//	protected boolean showTickMarks = true;
    protected boolean showLabel = true;
    public static boolean shiftRightLabel = false;
    double f = 1.2;
    int fixedThickness = -1;
    int fixedLabelPosition = -1;
    protected double fixedTick = Double.NaN;
    protected double fixedTickSpacing = Double.NaN;
    protected int fixedFractionalDigits = -1;
    public int tickMarkAllignment = HorizontalAxis.CENTER_TICK_ALLIGNMENT;
    public double multiplier = 1;

    protected double[] manualTicks = null;


    public void setSpacings(int markToLabelSpacing, int tickToMarkSpacing) {
        this.tickToMarkSpacing = tickToMarkSpacing;
        this.markToLabelSpacing = markToLabelSpacing;
        updateTicks();
    }


    /**
     * Use this function to override automatic tick spacing behavior
     * @author Matthew Grivich
     *
     * @param fixedMinimumTick double
     * @param fixedTickSpacing double
     * @param fixedMaximumFractionalDigits int
     */
    public void setFixedTickSpacing(double fixedMinimumTick, double fixedTickSpacing,
            int fixedMaximumFractionalDigits) {
        this.fixedTick = fixedMinimumTick;
        this.fixedTickSpacing = fixedTickSpacing;
        this.fixedFractionalDigits = fixedMaximumFractionalDigits;
    }

    
    /**
     * Set the number of ticks, by determining the correct fixedTickSpacing for the 
     * current axis range; only works for fixed linear ticks.  Rounds off to nearest
     * tens place, with slight bias to lower powers.
     * 
     * @author Peter H. Li
     * @param num
     */
    public void setNumRoundedFixedLinearTicks(int num) {
        fixedTickSpacing = Math.abs(x2 - x1) / num;
        roundTickSpacing();
    }

    
    public void roundTickSpacing() {
        int power = AxisUtil.getPowerOf10(fixedTickSpacing / 2); // Bias towards smaller powers
        fixedTickSpacing = AxisUtil.round(fixedTickSpacing, power);
    }
    

    public void setManualTicks(double[] ticks) {
        manualTicks = ticks;
    }

    public void clearManualTicks() {
        manualTicks = null;
    }


    public void setMarkToLabelSpacing(int markToLabelSpacing) {
        this.markToLabelSpacing = markToLabelSpacing;
        updateTicks();
    }


    public int getmarkToLabelSpacing() {
        return markToLabelSpacing;
    }


    public void setTickToMarkSpacing(int tickToMarkSpacing) {
        this.tickToMarkSpacing = tickToMarkSpacing;
        updateTicks();
    }


    public int getTickToMarkSpacing() {
        return tickToMarkSpacing;
    }


    public void setTickSize(int tickSize) {
        this.tickSize = tickSize;
        updateTicks();
    }


    public int getTickSize() {
        return tickSize;
    }


    public void setShowTicks(boolean showTicks) {
        this.showTicks = showTicks;
        updateTicks();
    }


//	public void setShowTickMarks(boolean showTickMarks) {
//	this.showTickMarks = showTickMarks;
//	updateTicks();
//	}


    public void setShowlabel(boolean showLabel) {
        this.showLabel = showLabel;
        updateTicks();
    }


    public boolean getShowTicks() {
        return showTicks;
    }


    public boolean getShowlabel() {
        return showLabel;
    }


    public Axis() {
        ticksAndLabels = new LinkedHashMap<Double, String> ();

        df = new DecimalFormat();
        df.setMinimumFractionDigits(0);
        df.setMinimumIntegerDigits(0);

        labelComponent = new JLabel("", JLabel.CENTER);
        labelComponent.setForeground(Color.black);
//		labelComponent.setBorder(BorderFactory.createLineBorder(Color.black,1));

        setLabelFont(new Font("Arial", Font.ITALIC, 12));
    }


    public void setAxisType(AxisType type) {
        this.type = type;
        updateTicks();
    }


    public AxisType getAxisType() {
        return type;
    }


    protected final double convert(double x) {
        switch (type) {
        case LINEAR:
            return x;
        default:
            return log10(x);
        }
    }


    public void setChangeListener(ChangeListener listener) {
        this.listener = listener;
    }


    public void setLabel(String label) {
        labelComponent.setText(label);
        labelComponent.validate();

        Dimension d = labelComponent.getPreferredSize();
        labelComponent.setSize( (int) (d.width * f), d.height);
        updateTicks();
    }


    public String getLabel() {
        return labelComponent.getText();
    }


    public void setScreenRange(int s1, int s2) {
        if (this.s1 != s1 || this.s2 != s2) {
            this.s1 = s1;
            this.s2 = s2;

            updateTicks();

            if (listener != null) {
                listener.stateChanged(new ChangeEvent(this));
            }
        }
    }


    public void setPlotRange(double x1, double x2) {
        this.x1 = x1;
        this.x2 = x2;
        updateTicks();

        if (listener != null) {
            listener.stateChanged(new ChangeEvent(this));
        }
    }


    public void setShowGridLines(boolean showGridLines) {
        this.showGridLines = showGridLines;
    }


    public boolean isShowingGridLines() {
        return showGridLines;
    }


    public double[] getPlotRange() {
        return new double[] {x1, x2};
    }


    public int[] getScreenRange() {
        return new int[] {s1, s2};
    }


    protected synchronized void createLinearFixedTicks() {
        ticksAndLabels.clear();

        int nFractDigits = fixedFractionalDigits == -1 ?
                AxisUtil.getFractionalDigits(multiplier * fixedTickSpacing) :
                fixedFractionalDigits;

        DecimalFormat df = new DecimalFormat();
        df.setMaximumFractionDigits(nFractDigits);
        df.setMaximumIntegerDigits(Integer.MAX_VALUE);

        for (double tick = fixedTick; tick <= x2; tick += fixedTickSpacing) {
            String label = df.format(multiplier * tick);
            ticksAndLabels.put(tick, label);
        }

        for (double tick = fixedTick - fixedTickSpacing; tick >= x1;
        tick -= fixedTickSpacing) {
            String label = df.format(multiplier * tick);
            ticksAndLabels.put(tick, label);
        }
    }

    //only works for linear axis.
    protected synchronized void createManualTicks() {

        if (manualTicks != null) {
            ticksAndLabels.clear();
            //find minimum spacing
            double spacing = Double.MAX_VALUE;
            for (int i = 0; i < manualTicks.length-1; i++) {
                if (manualTicks[i+1] - manualTicks[i] < spacing) {
                    
                    spacing = manualTicks[i+1] - manualTicks[i];
                }
            }
            int digits = 0;
            if (spacing >= 10) {
                digits = 0;
            } else {
                //get first non zero fractional digit

                for (int i = 0; i < Integer.MAX_VALUE; i++, spacing *= 10) {
                    if (spacing > 1) {
                        digits = i+1;
                        break;
                    }
                }
            }

            int nFractDigits = fixedFractionalDigits == -1 ?
                    digits : fixedFractionalDigits;

            DecimalFormat df = new DecimalFormat();
            df.setMaximumFractionDigits(nFractDigits);
            df.setMaximumIntegerDigits(Integer.MAX_VALUE);
            for (int i = 0; i < manualTicks.length; i++) {
                String label = df.format(manualTicks[i]);
                ticksAndLabels.put(manualTicks[i], label);
            }
        }
    }


    protected synchronized void createLinearTicks(int factor) {
        ticksAndLabels.clear();
        double[] autoScaledX = AxisUtil.performAutoScale(x1, x2);
        double dx = AxisUtil.calculateTickSpacing(autoScaledX[0], autoScaledX[1]);
        int nFractDigits = fixedFractionalDigits == -1 ?
                AxisUtil.getFractionalDigits(dx) : fixedFractionalDigits;
        
        dx *= factor;
        int nTicks = 2 + (int) ( (autoScaledX[1] - autoScaledX[0]) / dx);

        DecimalFormat df = new DecimalFormat();
        df.setMaximumFractionDigits(nFractDigits);
        df.setMaximumIntegerDigits(Integer.MAX_VALUE);

        for (int i = 0; i < nTicks; i++) {
            double mainTick = autoScaledX[0] + i * dx;
            // 10
            for (int j = 0; j < 1; j++) {
                double tick = mainTick + j * dx / 10;
                if (tick >= x1 && tick <= x2) {
                    String label;
                    if (j == 0) {
                        if (secondaryLabel == null) {
                            label = df.format(tick);
                        } else {
                            double secondaryXi = secondaryK * tick + secondaryB;
                            label = df.format(tick) + "(" + df.format(secondaryXi) + ")";
                        }
                    } else {
                        label = null;
                    }
                    ticksAndLabels.put(tick, label);
                }
            }
        }
    }


    protected synchronized void createLogTicks(int factor) {
        ticksAndLabels.clear();

        if (x1 == 0) {
            x1 = Math.min(0.01, x2 / 10);
        }

        double x1Log = log10(x1);
        double x2Log = log10(x2);

        int x1pow = (int) ceil(x1Log) - 1;
        int x2pow = (int) floor(x2Log);

        int nFractDigits = abs(x1pow) + 1;
        int nMainTicks = x2pow - x1pow + 1;

        DecimalFormat df = new DecimalFormat();
        df.setMaximumFractionDigits(nFractDigits);
        df.setMaximumIntegerDigits(Integer.MAX_VALUE);

//		System.err.println("LOG: " + x1Log + " -> " + x2Log);
//		System.err.println("POW: " + x1pow + " -> " + x2pow);
//		System.err.println("ticks: " + nMainTicks);

        for (int i = 0; i < nMainTicks; i++) {
            double mainTick = pow(10, x1pow + i);
            for (int j = 1; j <= 9; j++) {
                double tick = mainTick * j;
                if (tick >= x1 && tick <= x2) {
                    String label = null;
                    if (nMainTicks - 1 > 1) {
                        label = (j == 1) ? df.format(tick) : null;
                    } else {
                        label = df.format(tick);
                    }
                    ticksAndLabels.put(tick, label);
                }
            }
        }

    }


    public void setThickness(int fixedThickness) {
        this.fixedThickness = fixedThickness;
    }


    public void setLabelPosition(int fixedLabelPosition) {
        this.fixedLabelPosition = fixedLabelPosition;
    }


    public int getThickness() {
        if (fixedThickness > 0) {
            return fixedThickness;
        } else {
            return updateTicks();
        }
    }


    public void setSecondaryLabels(double secondaryK, double secondaryB,
            String secondaryLabel) {
        this.secondaryK = secondaryK;
        this.secondaryB = secondaryB;
        this.secondaryLabel = secondaryLabel;

        updateTicks();
    }


    /**
     * Must paint the axis.
     *
     * @param g Graphics2D
     * @param x0 int
     * @param y0 int
     * @param width int
     * @param height int
     */
    public abstract void paintAxis(Graphics2D g, int x0, int y0, int width, int height);


    /**
     * Must translate a real coordinate into a screen one.
     *
     * @param x double
     * @return double
     */
    public abstract double getScreenCoord(double x);


    /**
     * Must translate a screen coordinate into a real one.
     *
     * @param s double
     * @return double
     */
    public abstract double getPlotCoord(double s);


    /**
     * Must refresh the tick and tickmark info.
     *
     * @return int
     */
    public abstract int updateTicks();


    /**
     * Must return the scaling factor of the real->screen conversion.
     *
     * @return double
     */
    public abstract double getK();


    public void setLabelFont(Font font) {
        labelComponent.setFont(font);
        Dimension d = labelComponent.getPreferredSize();
        labelComponent.setSize( (int) (d.width * f), d.height);

        currentFont = font;
        BufferedImage image = new BufferedImage(10, 10, BufferedImage.TYPE_3BYTE_BGR);
        Graphics2D g2 = (Graphics2D) image.getGraphics();
        metrics = g2.getFontMetrics(currentFont);

//		descent = currentFont.getLineMetrics("H", new FontRenderContext(null, false, false)).getDescent()

        updateTicks();
    }


    public Font getLabelFont() {
        return currentFont;
    }

}
