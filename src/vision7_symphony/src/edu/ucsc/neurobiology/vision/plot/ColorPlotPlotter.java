package edu.ucsc.neurobiology.vision.plot;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.Raster;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

import edu.ucsc.neurobiology.vision.neuronviewer.STAPlotMaker;
import edu.ucsc.neurobiology.vision.plot.ColorPlotStyle.ColorStyle;
import edu.ucsc.neurobiology.vision.util.ParallelUtil;
import edu.ucsc.neurobiology.vision.util.ParallelUtil.StartEndIndices;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ColorPlotPlotter implements DataPlotter {
    private int drawThreads = 1;
    private final BlockingQueue<DrawThread> processingQueue;
    
    public ColorPlotPlotter() {
        if (STAPlotMaker.drawThreads > 0) this.drawThreads = STAPlotMaker.drawThreads;
        this.processingQueue = new ArrayBlockingQueue<DrawThread>(drawThreads);
    }
    
    public synchronized boolean accept(Object data, Object style) {
        if ( (data instanceof ColorPlotData) && (style instanceof ColorPlotStyle)) {
            return true;
        } else {
            return false;
        }
    }
    
    public synchronized void draw(Axis horAxis, Axis verAxis, Graphics graphics, Raster raster, 
            Object _data, Object _style) {
        long startTime = System.currentTimeMillis(); // For profiling if activated; see end of method
        
        ColorPlotData data = (ColorPlotData) _data;
        ColorStyle style = ((ColorPlotStyle) _style).getColorStyle();

        double dx = (data.getMaxX() - data.getMinX()) / data.getColumnsCount();
        double dy = (data.getMaxY() - data.getMinY()) / data.getRowsCount();
        
        int nThreads = drawThreads;
        DrawThread[] drawThreads = new DrawThread[nThreads]; 
        for (int i = 0; i < nThreads; i++) {
            drawThreads[i] = new DrawThread(data, style, horAxis, verAxis, dx, dy, i, nThreads, this);
            drawThreads[i].start();
        }
        
        int threadsProcessed = 0;
        while (threadsProcessed < nThreads) {
            try {
                processingQueue.take().process((Graphics2D) graphics);
            } catch (InterruptedException e) {return;}
            threadsProcessed++;
        }
        
        // Old single threaded version; may still pull this in to do one portion in the main thread
//      float[] colors = new float[3];
//      for (int i = 0; i < data.getColumnsCount(); i++) {
//          for (int j = 0; j < data.getRowsCount(); j++) {
//              colors = data.getCell(i, j, colors);
//              if ( (colors[0] < 0) || (colors[0] > 1) ||
//                  (colors[1] < 0) || (colors[1] > 1) ||
//                  (colors[2] < 0) || (colors[2] > 1)) {
//                  System.out.println("problem at " + i + ", " + j);
//                  System.out.println("" + colors[0] + ":" + colors[1] + ":" +
//                                     colors[2]);
//              }
//              float f;
//              Color c;
//
//              switch (style) {
//                  case RGB:
//                      c = new Color(colors[0], colors[1], colors[2]);
//                      break;
//                  case R:
//                      c = new Color(colors[0], colors[0], colors[0]);
//                      break;
//                  case G:
//                      c = new Color(colors[1], colors[1], colors[1]);
//                      break;
//                  case B:
//                      c = new Color(colors[2], colors[2], colors[2]);
//                      break;
//                  case R_AND_G:
//                      f = 0.5f * (colors[0] + colors[1]);
//                      c = new Color(f, f, f);
//                      break;
//                  case R_MINUS_G:
//                      f = 0.5f * (colors[0] - colors[1] + 1.f);
//                      c = new Color(f, f, f);
//                      break;
//                  default:
//                      throw new IllegalArgumentException(
//                          "Illegal Color Separation " + style);
//              }
//
//              double x1 = horAxis.getScreenCoord(data.getMinX() + i * dx);
//              double y1 = verAxis.getScreenCoord(data.getMaxY() - j * dy);
//              double x2 = horAxis.getScreenCoord(data.getMinX() + (i + 1) * dx);
//              double y2 = verAxis.getScreenCoord(data.getMaxY() - (j + 1) * dy);
//
//              graphics.setColor(c);
//              ( (Graphics2D) graphics).fill(
//                  new Rectangle2D.Double(x1, y1, x2 - x1, y2 - y1));
//          }
//      }
        
      graphics.setColor(Color.white);
      if (STAPlotMaker.profile) System.out.println("ColorPlorPlotter#draw: " + (System.currentTimeMillis() - startTime) + " ms.\n");
    }

    public void queueForProcessing(DrawThread toQueue) {
        try {
            processingQueue.put(toQueue);
        } catch (InterruptedException e) {return;}
    }

    class DrawThread extends Thread {
        private ColorPlotData  data;
        private ColorStyle style;
        private final Axis horAxis, verAxis;
        private final double dx, dy;
        private final StartEndIndices myCols;
        private final int numCols;
        private final ColorPlotPlotter plotter;
        private Color[][]       colors;
        private Rectangle2D[][] boxes;
        
        public DrawThread(ColorPlotData data, ColorStyle style, Axis horAxis, Axis verAxis, double dx, double dy, int thisThread, int nThreads, ColorPlotPlotter plotter) {
            this.data  = data;
            this.style = style;
            this.horAxis = horAxis;
            this.verAxis = verAxis;
            this.dx = dx;
            this.dy = dy;
            this.plotter = plotter;
            
            myCols = ParallelUtil.subIndex(thisThread, nThreads, data.getColumnsCount());
            numCols = myCols.end - myCols.start + 1;
            colors = new Color[numCols][data.getRowsCount()];
            boxes  = new Rectangle2D[numCols][data.getRowsCount()];
        }
        
        public void run() {
            float[] rgb = new float[3];
            float f;
            for (int i = 0; i < numCols; i++) {
                for (int j = 0; j < data.getRowsCount(); j++) {
                    
                    // FIXME: This is not technically thread-safe; it is possible for the underlying STA to be 
                    // changed in the process of being plotted.  I'm not sure if in practice there is ever another 
                    // thread going that would be doing this, but...  This probably shouldn't cause crashes as the 
                    // ways the STA could change are at least atomic (I think).  There may be issues with happens-before 
                    // relationships though.  Need to look at this more and see if it's possible to synchronize
                    // properly without losing performance.  Should be possible to get things right (but complicated) 
                    // using a lock shared by drawThreads.
                    rgb = data.getCell(i+myCols.start, j, rgb);
                    
                    switch (style) {
                        case RGB:
                            colors[i][j] = new Color(rgb[0], rgb[1], rgb[2]);
                            break;
                        case R:
                            colors[i][j] = new Color(rgb[0], rgb[0], rgb[0]);
                            break;
                        case G:
                            colors[i][j] = new Color(rgb[1], rgb[1], rgb[1]);
                            break;
                        case B:
                            colors[i][j] = new Color(rgb[2], rgb[2], rgb[2]);
                            break;
                        case R_AND_G:
                            f = 0.5f * (rgb[0] + rgb[1]);
                            colors[i][j] = new Color(f, f, f);
                            break;
                        case R_MINUS_G:
                            f = 0.5f * (rgb[0] - rgb[1] + 1.f);
                            colors[i][j] = new Color(f, f, f);
                            break;
                        default:
                            throw new IllegalArgumentException(
                                "Illegal Color Separation " + style);
                    }

                    double x1 = horAxis.getScreenCoord(data.getMinX() + (i + myCols.start)     * dx);
                    double y1 = verAxis.getScreenCoord(data.getMaxY() - j                  * dy);
                    double x2 = horAxis.getScreenCoord(data.getMinX() + (i + myCols.start + 1) * dx);
                    double y2 = verAxis.getScreenCoord(data.getMaxY() - (j + 1)            * dy);			
                    boxes[i][j] = new Rectangle2D.Double(x1, y1, x2 - x1, y2 - y1);
                }
            }
            
            plotter.queueForProcessing(this);
        }
        
        public void process(Graphics2D graphics) {
            for (int i = 0; i < numCols; i++) {
                for (int j = 0; j < data.getRowsCount(); j++) {
                    graphics.setColor(colors[i][j]);
                    graphics.fill(boxes[i][j]);
                }
            }
        }
        
    }
    
}