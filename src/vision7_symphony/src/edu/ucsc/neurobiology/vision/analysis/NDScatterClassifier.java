package edu.ucsc.neurobiology.vision.analysis;

import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import java.util.ArrayList;

import javax.swing.*;
import javax.swing.tree.*;

import edu.ucsc.neurobiology.vision.gui.GraphicsIO;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.neuronviewer.*;
import edu.ucsc.neurobiology.vision.plot.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Matthew Grivich, University of California, Santa Cruz
 */
public class NDScatterClassifier
extends JPanel {

    boolean showAxes = true;
    static int classIndex = 0;
    int[] neurons; //list of the neurons in plot - stored as Integer objects
    int[] classNumber;
    Color[] colors = new Color[]{Color.black, Color.magenta, Color.red, Color.blue, Color.green};
    double[][] data; //[component][eventNumber]
    double[][] rotatedData; //[component][eventNumber]
    double[][] basisVectors; //[component][basisNumber]  == fullRotationMatrix
    double[][] newBasisVectors;
    int axis1 = 0; //rotation axis 1
    int axis2 = 1; //rotation axis 2
    int setAxis = 0; //which axis to set next
    String[] dimensionNames;

    int nDims;
    int nNeurons;

    double screenWidth = 800.0; //width of jPanel becomes 792 due to border
    double screenHeight = 800.0; //height of jPanel becomes 773 due to border
    final double dataWidth = 3.5;
    final double dataHeight = 3.5;

    //coordinates for selection box
    double x1 = 0.0;
    double y1 = 0.0;
    double x2 = 0.0;
    double y2 = 0.0;
    int mouseClick = 0;

    PlotPanel pan;
    private final ClassifierListener classifierListener;
    TreePath classPath;

    final JPopupMenu menu = new JPopupMenu();



    public NDScatterClassifier(ClassifierListener classifierListener,
            TreePath classPath, int[] neurons, String[] classNames, double dataIn[][],
            String[] dimensionNames) {
        this.classifierListener = classifierListener;
        this.classPath = classPath;


        this.dimensionNames = dimensionNames;
        this.neurons = neurons;

        this.nNeurons = dataIn[0].length;
        this.nDims = dataIn.length;

        this.classNumber = new int[nNeurons];

//		sort neurons into class, by number, so that they can be plotted in color by class.
        ArrayList<String> classesList = new ArrayList<String>();
        boolean found = false;
        for(int i=0; i<nNeurons; i+=1) {
            for(int j=0; j<classesList.size(); j++) {
                if(classNames[i].equals(classesList.get(j))) {
                    found = true;
                    classNumber[i] = j;
                }
            }
            if(!found) {
                classesList.add(classNames[i]);
                classNumber[i] = classesList.size()-1;
            }
            found = false;
        }

//		Rescale data and paste into rotated data array.
//		System.out.println(nNeurons + "    " + nDims);
        data = new double[nDims][nNeurons];
        rotatedData = new double[nDims][nNeurons];

        for (int j = 0; j < nDims; j++) {
            double max = MathUtil.max(dataIn[j]);
            double min = MathUtil.min(dataIn[j]);
            max = Math.max(Math.max(Math.abs(max), Math.abs(min)), .0000001);
//			System.out.println(j + "    " + max);

            for (int i = 0; i < nNeurons; i++) {
                if (!Double.isNaN(dataIn[j][i])) {
                    rotatedData[j][i] = data[j][i] = dataIn[j][i] / max;
                }
//				System.out.println(j + "    " + data[j][i]);
            }

        }

        basisVectors = new double[nDims][nDims];
        newBasisVectors = new double[nDims][nDims];

        for (int i = 0; i < nDims; i++) {
            for (int j = 0; j < nDims; j++) {
                if (i == j) {
                    basisVectors[i][j] = 1.0;
                } else {
                    basisVectors[i][j] = 0.0;
                }
            }
        }

        menu.add(new AbstractAction("Save Image...") {
            public void actionPerformed(ActionEvent event) {
                GraphicsIO.saveImage(NDScatterClassifier.this, NDScatterClassifier.this, "Save Image");
            }});

        run();
    }


    public void run() {
        pan = new PlotPanel();

        JFrame f = new JFrame("NDScatterPlot");

        f.addKeyListener(new KeyAdapter() {
            public void keyPressed(KeyEvent e) {
                if (axis1 != axis2) {
                    if (e.getKeyCode() == KeyEvent.VK_LEFT) {
                        rotate(Math.PI / 8);
                        repaint();
                    } else if (e.getKeyCode() == KeyEvent.VK_RIGHT) {
                        rotate( -Math.PI / 8);
                        repaint();
                    }
                }
                if (e.getKeyCode() > 47 && e.getKeyCode() <= 47 + nDims) {
                    e.setKeyCode(e.getKeyCode() + 48);
                }
                if (e.getKeyCode() > 95 && e.getKeyCode() <= 95 + nDims) {
                    if (setAxis == 0) {
                        axis1 = e.getKeyCode() - 96;
                        setAxis = 1;
                        repaint();
//						repaint(0, 0, 200, 50);
                    } else {
                        axis2 = e.getKeyCode() - 96;
                        setAxis = 0;
                        repaint();
//						repaint(0, 0, 200, 50);
                    }
                }

                if (e.getKeyCode() == KeyEvent.VK_INSERT) {
                    showAxes = !showAxes;
                    repaint();
                }
            }
        });

        // add listener to jPanel so that coordinates come out correctly
        addMouseListener(new MouseAdapter() {
            public void mousePressed(MouseEvent e) {


                // one and two button mouse friendly.
                if (SwingUtilities.isRightMouseButton(e) || 
                        (e.getModifiersEx() & MouseEvent.ALT_DOWN_MASK) == MouseEvent.ALT_DOWN_MASK) {
                    menu.show(NDScatterClassifier.this, e.getX(), e.getY());
                    return;
                }


                if (mouseClick == 0) {
                    x1 = ( (double) e.getX() - screenWidth / 2.0) * dataWidth /
                    screenWidth;
                    y1 = ( (double) screenHeight / 2.0 - (double) e.getY()) *
                    dataHeight / screenHeight;
                    mouseClick = 1;

                } else {
                    mouseClick = 0;
                    x2 = ( (double) e.getX() - screenWidth / 2.0) * dataWidth /
                    screenWidth;
                    y2 = (screenHeight / 2.0 - (double) e.getY()) * dataHeight /
                    screenHeight;
                    double temp = x1;
                    x1 = Math.min(x2, temp);
                    x2 = Math.max(x2, temp);

                    temp = y1;
                    y1 = Math.min(y2, temp);
                    y2 = Math.max(y2, temp);

                    //                  System.out.println(x1 + "   " + x2 + "   " + y1 + "   " + y2);
                    repaint();
                    IntegerList selectedNeurons = new IntegerList();
                    for (int neuron = 0; neuron < nNeurons; neuron++) {
                        if (rotatedData[0][neuron] > x1 && rotatedData[0][neuron] < x2 &&
                                rotatedData[1][neuron] > y1 && rotatedData[1][neuron] < y2) {
                            selectedNeurons.add(neurons[neuron]);
                            //                              System.out.println(neurons.get(neuron));
                        } // if
                    }

                    classifierListener.classificationHappened(
                            selectedNeurons.toArray(),
                            classPath.pathByAddingChild(new DefaultMutableTreeNode("ND-New-" +
                                    classIndex, true)));
                    classIndex++;
                }
            }
        });

        f.setBackground(Color.white);
        f.add("Center", this);
        f.pack();
        f.setSize(new Dimension( (int) screenWidth, (int) screenHeight));
        f.setVisible(true);
    }






    public void paint(Graphics g) {

        Graphics2D g2 = (Graphics2D) g;
        g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
                RenderingHints.VALUE_ANTIALIAS_ON);
        Dimension d = getSize();
        screenWidth = d.width;
        screenHeight = d.height;
        g2.setColor(Color.white);
        g2.fillRect(0, 0, (int) screenWidth, (int) screenHeight);
        g2.setColor(Color.black);

        if (mouseClick == 0) {
            g2.draw(new Rectangle2D.Double(
                    x1 * screenWidth / dataWidth + screenWidth / 2.0,
                    screenHeight / 2.0 - y2 * screenHeight / dataHeight,
                    (x2 - x1) * screenWidth / dataWidth,
                    (y2 - y1) * screenHeight / dataHeight));

            // draw data
        }
        for (int i = 0; i < nNeurons; i++) {
            g2.setColor(colors[classNumber[i]%colors.length]);
            double x = rotatedData[0][i] * screenWidth / dataWidth + screenWidth / 2;
            double y = -rotatedData[1][i] * screenHeight / dataHeight + screenHeight / 2;
            if (!Double.isNaN(x) && !Double.isNaN(y)) {
                g2.drawRect( (int) x, (int) y, 3, 3);
            }
        }

        g2.setColor(Color.black);

        //draw basis vectors
        for (int i = 0; i < nDims; i++) {
            if (showAxes) {
                g2.drawLine(
                        (int) (screenWidth / 2), (int) (screenHeight / 2),
                        (int) (screenWidth / 2 + basisVectors[0][i] * screenWidth / 3),
                        (int) (screenHeight / 2 - basisVectors[1][i] * screenHeight / 3));
            }

            g2.drawString(
                    (new Integer(i)).toString(),
                    (int) (screenWidth / 2 + basisVectors[0][i] * screenWidth / 3),
                    (int) (screenHeight / 2 - basisVectors[1][i] * screenHeight / 3));
        }

        g2.setColor(Color.gray);
        g2.drawString("Rotation axes: " + axis1 + " " + axis2, 10, 20);
        for (int i = 0; i < dimensionNames.length; i++) {
            g2.drawString(i + ": " + dimensionNames[i], 10, 40 + 20 * i);
        }
        g2.setColor(Color.black);

    }


    public void rotate(double angle) {
        double rotationMatrix[][] = new double[nDims][nDims];

        double temp = 0;

        for (int i = 0; i < nDims; i++) {
            for (int j = 0; j < nDims; j++) {
                rotationMatrix[i][j] = Math.cos(angle) * basisVectors[i][axis1] *
                basisVectors[j][axis1]
                                +
                                Math.cos(angle) * basisVectors[i][axis2] *
                                basisVectors[j][axis2]
                                                -
                                                Math.sin(angle) * basisVectors[i][axis1] *
                                                basisVectors[j][axis2]
                                                                +
                                                                Math.sin(angle) * basisVectors[i][axis2] *
                                                                basisVectors[j][axis1];

                for (int k = 0; k < nDims; k++) {
                    if (axis1 != k && axis2 != k) {
                        rotationMatrix[i][j] += basisVectors[i][k] * basisVectors[j][k];
                    }
                }
            }
        }

        //determine newBasisVectors==totalRotationMatrix
        for (int i = 0; i < nDims; i++) {
            for (int j = 0; j < nDims; j++) {
                newBasisVectors[i][j] = 0;
            }

        }

        for (int i = 0; i < nDims; i++) {
            for (int j = 0; j < nDims; j++) {
                for (int k = 0; k < nDims; k++) {
                    newBasisVectors[i][j] += rotationMatrix[i][k] * basisVectors[k][j];
                }
            }
        }

        for (int i = 0; i < nDims; i++) {
            for (int j = 0; j < nDims; j++) {
                basisVectors[i][j] = newBasisVectors[i][j];
            }
        }

        //orthonormalize basisVectors == rotationMatrix, "fixes" rounding error
        //using Gram-Schmidt Process
        //for each basis vector

        for (int i = 0; i < nDims; i++) {
            //subtract off components parallel to other vectors
            for (int j = 0; j < i; j++) {
                double parallelComponent = 0.0;
                for (int k = 0; k < nDims; k++) {
                    parallelComponent += basisVectors[k][j] * basisVectors[k][i];
                }
                for (int k = 0; k < nDims; k++) {
                    basisVectors[k][i] = basisVectors[k][i] -
                    parallelComponent * basisVectors[k][j];
                }
            }

            //normalize basis vector
            double magnitude = 0.0;
            for (int j = 0; j < nDims; j++) {
                magnitude += basisVectors[j][i] * basisVectors[j][i];
            }
            magnitude = Math.sqrt(magnitude);
            for (int j = 0; j < nDims; j++) {
                basisVectors[j][i] /= magnitude;
            }
        }

        //rotate Data
        for (int i = 0; i < nDims; i++) {
            for (int j = 0; j < nNeurons; j++) {
                rotatedData[i][j] = 0;
            }
        }

        for (int i = 0; i < nNeurons; i++) {
            for (int row = 0; row < nDims; row++) {
                for (int column = 0; column < nDims; column++) {
                    rotatedData[row][i] += basisVectors[row][column] * data[column][i];

                }
            }
        }
    }





}
