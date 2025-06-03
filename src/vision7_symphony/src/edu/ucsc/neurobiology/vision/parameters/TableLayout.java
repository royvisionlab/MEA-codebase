package edu.ucsc.neurobiology.vision.parameters;

import java.awt.*;


/**
 * This class represents a two column table, each row containing one parameter.
 * The left column contains labels with parameters screen names, the right one
 * contains the component used to display and change the value of the parameter.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class TableLayout
    implements LayoutManager {
    final int hgap, vgap;
    int rows;
    final int cols;
    final int minimumWidth = 100;


    public TableLayout(int hgap, int vgap) {
        this.cols = 2;
        this.hgap = hgap;
        this.vgap = vgap;
    }


    public void addLayoutComponent(String name, Component comp) {
    }


    public void removeLayoutComponent(Component comp) {
    }


    public Dimension preferredLayoutSize(Container parent) {
        synchronized (parent.getTreeLock()) {
            Insets insets = parent.getInsets();
            int ncomponents = parent.getComponentCount();
            if (ncomponents == 0) {
                return new Dimension(0, 0);
            }

            rows = (ncomponents + cols - 1) / cols;

            int wMax = Integer.MIN_VALUE;
            int y = 0;
            for (int r = 0; r < rows; r++) {
                int maxH = Integer.MIN_VALUE;
                for (int c = 0; c < cols; c++) {
                    int i = r * cols + c;
                    Dimension prefSize = parent.getComponent(i).getPreferredSize();
                    if (c == 0) {
                        if (prefSize.width > wMax) {
                            wMax = prefSize.width;
                        }
                    }
                    if (i < ncomponents) {
                        int h = prefSize.height;
                        if (h > maxH) {
                            maxH = h;
                        }
                    }
                }

                y += maxH;
            }

//            return new Dimension(insets.left + insets.right + wMax + w2,
//                    insets.top + insets.bottom + y + (rows - 1)*vgap);
            return new Dimension(insets.left + insets.right,
                                 insets.top + insets.bottom + y + (rows - 1) * vgap);
        }
    }


    public Dimension minimumLayoutSize(Container parent) {
        Dimension size = preferredLayoutSize(parent);
        size.width = 0;
        return size;
    }


    public void layoutContainer(Container parent) {
        synchronized (parent.getTreeLock()) {
            Insets insets = parent.getInsets();
            int ncomponents = parent.getComponentCount();

            if (ncomponents == 0) {
                return;
            }
            rows = (ncomponents + cols - 1) / cols;

            // calculate the maximum width available
            int width = Integer.MIN_VALUE;
            for (int r = 0; r < rows; r++) {
                int i = r * cols;
                if (i < ncomponents) {
                    Dimension prefSize = parent.getComponent(i).getPreferredSize();
                    if (prefSize.width > width) {
                        width = prefSize.width;
                    }
                }
            }
            width = Math.max(parent.getWidth() - insets.left - insets.right - width,
                             minimumWidth);
            int[] w = {parent.getWidth() - insets.left - insets.right - width,
                      width - hgap};

            for (int r = 0, y = insets.top; r < rows; r++) {
                int maxH = Integer.MIN_VALUE;
                for (int c = 0; c < cols; c++) {
                    int i = r * cols + c;
                    if (i < ncomponents) {
                        int h = parent.getComponent(i).getPreferredSize().height;
                        if (h > maxH) {
                            maxH = h;
                        }
                    }
                }

                for (int c = 0, x = insets.left; c < cols; c++) {
                    int i = r * cols + c;
                    if (i < ncomponents) {
                        Component component = parent.getComponent(i);
                        int h = component.getPreferredSize().height;
                        parent.getComponent(i).setBounds(x, y + (maxH - h) / 2, w[c], h);
                        x += w[c] + hgap;
                    }
                }

                y += maxH + vgap;
            }
        }
    }

}
