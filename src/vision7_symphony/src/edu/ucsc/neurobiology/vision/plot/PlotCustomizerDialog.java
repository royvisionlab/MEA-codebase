package edu.ucsc.neurobiology.vision.plot;

import java.beans.*;
import java.util.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

import edu.ucsc.neurobiology.vision.parameters.*;
import static edu.ucsc.neurobiology.vision.plot.PlotPanel.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class PlotCustomizerDialog
    extends JDialog {

    private JButton okButton, defaultButton;

    private RangeParameter xRangeParam, yRangeParam;
    private BooleanParameter legendVisibleParam;
    private EnumeratorParameter legendLocationParam;

    private PlotPanel panel;
    private HashMap panelsMap = new HashMap();
    private Component oldComponent;
    private JList list;
//    private JPopupMenu menu = new JPopupMenu();
//    private int index;
    private JPanel emptyCustomizer = new JPanel();


    public PlotCustomizerDialog(DrawingControl drawingControl, final PlotPanel panel) {
        this.setTitle("Plot Customizer");
        this.panel = panel;

        final DefaultListModel model = new DefaultListModel();
        list = new JList(model);
        list.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        list.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));

        ParametersTable table = new ParametersTable();

        // axes controls, range
        PropertyChangeListener axesListener = new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent e) {
                updateAxes();
            }
        };

        double[] range = panel.getRange();
        double[] limits = panel.getLimits();
        xRangeParam = new RangeParameter("X Axis Range", null, null, range[0], range[1], limits[0], limits[1]);
        yRangeParam = new RangeParameter("Y Axis Range", null, null, range[2], range[3], limits[2], limits[3]);
        table.addParameter(xRangeParam, axesListener);
        table.addParameter(yRangeParam, axesListener);

        // the font parameters
        final String[] availFonts = GraphicsEnvironment.getLocalGraphicsEnvironment().
                                    getAvailableFontFamilyNames();
        Arrays.sort(availFonts);
        final EnumeratorParameter tickFontParam = new EnumeratorParameter("Tick Font", null, null);
        for (int i = 0; i < availFonts.length; i++) {
            tickFontParam.addChoice(i, availFonts[i]);
        }
        Font font = panel.getLabelFont();
        int index = Arrays.binarySearch(availFonts, font.getName());
        if (index >= 0) {
            tickFontParam.setValue(index);
        } else {
            throw new Error("wrong font");
        }
        final IntegerParameter fontSizeParam = new IntegerParameter("Font Size", null, null,
            font.getSize(), 1, 1000);
        final EnumeratorParameter fontStyle = new EnumeratorParameter("Font Style", null, null);
        fontStyle.addChoice(Font.PLAIN, "Plain");
        fontStyle.addChoice(Font.ITALIC, "Italic");
        fontStyle.addChoice(Font.BOLD, "Bold");
        fontStyle.setValue(font.getStyle());

        PropertyChangeListener fontListener = new PropertyChangeListener() {
            public void propertyChange(PropertyChangeEvent e) {
                panel.setLabelFont(new Font(
                    availFonts[ (int) tickFontParam.getValue()],
                    (int) fontStyle.getValue(),
                    fontSizeParam.getValue()));
            }
        };
        table.addParameter(tickFontParam, fontListener);
        table.addParameter(fontSizeParam, fontListener);
        table.addParameter(fontStyle, fontListener);

        legendVisibleParam = new BooleanParameter(
            "Legend Visible", null, null, panel.isLegendVisible());
        table.addParameter(legendVisibleParam, axesListener);

        legendLocationParam = new EnumeratorParameter("Legend Location", null, null);
        legendLocationParam.addChoice(Location.UPRIGHT.ordinal(), "UpRight");
        legendLocationParam.addChoice(Location.DOWNRIGHT.ordinal(), "DownRight");
        legendLocationParam.addChoice(Location.UPLEFT.ordinal(), "UpLeft");
        legendLocationParam.addChoice(Location.DOWNLEFT.ordinal(), "DownLeft");
        legendLocationParam.setValue(panel.getLegendLocation().ordinal());
        table.addParameter(legendLocationParam, axesListener);

        panelsMap.put("General", new JScrollPane(table));
        model.addElement("General");

        ArrayList uniqueStyles = new ArrayList();
        for (int i = 0; i < panel.plotData.size(); i++) {
            Object style = panel.styles.get(i);
            if (!uniqueStyles.contains(style)) {
                uniqueStyles.add(style);
            }
        }

        // Axes
        try {
            panelsMap.put("X Axis", new StyleCustomizer(drawingControl, panel.getXAxis()));
            model.addElement("X Axis");

            panelsMap.put("Y Axis", new StyleCustomizer(drawingControl, panel.getYAxis()));
            model.addElement("Y Axis");
        } catch (Exception e) {
            e.printStackTrace();
        }

        try {
            panelsMap.put("Padding", new StyleCustomizer(drawingControl, panel.getBorder()));
            model.addElement("Padding");
        } catch (Exception e) {
            e.printStackTrace();
        }

        for (int i = 0; i < uniqueStyles.size(); i++) {
            PlotStyle style = (PlotStyle) uniqueStyles.get(i);

            JPanel customizer = null;
            try {
                customizer = new StyleCustomizer(drawingControl, style);
            } catch (Exception e) {
                e.printStackTrace();
            }

            if (customizer == null) {
                customizer = emptyCustomizer;
            }

            String desc = style.getDescription().trim();
            if (desc == null || panelsMap.containsKey(desc)) {
                throw new InternalError(
                    "Repeating " + style.getClass().getName() + " description \"" + desc +
                    "\"");
            } else {
                panelsMap.put(desc, customizer);
                model.addElement(desc);
            }
        }

        list.addListSelectionListener(new ListSelectionListener() {
            public void valueChanged(ListSelectionEvent e) {
                Object selected = list.getSelectedValue();
                if (selected == null) {
                    return;
                }
                Component c = (Component) panelsMap.get(selected);

                if (oldComponent != null) {
                    getContentPane().remove(oldComponent);
                }
                add(c, BorderLayout.CENTER);
                oldComponent = c;

                validate();
                repaint();
            }
        });

        JScrollPane scroll = new JScrollPane(list);
        scroll.setPreferredSize(new Dimension(100, 0));
        add(scroll, BorderLayout.WEST);
        list.setSelectedIndex(0);
        add(createButtonsPanel(), BorderLayout.SOUTH);

        JFrame f = (JFrame) SwingUtilities.getAncestorOfClass(JFrame.class, panel);
        int w = 400;
        int h = 400;
        Point p = f.getLocation();
        Dimension d = f.getSize();
        this.setBounds(p.x + d.width / 2 - w / 2, p.y + d.height / 2 - h / 2, w, h);

        this.setVisible(true);
        this.setAlwaysOnTop(true);
    }


    private void updateAxes() {
        panel.setRange(xRangeParam.getMin(), xRangeParam.getMax(),
                       yRangeParam.getMin(), yRangeParam.getMax());
        panel.setLegendVisible(legendVisibleParam.getValue());
        panel.setLegendLocation(Location.values()[ (int) legendLocationParam.getValue()]);
        panel.replotAllData();
    }


    private JPanel createButtonsPanel() {
        JPanel buttonPanel = new JPanel(new FlowLayout(FlowLayout.RIGHT));

        okButton = new JButton("Ok");
        okButton.setMnemonic(KeyEvent.VK_O);
        okButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                updateAxes();
                PlotCustomizerDialog.this.dispose();
            }
        });
        buttonPanel.add(okButton);
        /*
                cancelButton = new JButton("Cancel");
                cancelButton.setMnemonic(KeyEvent.VK_C);
                cancelButton.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        PlotCustomizerDialog.this.dispose();
                    }
                });
                buttonPanel.add(cancelButton);
         */
        defaultButton = new JButton("Default");
        defaultButton.setMnemonic(KeyEvent.VK_D);
        buttonPanel.add(defaultButton);

        return buttonPanel;
    }

}
