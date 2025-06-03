package edu.ucsc.neurobiology.vision.plot;

import java.beans.*;
import java.lang.reflect.*;
import java.util.*;

import java.awt.*;
import javax.swing.*;

import edu.ucsc.neurobiology.vision.parameters.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class StyleCustomizer
    extends LazyPanel {

    ParametersTable symbolTable = new ParametersTable();


    public StyleCustomizer(
        final DrawingControl drawingControl, final Object style) throws
        NoSuchMethodException, IllegalAccessException, InvocationTargetException {

        this.setLayout(new GridLayout(0, 1));
        this.setOpaque(false);

        HashMap<String, Method> map = new HashMap<String, Method> ();
        final Method[] methods = style.getClass().getMethods();
        String[] methodNames = new String[methods.length];
        for (int i = 0; i < methods.length; i++) {
            methodNames[i] = methods[i].getName();
            map.put(methodNames[i], methods[i]);
        }
        Arrays.sort(methodNames);

        for (int n = 0; n < methodNames.length; n++) {
            final Method getMethod = map.get(methodNames[n]);
            String getter = getMethod.getName();
            String propertyName = getter.substring(3);

            if (getter.startsWith("get") && map.containsKey("set" + propertyName)) {
                try {
                    final Method setMethod = style.getClass().getMethod(
                        "set" + propertyName, new Class[] {getMethod.getReturnType()});

                    String returnType = getMethod.getReturnType().getName();
                    final Object value = getMethod.invoke(style);

                    final Enum[] enums = (Enum[]) getMethod.getReturnType().
                                         getEnumConstants();
                    if (enums != null) {
                        // EnumeratorParameter
                        final EnumeratorParameter p = new EnumeratorParameter(
                            propertyName, null, null);
                        for (Enum e : enums) {
                            p.addChoice(e.ordinal(), e.toString().replaceAll("_", " "));
                        }

                        symbolTable.addParameter(p, new PropertyChangeListener() {
                            public void propertyChange(PropertyChangeEvent event) {
                                try {
                                    setMethod.invoke(
                                        style, new Object[] {enums[ (int) p.getValue()]});
                                    drawingControl.updateNeeded(style);
                                } catch (Exception e) {
                                    e.printStackTrace();
                                }
                            }
                        });

                        p.setValue( ( (Enum) value).ordinal());
                    } else if (returnType.equals("int")) {
                        // IntegerParameter
                        final IntegerParameter p = new IntegerParameter(
                            propertyName, null, null, ( (Integer) value).intValue(),
                            Integer.MIN_VALUE,
                            Integer.MAX_VALUE);
                        symbolTable.addParameter(p, new PropertyChangeListener() {
                            public void propertyChange(PropertyChangeEvent event) {
                                try {
                                    setMethod.invoke(
                                        style, new Object[] {new Integer(p.getValue())});
                                    drawingControl.updateNeeded(style);
                                } catch (Exception e) {
                                    e.printStackTrace();
                                }
                            }
                        });
                    } else if (returnType.equals("float")) {
                        // IntegerParameter
                        final DoubleParameter p = new DoubleParameter(
                            propertyName, null, null, ( (Float) value).floatValue());
                        symbolTable.addParameter(p, new PropertyChangeListener() {
                            public void propertyChange(PropertyChangeEvent event) {
                                try {
                                    setMethod.invoke(
                                        style, new Object[] {new Float(p.getValue())});
                                    drawingControl.updateNeeded(style);
                                } catch (Exception e) {
                                    e.printStackTrace();
                                }
                            }
                        });
                    } else if (returnType.equals("boolean")) {
                        final BooleanParameter p = new BooleanParameter(propertyName, null, null,
                            ( (Boolean) value).booleanValue());
                        symbolTable.addParameter(p, new PropertyChangeListener() {
                            public void propertyChange(PropertyChangeEvent event) {
                                try {
                                    setMethod.invoke(
                                        style, new Object[] {new Boolean(p.getValue())});
                                    drawingControl.updateNeeded(style);
                                } catch (Exception e) {
                                    e.printStackTrace();
                                }
                            }
                        });
                    } else if (returnType.equals("java.awt.Color")) {
                        final ColorParameter p = new ColorParameter(propertyName, null, null,
                            (Color) value);
                        symbolTable.addParameter(p, new PropertyChangeListener() {
                            public void propertyChange(PropertyChangeEvent event) {
                                try {
                                    setMethod.invoke(style, new Object[] {p.getColor()});
                                    drawingControl.updateNeeded(style);
                                } catch (Exception e) {
                                    e.printStackTrace();
                                }
                            }
                        });
                    } else if (returnType.equals("java.lang.String")) {
                        final StringParameter p = new StringParameter(propertyName, null, null,
                            (String) value);
                        symbolTable.addParameter(p, new PropertyChangeListener() {
                            public void propertyChange(PropertyChangeEvent event) {
                                try {
                                    setMethod.invoke(style, new Object[] {p.getValue()});
                                    drawingControl.updateNeeded(style);
                                } catch (Exception e) {
                                    e.printStackTrace();
                                }
                            }
                        });
                    }
                } catch (NoSuchMethodException e) {
                    continue;
                }
            }
        }
    }


    public void lazyConstructor() {
        JPanel symbolPanel = new JPanel(new BorderLayout());
//        symbolPanel.setBorder(BorderFactory.createTitledBorder(
//            BorderFactory.createEtchedBorder(EtchedBorder.RAISED), "  Properties "));
        symbolPanel.add(new JScrollPane(symbolTable));

        add(symbolPanel);
    }

}
