package edu.ucsc.neurobiology.vision.parameters;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;


/**
 * This class implements the <code>Parameter</code> interface to represent a
 * range in double numbers represented graphically by two text fields and a set
 * of navigation buttons. It also provides getter and setter functions for the
 * minimum and maximum range values stored in the parameter. Every time the
 * minimum or maximum changes a <code>PropetryChangeEvent</code> is dispatched
 * to all listeners.
 *
 * @author Matthew Grivich, The Salk Institute
 */
public class RangeParameter
extends AbstractParameter {


    private JPanel panel;
    public RangeSpinner spin;
    private double _x1, _x2, _min, _max;




    /**
     * Creates an instance of RangeParameter which will display by default
     * the range given by x1 (lower limit) and x2 (upper limit). x2 <b>must</b>
     * be greater than x1.
     *
     * @param name the internal name of the parameter
     * @param screenName the screen name associated with the parameter
     * @param toolTip the explanatory tooltip associated with the parameter
     * @param x1 the default lower limit of the range
     * @param x2 the default upper limit of the range
     */
    public RangeParameter(String name, String screenName, String toolTip, double x1,
            double x2, double min, double max) {
        super(name, screenName, toolTip);

        _x1 = x1;
        _x2 = x2;
        _min = min;
        _max = max;

        initialize();
    }




    public void initialize() {
        super.initialize();
        
        spin = new RangeSpinner(_x1, _x2);
        //Tells any listener for the range parameter when it changes.
        

        spin.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent e) {
                JSpinner spinner = (JSpinner)e.getSource();
                _x1 = ((double[]) spinner.getValue())[0];
                _x2 = ((double[]) spinner.getValue())[1];
                fireParameterChange();
                //System.out.println("change");
            }
        });

    }


    public JComponent getValueComponent() {
        return spin;
    }


    /**
     * Returns the lower limit of the range.
     */
//	public double getMin() {
//	try {
//	return Double.parseDouble(x1Field.getText());
//	} catch (NumberFormatException e) {
//	return Double.NaN;
//	}
//	}


    /**
     */
//	public void setMin(double min) {
//	x1Field.setText("" + min);
//	}


    /**
     * Returns the upper limit of the range.
     */
//	public double getMax() {
//	try {
//	return Double.parseDouble(x2Field.getText());
//	} catch (NumberFormatException e) {
//	return Double.NaN;
//	}
//	}


    /**
     */
//	public void setMax(double max) {
//	x2Field.setText("" + max);
//	}

    public double getMin() {
        return _x1;
    }

    public double getMax() {
        return _x2;
    }


    public String valueAsString() {
        return null;
    }


    public Object valueAsObject() {
        return null;
    }





    public class RangeSpinner extends JSpinner {

        private JFormattedTextField x1Field, x2Field;
        public double[] range = {0,0}; 


        public RangeSpinner(double x1, double x2) {
            super();
            range = new double[2];
            range[0] = x1;
            range[1] = x2;
            setModel(new SpinnerRangeModel());
            setEditor(new Editor(this));
            fireStateChanged();

        }

        protected JComponent createEditor(SpinnerModel model) {
            return  (new Editor(this));
        }

        public class SpinnerRangeModel extends AbstractSpinnerModel {

            public Object getNextValue() {
                double increment = range[1] - range[0];
                if(range[1] + increment >= _max) {
                    return range;
                }
                range[0]+=increment;
                range[1]+=increment;
                return range;
            }

            public Object getPreviousValue() {
                double increment = range[1] - range[0];
                if(range[0] - increment <= _min) {
                    return range;
                }
                range[0]-=increment;
                range[1]-=increment;
                return range;
            }

            public Object getValue() {
                return range;
            }

            public void setValue(Object value) {
                range[0] = ((double[]) value)[0];
                range[1] = ((double[]) value)[1];
                fireStateChanged();
            }

        }

        public class Editor extends JPanel implements ChangeListener {



            Editor(RangeSpinner spinner) {
                this.setLayout(new GridLayout(1, 2));
                
                x1Field = new JFormattedTextField();

                x1Field.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent event) {
                            double[] r = (double[]) getValue();
                            double newVal = Double.parseDouble(x1Field.getText());

                            if(newVal >= _min && (newVal + r[1]-r[0]) <= _max) {
                                //if(newVal > r[1]) newVal = r[1];
                                setValue(new double[]{newVal,newVal + r[1] - r[0]});
                                
                            //prevent change
                            } else {

                                setValue(new double[]{r[0], r[1]});
                            }

                    }
                });
                
                x1Field.addFocusListener(new FocusListener() {
                    public void focusLost(FocusEvent e) {
                    
                            double[] r = (double[]) getValue();
                            double newVal = Double.parseDouble(x1Field.getText());

                            if(newVal >= _min && (newVal + r[1]-r[0]) <= _max) {
                                //if(newVal > r[1]) newVal = r[1];
                                setValue(new double[]{newVal,newVal + r[1] - r[0]});
                                
                            //prevent change
                            } else {

                                setValue(new double[]{r[0], r[1]});
                            }

                    }
                    
                    public void focusGained(FocusEvent e) {}
                    
                });


                x2Field = new JFormattedTextField();			
                x2Field.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent event) {
                            double[] r = (double[]) getValue();
                            double newVal = Double.parseDouble(x2Field.getText());
                            if(newVal <= _max) {
                                if(newVal < r[0]) newVal = r[0];
                                setValue(new double[] {r[0], newVal});
                                
                            //prevent change
                            }else {
                                setValue(new double[]{r[0], r[1]});
                            }
    
                    }
                });
                
                x2Field.addFocusListener(new FocusListener() {
                    public void focusLost(FocusEvent e) {
                    
                        double[] r = (double[]) getValue();
                        double newVal = Double.parseDouble(x2Field.getText());
                        if(newVal <= _max) {
                            if(newVal < r[0]) newVal = r[0];
                            setValue(new double[] {r[0], newVal});
                            
                        //prevent change
                        }else {
                            setValue(new double[]{r[0], r[1]});
                        }

                    }
                    
                    public void focusGained(FocusEvent e) {}
                    
                });



                this.add(x1Field);
                this.add(x2Field);
                // Add the listener
                spinner.addChangeListener(this);

                // Set the preferred size
//				setPreferredSize(new Dimension(preferredWidth, preferredHeight));

                // Display the current color
                //x1Field.setV
                //setColor((String)spinner.getValue());

            }

            // Handles changes to the value stored in the model
            public void stateChanged(ChangeEvent evt) {
                JSpinner spinner = (JSpinner)evt.getSource();




                x1Field.setText("" + ((double[]) spinner.getValue())[0]);
                x2Field.setText("" + ((double[]) spinner.getValue())[1]);


            }


        }



    }
    


    public static void main(String[] args) {
        JFrame jFrame = new JFrame();
        jFrame.setSize(100,100);

        RangeParameter rp = new RangeParameter("t", "n", "s", 0, 1, 0, 1000);


        jFrame.add(rp.spin);
        jFrame.setVisible(true);


    }
}
