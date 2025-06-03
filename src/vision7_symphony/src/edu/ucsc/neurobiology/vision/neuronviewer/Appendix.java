package edu.ucsc.neurobiology.vision.neuronviewer;


public class Appendix {

    /*
            public void showResponseCurve(String text, int ...id) throws IOException {
                int tPeriod = 0;

                PlotPanel p = new PlotPanel();
                ScatterPlotStyle style = new ScatterPlotStyle(
                    SymbolType.FILLED_SQUARE, 4, black, true, black, 1);

                double[] x = new double[nSpatialPeriods];
                for (int i = 0; i < x.length; i++) {
                    x[i] = getFrequency(i);
                }
                double[] tunningAv = new double[nSpatialPeriods];

                for (int i = 0; i < id.length; i++) {
                    this.setCurrentNeuron(id[i]);
                    DoubleHistogram[] hh = getOnePeriodHistograms(1);
                    double[] tunning = new double[nSpatialPeriods];

                    for (int sPeriod = 0; sPeriod < nSpatialPeriods; sPeriod++) {
                        double max = Double.NEGATIVE_INFINITY;
                        for (int sPhase = 0; sPhase < nPhases; sPhase++) {
     int runID = getRunID(spatialPeriods[sPeriod], temporalPeriods[tPeriod],
                                                 phases[sPhase]);

                            int n = hh[runID].getMaxValueBin();
                            double v = 0;
                            for (int k = n - 2; k <= n + 2; k++) {
                                if (k > 0 && k < hh[runID].getBinCount()) {
                                    v += hh[runID].getBin(k);
                                }
                            }

                            if (v > max) {
                                max = v;
                            }
                        }
                        tunning[sPeriod] = max;
                    }

                    MathUtil.divide(tunning, MathUtil.sum(tunning));
                    MathUtil.add(tunningAv, tunning);
                    p.addData(new ScatterPlot(x, tunning, null), style);
                }

                MathUtil.divide(tunningAv, id.length);

                PlotUtil.showData(text, p);
                p.setAxesType(AxisType.LOG10, AxisType.LINEAR);
     p.addData(new ScatterPlot(x, tunningAv, null), "DISK 8 red, SOLID 3 red");
                p.autoscale();
            }


            public double[] getResponseCurve() throws IOException {
                int tPeriod = 0;
                DoubleHistogram[] hh = getOnePeriodHistograms(1);
                double[] r = new double[nSpatialPeriods];
                for (int sPeriod = 0; sPeriod < nSpatialPeriods; sPeriod++) {
                    double max = Double.NEGATIVE_INFINITY;
                    for (int sPhase = 0; sPhase < nPhases; sPhase++) {
     int runID = getRunID(spatialPeriods[sPeriod], temporalPeriods[tPeriod],
                                             phases[sPhase]);

                        int n = hh[runID].getMaxValueBin();
                        double v = 0;
                        for (int k = n - 2; k <= n + 2; k++) {
                            if (k > 0 && k < hh[runID].getBinCount()) {
                                v += hh[runID].getBin(k);
                            }
                        }

                        if (v > max) {
                            max = v;
                        }
                    }
                    r[sPeriod] = max;
                }

                return r;
            }

     */

    /*
        CalculationAction makeRefractoryPeriodHistogramAction =
            new CalculationAction("Make Refractory Period Histogram",
                                  CalculationAction.CLASS_ACTION) {

            public void doAction(IntegerList list, final TreePath classTreePath) {
                DoubleHistogram h = new DoubleHistogram("", 0, 200, 1);
                final HashMap<Integer, Double> valueMap = new HashMap();

                try {
                    for (int n = 0; n < list.size(); n++) {
                        int[] t = neuronFile.getSpikeTimes(list.get(n));
     DoubleHistogram isiHist = AutocorrelationCalculator.calculate(t, 200,
                            0.05);
                        double[] acf = isiHist.toArray();

                        double p = 0;
                        double threshold = 0.0004 * t.length;

                        for (int i = 10; i < acf.length; i++) {
                            p += acf[i];
                            if (p > threshold) {
                                h.fill(i, 1);
                                valueMap.put(new Integer(list.get(n)), new Double(i));
                                break;
                            }
                        }

                    }
                } catch (IOException e) {
                    VisionUtil.reportException(e);
                }

                PlotPanel p = new PlotPanel();
                p.addSelectionAction(new SelectionAction("New Class") {
     public void selectionPerformed(JComponent source, Selection selection) {
                        SelectionRegion r = selection.getSelection();
                        IntegerList idList = new IntegerList();
                        Rectangle2D b = r.getBounds2D();
                        double x1 = b.getX();
                        double x2 = x1 + b.getWidth();

                        for (Integer id : valueMap.keySet()) {
                            Double ox = valueMap.get(id);
                            if (ox != null) {
                                double x = ox.doubleValue();
                                if (x > x1 && x < x2) {
                                    idList.add(id.intValue());
                                }
                            }
                        }

                        String className = "nc" + newClassIndex++;
                        NeuronViewer.this.classificationHappened(idList.toArray(),
                            classTreePath.pathByAddingChild(new DefaultMutableTreeNode(
                                className, true)));
                    }
                });
                p.addData(h, new HistogramStyle("Refractory Period"));
                p.autoscale();
                showData("", p);
            }
        };
     */


    /*
        public static void profileCreate() {
            Thread t = new Thread() {
                public void run() {
                    double t1 = System.currentTimeMillis();
                    for (int i = 0; i < 10; i++) {
                        try {
                            NeuronViewer s = new NeuronViewer("c:\\data003", false);
                        } catch (IOException e) {
                            e.printStackTrace();
                        }
                    }
                    double t2 = System.currentTimeMillis();
                    System.out.println( (t2 - t1) / 1000.0);
                }
            };
            t.start();
        }
     */

    /*
                if (action.equals("Export Time Courses (Right)")) {
                    if (rightSelection instanceof Integer) {
                        int neuronID = ( (Integer) rightSelection).intValue();
                        System.out.println("Export " + neuronID);

                        double[] redTimeFilter = paramsFile.geArrayCell(
                            neuronID, "RedTimeCourse");
                        double[] greenTimeFilter = paramsFile.geArrayCell(
                            neuronID, "GreenTimeCourse");
                        double[] blueTimeFilter = paramsFile.geArrayCell(
                            neuronID, "BlueTimeCourse");

                        try {
                            PrintWriter w = new PrintWriter(new FileWriter("tc" +
                                neuronID + ".txt"));
                            for (int i = 0; i < staDepth; i++) {
                                w.println( (i - staDepth + 1) * refreshInterval
                                          + "\t" + redTimeFilter[i]
                                          + "\t" + greenTimeFilter[i]
                                          + "\t" + blueTimeFilter[i]);
                            }
                            w.close();
                        } catch (IOException e) {
                            VisionUtil.reportException(e);
                        }
                    } else {
                        JOptionPane.showMessageDialog(
     mainFrame, "Wrong Selection. You should select a neuron !", "Error",
                            JOptionPane.ERROR_MESSAGE);
                    }
     */

    /*
        CalculationAction showF1F2intersectionPlot =
     new CalculationAction("Show F1 F2 Intersection", CalculationAction.CLASS_ACTION) {

            public void doAction(final IntegerList list, final TreePath classTreePath) {
                //            double f = 1 / spatialPeriods[sPeriod] / mmPerPixel;
//            double[] freqs = paramsFile.getDoubleArrayCell(id[0], "reversingFrequencies");

                DoubleHistogram n1 = new DoubleHistogram("", 0, 10, 0.15);

                for (int k = 0; k < list.size(); k++) {
                    int id = list.get(k);

                    double[] f1 = paramsFile.geArrayCell(id, "T1reversingF1");
                    double[] f2 = paramsFile.geArrayCell(id, "T1reversingF2");

                    for (int i = f1.length - 1; i > 0; i--) {
                        Point2D.Double p = new Point2D.Double();
                        MathUtil.getLineIntersection(
                            i, f1[i], i - 1, f1[i - 1],
                            i, f2[i], i - 1, f2[i - 1], p, true);
                        if (!Double.isNaN(p.x)) {
                            n1.fill(p.x, 1);
                            break;
                        }
                    }
                }

                PlotUtil.showData(classTreePath.toString(), n1, new HistogramStyle());
            }
        };
     */
    /*
        CalculationAction showAxonSpeedPlot =
     new CalculationAction("Show Axon Speed Plot", CalculationAction.CLASS_ACTION) {

            public void doAction(final IntegerList list, final TreePath classTreePath) {
                             try {
                    EITest t = new EITest(filePathRoot + ".ei");
                    DoubleHistogram h = new DoubleHistogram("", 0, 5, 0.05);
                    for (int i = 0; i < list.size(); i++) {
                        h.fill(t.getSpeed(list.get(i)), 1);
                    }
                    PlotUtil.showData(classTreePath.toString(), h, new HistogramStyle());
                             } catch (IOException ex) {
                    ex.printStackTrace();
                             }
            }
        };
     */
    /*
        CalculationAction savePlots =
            new CalculationAction("Save Plots", CalculationAction.CLASS_ACTION) {

            public void doAction(final IntegerList _list, final TreePath classTreePath) {
                double d = 4.5;

                if (rightSelection instanceof String) {
                    String path = InteractiveTree.pathToString(classTreePath);
                    path = path.substring(4);
                    path = path.replace("/", " ");

                    PlotPanel mosaic = (PlotPanel) rightSplitPane.getLeftComponent();
                    try {
                        mosaic.setLegendVisible(false);
                        mosaic.saveAsEPS(new File(experimentName + " " + path + ", " +
                                                  "mosaic.eps"), d, d / 2);
                        mosaic.setLegendVisible(true);
                    } catch (IOException ex) {
                        ex.printStackTrace();
                    }

                    JPanel p = (JPanel) rightSplitPane.getRightComponent();
                    for (int i = 0; i < p.getComponentCount(); i++) {
                        if (p.getComponent(i) instanceof PlotPanel) {
                            PlotPanel panel = (PlotPanel) p.getComponent(i);
                            try {
                                panel.setLegendVisible(false);
                                String name = panel.getName();
                                if (name == null || name.length() == 0) {
                                    name = "" + i;
                                }
                                panel.saveAsEPS(new File(
     experimentName + " " + path + ", " + name + ".eps"), d,
                                                d * 0.4);
                                panel.setLegendVisible(true);
                            } catch (IOException ex) {
                                ex.printStackTrace();
                            }
                        }
                    }
                }
            }
        };
     */
    /*
        CalculationAction removeDuplicates =
            new CalculationAction("Remove Duplicates", CalculationAction.CLASS_ACTION) {

            public void doAction(final IntegerList list, final TreePath classTreePath) {
                SwingWorker worker = new SwingWorker() {
                    public Object construct() {
                        DoubleParameter minRedChiSq = new DoubleParameter(
                            "Minimum Reduced Chi-Squared", null, null, 4.0);
                        DoubleParameter maxSeparation = new DoubleParameter(
                            "Max Seperation (microns)", null, null, 500.0);
                        ParametersTable table = new ParametersTable();
                        table.addParameter(minRedChiSq);
                        table.addParameter(maxSeparation);

                        ParametersDialog dialog = new ParametersDialog(
                            mainFrame, "Choose Parameters", table);

                        dialog.setVisible(true);
                        if (dialog.isCancelSelected()) {
                            return null;
                        }

                        try {
                            return removeDuplicates(list, maxSeparation.getValue(),
                                                    minRedChiSq.getValue());
                        } catch (IOException e) {
                            VisionUtil.reportException(e);
                            return null;
                        }
                    }


                    public void finished() {
                        int[] neurons = (int[]) getValue();
                        if (neurons != null && neurons.length != 0) {
                            NeuronViewer.this.classificationHappened(
                                neurons, classTreePath.pathByAddingChild(new
                                DefaultMutableTreeNode("Duplicates", true)));
                        }
                    }
                };
                worker.start();
            }
        };
     */
    /*
        CalculationAction removeDuplicatesByEI =
     new CalculationAction("Remove Duplicates By EI", CalculationAction.CLASS_ACTION) {
            public void doAction(IntegerList list, final TreePath classTreePath) {
                ArrayList neuronsInClass = new ArrayList();
                ArrayList duplicateNeurons = new ArrayList();
                ArrayList imgList = new ArrayList();
                for (int i = 0; i < list.size(); i++) {
                    neuronsInClass.add(new Integer(list.get(i)));
                    imgList.add(imgFile.getImage(list.get(i)));
                }
                System.out.println("Compare");
                DoubleParameter minRedChiSq = new DoubleParameter(
                    "Minimum Reduced Chi-Squared", null, null, 4.0);
                DoubleParameter ChiSq = new DoubleParameter(
                    "Chi-Squared", null, null, 10);
                ParametersTable table = new ParametersTable();
                table.addParameter(minRedChiSq);
                table.addParameter(ChiSq);
                ParametersDialog dialog = new ParametersDialog(
                    frame, "Choose Parameters", table);
                dialog.setVisible(true);
                if (dialog.getResult() == ParametersDialog.CANCEL_OPTION) {
                    return;
                }
                double thresholdChi2 = minRedChiSq.getValue();
                DoubleHistogram h = new DoubleHistogram("", 0, ChiSq.getValue(), 0.01);
                for (int n1 = 0; n1 < neuronsInClass.size(); n1++) {
                    int neuron1 = ( (Integer) neuronsInClass.get(n1)).intValue();
                    System.out.println("Comparing Neuron " + neuron1);
                    float[][][] image1 = (float[][][]) imgList.get(n1);
         double nSpikes1 = ( (Double) (paramsFile.getCell(neuron1, "nSpikes")).
                                       getEmbeddedObject()).doubleValue();
                    for (int n2 = n1 + 1; n2 < neuronsInClass.size(); n2++) {
                        int neuron2 = ( (Integer) neuronsInClass.get(n2)).intValue();
                        float[][][] image2 = (float[][][]) imgList.get(n2);
         double nSpikes2 = ( (Double) (paramsFile.getCell(neuron2, "nSpikes")).
                                           getEmbeddedObject()).doubleValue();
                        double redChi2 = 0;
                        int n = 0;
                        for (int e = 1; e < imgFile.nElectrodes; e++) {
                            int i1 = minIndex(image1[0][e]);
                            int i2 = minIndex(image2[0][e]);
                            double a1 = -image1[0][e][i1];
                            double a2 = -image2[0][e][i2];
                            double err1 = image1[1][e][i1];
                            double err2 = image2[1][e][i2];
                            if (a1 > 4 * err1 || a2 > 4 * err2) {
                                for (int i = 0; i < imgFile.nSamples; i++) {
                                    double dx = (image1[0][e][i] - image2[0][e][i]);
                                    double e1 = image1[1][e][i];
                                    double e2 = image2[1][e][i];
                                    redChi2 += dx * dx / (e1 * e1 + e2 * e2);
                                    n++;
                                }
                            }
                        }
                        redChi2 /= n;
                        h.fill(redChi2, 1);
                        if (redChi2 < thresholdChi2) {
                            if (nSpikes1 < nSpikes2) {
                                duplicateNeurons.add(neuronsInClass.get(n1));
                                System.out.println(
                                    "Removing Neuron " + neuron1 + " Compared to " +
                                    neuron2 + " Reduced Chi Squared: " + redChi2);
                            } else {
                                duplicateNeurons.add(neuronsInClass.get(n2));
                                System.out.println(
                                    "Removing Neuron " + neuron2 + " Compared to " +
                                    neuron1 + " Reduced Chi Squared: " + redChi2);
                            }
                        }
                    }
                }
//            String className = "Duplicates";
//            int[] duplicates = new int[duplicateNeurons.size()];
//            for (int i = 0; i < duplicateNeurons.size(); i++) {
//                duplicates[i] = ( (Integer) duplicateNeurons.get(i)).intValue();
//            }
//            STAViewer.this.classificationHappened(
//                duplicates, classTreePath.pathByAddingChild(new
//                DefaultMutableTreeNode(className, true)));
                showData("EI Comparioson Reduced Chi Squared Histogram", h,
                                       new HistogramStyle());
            }
        };
     */

    /*
         CalculationAction spacingHistogram =
        new CalculationAction("Spacing Histogram", CalculationAction.CLASS_ACTION) {

        public void doAction(IntegerList list, final TreePath classTreePath) {
            ParametersTable table = configuration.showDialog(
                "BasicHistogram", "Spacing Histogram", mainFrame);
            if (table == null) {
                return;
            }

            double xMin = ( (DoubleParameter) table.getParameter("Min")).getValue();
            double xMax = ( (DoubleParameter) table.getParameter("Max")).getValue();
            double binWidth = ( (DoubleParameter) table.getParameter("BinWidth")).
                              getValue();
            String classPath = InteractiveTree.pathToString(classTreePath);
            DefaultMutableTreeNode currentFolder = (DefaultMutableTreeNode)
                classTreePath.getLastPathComponent();
            DoubleHistogram h = new DoubleHistogram("", xMin, xMax, binWidth);
            DoubleHistogram hError = new DoubleHistogram("", xMin, xMax, binWidth);
            int count = 0;
            for (int k = 0; k < currentFolder.getChildCount(); k++) {
//            if (currentFolder.getChildAt(i).getAllowsChildren()) {}

                classPath = InteractiveTree.pathToString(classTreePath.
                    pathByAddingChild(currentFolder.getChildAt(k)));

                int[] neuronsInClass = getNeuronsInClass(classPath, paramsFile);

                double[] x0 = new double[neuronsInClass.length];
                double[] y0 = new double[neuronsInClass.length];

                for (int i = 0; i < neuronsInClass.length; i++) {
                    x0[i] = ( (Double) (paramsFile.getCell(neuronsInClass[i], "x0")).
                             getEmbeddedObject()).doubleValue();
                    y0[i] = ( (Double) (paramsFile.getCell(neuronsInClass[i], "y0")).
                             getEmbeddedObject()).doubleValue();
//System.out.println(x0[i] + "   " + y0[i]);
                }

                for (int i = 0; i < neuronsInClass.length; i++) {
                    for (int j = i + 1; j < neuronsInClass.length; j++) {
                        h.fill(
                            sqrt( (x0[i] - x0[j]) * (x0[i] - x0[j]) +
                                      (y0[i] - y0[j]) * (y0[i] - y0[j])), 1.0);
                        count++;
                    }
                }
            }
            for (int i = 0; i < h.getBinCount(); i++) {
                double divisor = count; //2*PI*i*h.getBinInterval();
                //    System.out.println(i +  "   " + h.getBin(i) + "   " +divisor);
                hError.setBin(i, sqrt(h.getBin(i)));
                h.setBin(i, h.getBin(i) / divisor);
                hError.setBin(i, (hError.getBin(i) / divisor) * 2);
                //    System.out.println(h.getBin(i));

            }

            DoubleErrorHistogram h2 = new DoubleErrorHistogram(h, hError);

            DoubleHistogram monte = new DoubleHistogram("", xMin, xMax, binWidth);
            DoubleHistogram monteError = new DoubleHistogram("", xMin, xMax, binWidth);
            int length = 1000; //137;
            double[] x0 = new double[length];
            double[] y0 = new double[length];
            double hypotenuese = 36;
            double aspectRatio = 2;
            double shortSide = sqrt(hypotenuese * hypotenuese /
                                         (aspectRatio * aspectRatio + 1));
            double longSide = shortSide * aspectRatio;
            System.out.println("Short: " + shortSide + " Long: " + longSide);
            for (int i = 0; i < length; i++) {
                x0[i] = random() * longSide;
                y0[i] = random() * shortSide;
            }

            count = 0;
            for (int i = 0; i < length; i++) {
                for (int j = i + 1; j < length; j++) {
                    monte.fill(
                        sqrt( (x0[i] - x0[j]) * (x0[i] - x0[j]) +
                                  (y0[i] - y0[j]) * (y0[i] - y0[j])), 1.0);
                    count++;
                }
            }

            for (int i = 0; i < h.getBinCount(); i++) {
                double divisor = count; //2*PI*i*h.getBinInterval();
                monteError.setBin(i, sqrt(monte.getBin(i)));

//                monte.setBin(i, monte.getBin(i) / divisor);
//                monteError.setBin(i, (monteError.getBin(i) / divisor) * 2);

                monte.setBin(i, monte.getBin(i) / divisor);
                monteError.setBin(i, (monteError.getBin(i) / divisor) * 2);

            }

            DoubleErrorHistogram monte2 = new DoubleErrorHistogram(monte, monteError);
            PlotPanel p = new PlotPanel();

            p.setLabels("Distance", "N");
            p.addData(h2, new HistogramStyle("Spacing"));
            p.addData(monte2, new HistogramStyle(red, 1));
            p.autoscale();
            showData(classPath, p);
        }
         };
     */

    /*
        CalculationAction crossSpacingHistogram =
     new CalculationAction("Cross Spacing Histogram", CalculationAction.CLASS_ACTION) {

            public void doAction(IntegerList list, final TreePath classTreePath) {
                ParametersTable table = configuration.showDialog(
                    "BasicHistogram", "Spacing Histogram", mainFrame);
                if (table == null) {
                    return;
                }

                double xMin = ( (DoubleParameter) table.getParameter("Min")).getValue();
                double xMax = ( (DoubleParameter) table.getParameter("Max")).getValue();
                double binWidth = ( (DoubleParameter) table.getParameter("BinWidth")).
                                  getValue();
                String classPath = InteractiveTree.pathToString(classTreePath);

                int[] neuronsInClass = getNeuronsInClass(classPath, paramsFile);
                int[] neuronsInReferenceClass = getNeuronsInClass(referenceClass,
                    paramsFile);
                DoubleHistogram h = new DoubleHistogram("", xMin, xMax, binWidth);
                DoubleHistogram hError = new DoubleHistogram("", xMin, xMax, binWidth);

                double[] x0 = new double[neuronsInClass.length];
                double[] y0 = new double[neuronsInClass.length];
                double[] x1 = new double[neuronsInReferenceClass.length];
                double[] y1 = new double[neuronsInReferenceClass.length];

                for (int i = 0; i < neuronsInClass.length; i++) {
                    x0[i] = ( (Double) (paramsFile.getCell(neuronsInClass[i], "x0")).
                             getEmbeddedObject()).doubleValue();
                    y0[i] = ( (Double) (paramsFile.getCell(neuronsInClass[i], "y0")).
                             getEmbeddedObject()).doubleValue();
                }

                for (int i = 0; i < neuronsInReferenceClass.length; i++) {
     x1[i] = ( (Double) (paramsFile.getCell(neuronsInReferenceClass[i], "x0")).
                             getEmbeddedObject()).doubleValue();
     y1[i] = ( (Double) (paramsFile.getCell(neuronsInReferenceClass[i], "y0")).
                             getEmbeddedObject()).doubleValue();
                }

                int count = 0;
                for (int i = 0; i < neuronsInClass.length; i++) {
                    for (int j = 0; j < neuronsInReferenceClass.length; j++) {
                        h.fill(
                            sqrt( (x0[i] - x1[j]) * (x0[i] - x1[j]) +
                                      (y0[i] - y1[j]) * (y0[i] - y1[j])), 1.0);
                        count++;
                    }
                }

                for (int i = 0; i < h.getBinCount(); i++) {
                    double divisor = count; //2*PI*i*h.getBinInterval();
                    //    System.out.println(i +  "   " + h.getBin(i) + "   " +divisor);
                    hError.setBin(i, sqrt(h.getBin(i)));
                    h.setBin(i, h.getBin(i) / divisor);
                    hError.setBin(i, (hError.getBin(i) / divisor) * 2);
                    //    System.out.println(h.getBin(i));

                }

                DoubleErrorHistogram h2 = new DoubleErrorHistogram(h, hError);

                DoubleHistogram monte = new DoubleHistogram("", xMin, xMax, binWidth);
     DoubleHistogram monteError = new DoubleHistogram("", xMin, xMax, binWidth);
                int length = 5000; //137;
                x0 = new double[length];
                y0 = new double[length];
                double hypotenuese = 36;
                double aspectRatio = 2;
                double shortSide = sqrt(hypotenuese * hypotenuese /
                                             (aspectRatio * aspectRatio + 1));
                double longSide = shortSide * aspectRatio;
                System.out.println("Short: " + shortSide + " Long: " + longSide);
                for (int i = 0; i < length; i++) {
                    x0[i] = random() * longSide + 32 - (longSide / 2);
                    y0[i] = random() * shortSide + 16 - (shortSide / 2);
                }

                count = 0;
                for (int i = 0; i < length; i++) {
                    for (int j = 0; j < neuronsInReferenceClass.length; j++) {
                        monte.fill(
                            sqrt( (x0[i] - x1[j]) * (x0[i] - x1[j]) +
                                      (y0[i] - y1[j]) * (y0[i] - y1[j])), 1.0);
                        count++;
                    }
                }

                for (int i = 0; i < monte.getBinCount(); i++) {
                    double divisor = count; //2*PI*i*h.getBinInterval();
                    monteError.setBin(i, sqrt(monte.getBin(i)));
                    monte.setBin(i, monte.getBin(i) / divisor);
                    monteError.setBin(i, (monteError.getBin(i) / divisor) * 2);
                }

     DoubleErrorHistogram monte2 = new DoubleErrorHistogram(monte, monteError);
                PlotPanel p = new PlotPanel();

                p.setLabels("Distance", "N");
                p.addData(h2, new HistogramStyle());
                p.addData(monte2, new HistogramStyle(red, 1));
                p.autoscale();
                showData(classPath, p);
            }
        };
     */

    /*
        CalculationAction overlayedPlotAction =
            new CalculationAction("Plot Overlayed", CalculationAction.CLASS_ACTION) {

            public void doAction(IntegerList list, TreePath classTreePath) {
                System.out.println(list.size() + " neurons");
                ParametersTable table = configuration.showDialog(
                    "PlotOverlayed", "Plot Overlayed", mainFrame);
                if (table == null) {
                    return;
                }

                String expression = ( (StringParameter) table.getParameter(
                    "Expression")).getValue();
                String classPath = InteractiveTree.pathToString(classTreePath);
                HashMap<Integer, double[]> xMap = null;
                try {
     xMap = paramsFile.evaluate(expression, "classID==\"" + classPath + "\"");
                } catch (CannotEvaluateException ex) {
                    VisionUtil.reportException(ex);
                    return;
                }
                PlotPanel p = new PlotPanel();
                for (Integer id : xMap.keySet()) {
                    double[] x = xMap.get(id);
                    if (x != null) {
                        p.addData(new ScatterPlot(x), averageStyle);
                    }
                }
                p.setLabels(expression, "");
                p.autoscale();
                showData(classPath, p);
            }
        };
     */
}
