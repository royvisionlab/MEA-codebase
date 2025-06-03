package edu.ucsc.neurobiology.vision.neuronviewer;

import java.awt.Component;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

import javax.swing.JPanel;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreePath;

import edu.ucsc.neurobiology.vision.Vision;
import edu.ucsc.neurobiology.vision.gui.GraphicsIO;
import edu.ucsc.neurobiology.vision.math.CannotEvaluateException;
import edu.ucsc.neurobiology.vision.math.MeanVarianceCalculator;
import edu.ucsc.neurobiology.vision.util.IntegerList;
import edu.ucsc.neurobiology.vision.util.StringUtil;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class GenerateReportAction
    extends CalculationAction {

    class Expression {
        String name;
        String expression;
        int nDigits;

        public Expression(String name, String expression, int nDigits) {
            if (expression == null || expression.trim().length() == 0) {
                throw new Error("Connot use expression " + expression);
            }
            this.expression = expression;
            if (name == null || name.trim().length() == 0) {
                this.name = expression;
            } else {
                this.name = name;
            }
            this.nDigits = nDigits;
        }
    }


    ArrayList<Expression> expressions = new ArrayList<Expression>();

    String filePathRoot;
    String experimentName;
    String datasetName;


    public GenerateReportAction() {
        super("Generate Report", CalculationAction.CLASS_ACTION);

        expressions.add(new Expression("RF Diam.(um)", "2*116*((SigmaX*SigmaY)^0.5)", 0));
        expressions.add(new Expression("t1 (ms)", "-t1", 1));
        expressions.add(new Expression("t2 (ms)", "-t2", 1));
        expressions.add(new Expression("DOT", "dot2", 2));
        expressions.add(new Expression("Latency (ms)", "-rl", 1));
        expressions.add(new Expression("n-index", "extreme(T1reversingF2/T1reversingF1)",
                                       2));
    }


    public void initialize(NeuronViewer viewer) {
        super.initialize(viewer);

        filePathRoot = viewer.filePathRoot;
        experimentName = viewer.experimentName;
        datasetName = viewer.datasetName;
    }


    public void doAction(final IntegerList list, final TreePath classTreePath) {
//        String base = new File(filePathRoot).getParent() + File.separator +
//                      experimentName + "-" + datasetName.substring(4) + "-";
        String base = new File(filePathRoot).getParent() + File.separator + 
        datasetName + "-";

        String htmlFilePath = base + "report.html";
        String htmlImagesPath = base + "img";

        File imagesDirectory = new File(htmlImagesPath);
        String htmlImagesPartialPath = imagesDirectory.getName();
        imagesDirectory.mkdir();

        System.out.println("htmlFilePath = " + htmlFilePath);
        System.out.println("htmlImagesPath = " + htmlImagesPath);
        System.out.println("htmlImagesPartialPath = " + htmlImagesPartialPath);

        try {
            String NEWLINE = System.getProperty("line.separator");

            PrintWriter html = new PrintWriter(new FileWriter(htmlFilePath), false);
            html.write("<html>" + NEWLINE);
            html.write("<head>" + NEWLINE);
            html.write("<title>Title</title>" + NEWLINE);
            html.write("</head>" + NEWLINE);
            html.write("<body>" + NEWLINE);

            //html.write("<h1>" + experimentName + "/" + datasetName + "</h1>");
            html.write("<h1>" + filePathRoot + ".params" + "</h1>");
            
            String a = experimentName + "-" + datasetName.substring(4) +
                       "-classification.html";
            html.write("<a href=\"" + a + "\"> Classification </a>");
            html.write("<br><br>");

            html.write("<table border=1 cellpadding=4 cellspacing=0><tr>" +
                       "<td>Class");
            for (int i = 0; i < expressions.size(); i++) {
                html.write("<td>" + expressions.get(i).name);
            }
            /*
                        html.write("<table border=1 cellpadding=4 cellspacing=0><tr>" +
                                   "<td>Class" +
                                   "<td>Radius" +
                                   "<td>Spikes" +
                                   "<td>t1" +
                                   "<td>t2" +
                                   "<td>a1" +
                                   "<td>a2" +
                                   "<td>n1" +
                                   "<td>n2"
                            );
             */

            DefaultMutableTreeNode root = (DefaultMutableTreeNode) classTreePath.
                                          getLastPathComponent();
            writeFolder(root, htmlImagesPath, htmlImagesPartialPath, html);

            html.write("</table>");

            html.write("</body>" + NEWLINE);
            html.flush();
            html.close();
        } catch (IOException e) {
            Vision.reportException(e);
        }
    }


    public void writeFolder(DefaultMutableTreeNode folder, String imagesPath,
                            String imagesPartialPath, PrintWriter indexHTML) throws
        IOException {

        // Generate full path from current position in tree
        String folderName = folder.toString();
        DefaultMutableTreeNode currentFolder = folder;
        boolean done = folder.isRoot();
        while (!done) {
            currentFolder = (DefaultMutableTreeNode) currentFolder.getParent();
            folderName = currentFolder.toString() + "/" + folderName;
            done = currentFolder.isRoot();
        }
        int[] neuronsInClass = InteractiveTree.getNeuronsInClass(
            folder, new IntegerList(), false).toArray();

        if (neuronsInClass.length != 0) {
            if (folderName.equals("All")) {
                System.err.println("You have to classify the " + neuronsInClass.length +
                                   " neurons under All before a report can be written.");
                return;
            }

            MeanVarianceCalculator mvc = new MeanVarianceCalculator();

            // write the class name
            String cName = folderName.substring(4).replace('/', ' ');
            indexHTML.write("<tr>");
            // Class
            indexHTML.write("<td> <b><a href=\"" + imagesPartialPath + "/" + cName +
                            ".html" + "\">" + cName + "</a></b>");

            // save all the statistics
            for (int i = 0; i < expressions.size(); i++) {
                Expression e = expressions.get(i);
                calculateStatistics(e.expression, folderName, mvc);
                indexHTML.write(
                    "<td>" +
                    StringUtil.format(mvc.getMean(), e.nDigits, 1) + " &#177; " +
                    StringUtil.format(mvc.getMeanVariance(), e.nDigits, 1));
            }

            //
            PrintWriter classHTML = new PrintWriter(imagesPath + "/" + cName + ".html");
            classHTML.write("<html><head><title>" + cName + "</title></head><body>");
            classHTML.write("<h2>Class: " + cName + "</h2>");

            for (int i = 0; i < neuronsInClass.length; i++) {
                // save the neuron plot
                viewer.rightTree.setSelectionPath(
                    new TreePath(folder.getPath()).pathByAddingChild(
                        new DefaultMutableTreeNode(new Integer(neuronsInClass[i]), false)
                    ));
                String iName = neuronsInClass[i] + ".png";
                String gifName = neuronsInClass[i] + ".gif";

                GraphicsIO.saveComponentToPNG(viewer.rightSplitPane,
                        imagesPath + File.separator + iName);
                
                //find ei and save it out as animated gif.  This is very slow.
                for(int j=0; j<viewer.rightSplitPane.getComponentCount(); j++) {
                    Component comp = viewer.rightSplitPane.getComponent(j);
                    if(comp instanceof JPanel){
                        JPanel jPanel = (JPanel) comp;
                        for(int k=0; k<jPanel.getComponentCount(); k++) {
                            Component comp2 = jPanel.getComponent(k);
                            if(comp2 instanceof PhysiologicalImagePanel) {
                                PhysiologicalImagePanel eiPanel = (PhysiologicalImagePanel) comp2;
                                eiPanel.makeGIF(imagesPath + File.separator + gifName);
                            }
                        }
                    }
                }

                classHTML.write(
                    "<a href=\"" + iName + "\">" + neuronsInClass[i] + "</a>, ");
                
                classHTML.write(
                        "<a href=\"" + gifName + "\">" + "ei" + "</a>, ");
            }

            classHTML.write("<hr>");
            classHTML.write("<img src=\"" + cName + ".png" + "\"> </img>");
            classHTML.write("</body></html>");
            classHTML.close();

            // save the class plot
            viewer.rightTree.setSelectionPath(new TreePath(folder.getPath()));
            GraphicsIO.saveComponentToPNG(viewer.rightSplitPane,
                                          imagesPath + File.separator + cName + ".png");
        }

        for (int i = 0; i < folder.getChildCount(); i++) {
            if (folder.getChildAt(i).getAllowsChildren()) {
                writeFolder( (DefaultMutableTreeNode) folder.getChildAt(i),
                            imagesPath, imagesPartialPath, indexHTML);
            }
        }
    }


    private void calculateStatistics(String expression, String classPath,
                                     MeanVarianceCalculator mvc) {

        mvc.reset();

        HashMap<Integer, Double> valueMap = null;
        try {
            valueMap = paramsFile.evaluate(
                expression, "classID==\"" + classPath + "\"&" + "true");
        } catch (CannotEvaluateException ex) {
            return;
        }

        for (Integer id : valueMap.keySet()) {
            Double v = valueMap.get(id);
            if (v != null) {
                mvc.add(v);
            }
        }
    }

}
