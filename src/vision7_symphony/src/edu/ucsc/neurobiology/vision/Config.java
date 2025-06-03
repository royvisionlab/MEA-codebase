package edu.ucsc.neurobiology.vision;

import java.io.*;
import java.util.*;
import javax.xml.parsers.*;
import javax.xml.transform.*;
import javax.xml.transform.dom.*;
import javax.xml.transform.stream.*;

import java.awt.*;

import org.w3c.dom.*;
import org.xml.sax.*;

import edu.ucsc.neurobiology.vision.parameters.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * This class handles the configuration information stored in the <tt>config.xml</tt>
 * file. The various properties of the config file can be accessed using corresponding
 * methods of this class.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class Config
    implements ErrorHandler {

    private Document document;
    private Element configuration;
    private File file;
    public String fileName;


    /**
     * Creates a new <tt>Config</tt> object and initializes it with the information from
     * the given configuration file.
     */
    public Config(String fn) throws IOException, SAXException,
        ParserConfigurationException {

        fileName = fn;
        file = new File(fileName);
        if ( !file.exists() ) {
            throw new IllegalArgumentException("Config file " + fileName + " does not exist.");
        }
        if ( !file.canRead() ) {
            throw new IllegalArgumentException("Cannot read the config file: " + fileName);
        }

        DocumentBuilder documentBuilder = makeDocBuilder();
        document = documentBuilder.parse(file);
        configuration = document.getDocumentElement();
    }
    
    public Config(InputStream configStream) throws ParserConfigurationException, SAXException, IOException {
        this(configStream, "");
    }
    
    public Config(InputStream configStream, String fn) throws ParserConfigurationException, SAXException, IOException {
        fileName = fn;
        DocumentBuilder documentBuilder = makeDocBuilder();
        document = documentBuilder.parse(configStream);
        configuration = document.getDocumentElement();    
    }
    
    private DocumentBuilder makeDocBuilder() throws ParserConfigurationException {
        DocumentBuilderFactory documentFactory = DocumentBuilderFactory.newInstance();
        documentFactory.setIgnoringComments(true);
        documentFactory.setIgnoringElementContentWhitespace(true);
        documentFactory.setValidating(false);
        
        DocumentBuilder documentBuilder = documentFactory.newDocumentBuilder();
        documentBuilder.setErrorHandler(this);
        return documentBuilder;
    }


    /** Receive notification of a recoverable error, for internal use, Do NOT call. */
    public void error(SAXParseException e) {
        Vision.reportFatalException(e);
    }


    /** Receive notification of a non-recoverable error, for internal use, Do NOT call. */
    public void fatalError(SAXParseException e) {
        Vision.reportFatalException(e);
    }


    /** Receive notification of a warning, for internal use, Do NOT call. */
    public void warning(SAXParseException e) {
        Vision.reportFatalException(e);
    }


    public void save() throws Exception {
        TransformerFactory xformFactory = TransformerFactory.newInstance();
        Transformer idTransform = xformFactory.newTransformer();
        idTransform.setOutputProperty("indent", "yes");
        idTransform.setOutputProperty("{http://xml.apache.org/xalan}indent-amount", "4");
        document.normalize();

        FileWriter w = new FileWriter(file);
        idTransform.transform(new DOMSource(document), new StreamResult(w));
        w.close();
    }


    private String getAtributeForNode(Node node, String atributeName) {
        NamedNodeMap attributes = node.getAttributes();
        return attributes.getNamedItem(atributeName).getNodeValue();
    }


    public Node getGroupNode(String groupName) {
        String[] names = StringUtil.decomposeString(groupName, ".");
        Node currentNode = configuration;

        for (int k = 0; k < names.length; k++) {
            NodeList currentNodeList = currentNode.getChildNodes();

            for (int i = 0; i < currentNodeList.getLength(); i++) {
                Node groupNode = currentNodeList.item(i);
                if (groupNode.getNodeType() == Node.ELEMENT_NODE &&
                    getAtributeForNode(groupNode, "name").equals(names[k])) {

                    currentNode = groupNode;
                }
            }
//            if (currentNode == null) {
//                return null;
//            }
        }

        return currentNode;
    }


    public ParametersTable showDialog(String groupName, String title, Component parent) {
        ParametersTable table = getParameterGroup(groupName);
        if (table == null)
            throw new IllegalArgumentException( "Sorry, cannot find parameter group: " + groupName);
        table.expandAll();
        
        ParametersDialog dialog = new ParametersDialog(parent, title, table);
        dialog.setVisible(true);
        
        if (dialog.isOkSelected() || dialog.isSaveSelected()) {
            try {
                setParameterGroup(groupName, table);
            } catch (Exception e) {
                Vision.reportException("Cannot save parameters", e);
            }
        }
        
        dialog.getContentPane().removeAll();
        dialog.dispose();
        
        if (dialog.isOkSelected()) {
            return table;
        } else {
            return null;
        }
    }


    public ParametersTable getParameterGroup(String groupName) {
        Node node = getGroupNode(groupName);
        if (node == null) return null;
        return new ParametersTable(node);
    }
    

    public LinkedHashMap<String, String> getParameterList(String groupName) {
        Node node = getGroupNode(groupName);
        if (node == null) return null;
        return getParameterList(node, new LinkedHashMap<String, String>(), "");
    }
    

    /**
     *
     * @param groupName String
     * @return LinkedHashMap
     */
    private LinkedHashMap<String, String> getParameterList(Node node, LinkedHashMap<String, String> params,
                                           String longName) {
        NodeList childNodeList = node.getChildNodes();
        for (int i = 0; i < childNodeList.getLength(); i++) {
            Node childNode = childNodeList.item(i);
            if (childNode.getNodeType() != Node.ELEMENT_NODE) {
                continue;
            }
            NamedNodeMap attributes = childNode.getAttributes();
            String name = attributes.getNamedItem("name").getNodeValue().trim();
            String value = attributes.getNamedItem("value").getNodeValue().trim();

            params.put(longName + name, value);

            // go into subgroups
            if (childNode.getNodeName().equals("ParametersGroup")) {
                getParameterList(childNode, params, longName + name + ".");
            }
        }

        return params;
    }


    public void setParameterGroup(String groupName, HashMap<String, String> params) throws Exception {
        Node node = getGroupNode(groupName);
        if (node == null) {
            return;
        }
        setParameterGroup(node, params, "");
    }


    private void setParameterGroup(Node node, HashMap<String, String> params, String longName) throws
        Exception {

        NodeList childNodeList = node.getChildNodes();
        for (int i = 0; i < childNodeList.getLength(); i++) {
            Node childNode = childNodeList.item(i);
            if (childNode.getNodeType() != Node.ELEMENT_NODE) {
                continue;
            }
            NamedNodeMap attributes = childNode.getAttributes();
            String name = attributes.getNamedItem("name").getNodeValue().trim();

            attributes.getNamedItem("value").setNodeValue(
                (String) params.get(longName + name));

            // go into subgroups
            if (childNode.getNodeName().equals("ParametersGroup")) {
                setParameterGroup(childNode, params, longName + name + ".");
            }
        }

        save();
    }


    public void setParameterGroup(String groupName, ParametersTable params) throws
        Exception {

        Node node = getGroupNode(groupName);
        if (node == null) {
            return;
        }
        setParameterGroup(node, params, "");
    }


    private void setParameterGroup(Node node, ParametersTable params, String longName) throws
        Exception {

        NodeList childNodeList = node.getChildNodes();
        for (int i = 0; i < childNodeList.getLength(); i++) {
            Node childNode = childNodeList.item(i);
            if (childNode.getNodeType() != Node.ELEMENT_NODE) {
                continue;
            }
            NamedNodeMap attributes = childNode.getAttributes();
            String name = attributes.getNamedItem("name").getNodeValue().trim();

            Parameter p = params.getParameter(longName + name);
            if (p == null) {
                throw new Error("Parameter missing: " + longName + name);
            }
            attributes.getNamedItem("value").setNodeValue(p.valueAsString());

            // go into subgroups
            if (childNode.getNodeName().equals("ParametersGroup")) {
                setParameterGroup(childNode, params, longName + name + ".");
            }
        }

        save();
    }


}
