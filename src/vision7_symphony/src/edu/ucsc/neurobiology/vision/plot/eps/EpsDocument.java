/*
 Copyright Paul James Mutton, 2001-2004, http://www.jibble.org/

 This file is part of EpsGraphics2D.

 This software is dual-licensed, allowing you to choose between the GNU
 General Public License (GPL) and the www.jibble.org Commercial License.
 Since the GPL may be too restrictive for use in a proprietary application,
 a commercial license is also provided. Full license information can be
 found at http://www.jibble.org/licenses/

 @author Matthew Grivich, University of California, Santa Cruz
 @author Dumitru Petrusca, University of California, Santa Cruz
 */
package edu.ucsc.neurobiology.vision.plot.eps;


import java.io.*;
import java.util.*;


/**
 * This represents an EPS document. Several EpsGraphics2D objects may point
 * to the same EpsDocument.
 *  <p>
 * Copyright Paul Mutton,
 *           <a href="http://www.jibble.org/">http://www.jibble.org/</a>
 *
 */
public class EpsDocument {

    private float minX, _minX;
    private float minY, _minY;
    private float maxX, _maxX;
    private float maxY, _maxY;

    private boolean _isClipSet = false;

    private String _title;
    private StringWriter _stringWriter;
    private BufferedWriter _bufferedWriter = null;

    // We need to remember which was the last EpsGraphics2D object to use
    // us, as we need to replace the clipping region if another EpsGraphics2D
    // object tries to use us.
    private EpsGraphics2D _lastG = null;


    /**
     * Constructs an empty EpsDevice that writes directly to a file.
     * Bounds must be set before use.
     */
    public EpsDocument(String title, int minX, int minY, int maxX, int maxY) {
        _title = title;
        this.minX = this._minX = minX;
        this.minY = this._minY = minY;
        this.maxX = this._maxX = maxX;
        this.maxY = this._maxY = maxY;
        _stringWriter = new StringWriter();
        _bufferedWriter = new BufferedWriter(_stringWriter);
    }


    /**
     * Returns the title of the EPS document.
     */
    public synchronized String getTitle() {
        return _title;
    }


    /**
     * Updates the bounds of the current EPS document.
     */

    public synchronized void updateBounds(double x, double y) {
        if (x > maxX) {
            maxX = (float) x;
        }
        if (x < minX) {
            minX = (float) x;
        }
        if (y > maxY) {
            maxY = (float) y;
        }
        if (y < minY) {
            minY = (float) y;
        }
    }


    /**
     * Appends a line to the EpsDocument.  A new line character is added
     * to the end of the line when it is added.
     */
    public synchronized void append(EpsGraphics2D g, String line) {
        if (_lastG == null) {
            _lastG = g;
        } else if (g != _lastG) {
            EpsGraphics2D lastG = _lastG;
            _lastG = g;
            // We are being drawn on with a different EpsGraphics2D context.
            // We may need to update the clip, etc from this new context.
            if (g.getClip() != lastG.getClip()) {
                g.setClip(g.getClip());
            }
            if (!g.getColor().equals(lastG.getColor())) {
                g.setColor(g.getColor());
            }
            if (!g.getBackground().equals(lastG.getBackground())) {
                g.setBackground(g.getBackground());
            }
            // We don't need this, as this only affects the stroke and font,
            // which are dealt with separately later on.
            //if (!g.getTransform().equals(lastG.getTransform())) {
            //    g.setTransform(g.getTransform());
            //}
            if (!g.getPaint().equals(lastG.getPaint())) {
                g.setPaint(g.getPaint());
            }
            if (!g.getComposite().equals(lastG.getComposite())) {
                g.setComposite(g.getComposite());
            }
            if (!g.getComposite().equals(lastG.getComposite())) {
                g.setComposite(g.getComposite());
            }
            if (!g.getFont().equals(lastG.getFont())) {
                g.setFont(g.getFont());
            }
            if (!g.getStroke().equals(lastG.getStroke())) {
                g.setStroke(g.getStroke());
            }
        }
        _lastG = g;

        try {
            if(!line.contains("NaN")) //Do not add NaNs to file.  Added by Matthew Grivich, The Salk Institute
                _bufferedWriter.write(line + "\n");
        } catch (IOException e) {
            throw new EpsException("Could not write to the output file: " + e);
        }
    }


    /**
     * Outputs the contents of the EPS document to the specified
     * Writer, complete with headers and bounding box.
     */
    public synchronized void write(Writer writer) throws IOException {
        float offsetX = -minX;
        float offsetY = -minY;

        writer.write("%!PS-Adobe-3.0 EPSF-3.0\n");
        writer.write("%%Creator: EpsGraphics2D " + EpsGraphics2D.VERSION +
                     " by Paul Mutton, http://www.jibble.org/\n");
        writer.write("%%Title: " + _title + "\n");
        writer.write("%%CreationDate: " + new Date() + "\n");

        writer.write("%%BoundingBox: 0 0 " +
                     ( (int) Math.ceil(_maxX - _minX + 1)) + " " +
                     ( (int) Math.ceil(_maxY - _minY + 1)) + "\n");

        writer.write("%%DocumentData: Clean7Bit\n");
        writer.write("%%DocumentProcessColors: Black\n");
        writer.write("%%ColorUsage: Color\n");
        writer.write("%%Origin: 0 0\n");
        writer.write("%%Pages: 1\n");
        writer.write("%%Page: 1 1\n");
        writer.write("%%EndComments\n\n");

        writer.write("gsave\n");

        if (_stringWriter != null) {
            writer.write( (_minX + 0.5) + " " + (_maxY + 0.5) + " translate\n");

            _bufferedWriter.flush();
            StringBuffer buffer = _stringWriter.getBuffer();
            for (int i = 0; i < buffer.length(); i++) {
                writer.write(buffer.charAt(i));
            }

            writeFooter(writer);
        } else {
            writer.write(offsetX + " " + ( (maxY - minY) - offsetY) + " translate\n");
        }

        writer.flush();
    }


    private void writeFooter(Writer writer) throws IOException {
        writer.write("grestore\n");
        if (isClipSet()) {
            writer.write("grestore\n");
        }
        writer.write("showpage\n");
        writer.write("\n");
        writer.write("%%EOF");
        writer.flush();
    }


    public synchronized void flush() throws IOException {
        _bufferedWriter.flush();
    }


    public synchronized void close() throws IOException {
        if (_stringWriter == null) {
            writeFooter(_bufferedWriter);
            _bufferedWriter.flush();
            _bufferedWriter.close();
        }
    }


    public boolean isClipSet() {
        return _isClipSet;
    }


    public void setClipSet(boolean isClipSet) {
        _isClipSet = isClipSet;
    }

}
