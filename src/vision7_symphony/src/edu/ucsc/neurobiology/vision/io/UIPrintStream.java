package edu.ucsc.neurobiology.vision.io;

import java.io.*;

import javax.swing.*;


/**
 * This class acts as a PrintStream. It displays all output written to it in a text area.
 * Used by Vision to provide a graphical output.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class UIPrintStream
    extends PrintStream {

    JTextArea area;
    String separator;


    public UIPrintStream(OutputStream os, JTextArea area) {
        super(os, true);
        this.area = area;
        this.area.setEditable(false);
        separator = System.getProperty("line.separator");
    }


    public void print(char c) {
        super.print(c);
        area.append(String.valueOf(c));
        scroll();
    }


    public void print(char[] s) {
        super.print(s);
        area.append(String.valueOf(s));
        scroll();
    }


    public void print(double d) {
        super.print(d);
        area.append(String.valueOf(d));
        scroll();
    }


    public void print(float f) {
        super.print(f);
        area.append(String.valueOf(f));
        scroll();
    }


    public void print(int i) {
        super.print(i);
        area.append(String.valueOf(i));
        scroll();
    }


    public void print(long l) {
        super.print(l);
        area.append(String.valueOf(l));
        scroll();
    }


    public void print(Object obj) {
        super.print(obj);
        area.append(String.valueOf(obj));
        scroll();
    }


    public void print(String s) {
        super.print(s);
        area.append(s);
        scroll();
    }


    public void println() {
        super.println();
        area.append(separator);
        scroll();
    }


    public void println(char c) {
        super.println(c);
        area.append(separator);
        scroll();
    }


    public void println(char[] s) {
        super.println(s);
        area.append(separator);
        scroll();
    }


    public void println(double d) {
        super.println(d);
        area.append(separator);
        scroll();
    }


    public void println(float f) {
        super.println(f);
        area.append(separator);
        scroll();
    }


    public void println(int i) {
        super.println(i);
        area.append(separator);
        scroll();
    }


    public void println(long l) {
        super.println(l);
        area.append(separator);
        scroll();
    }


    public void println(Object obj) {
        super.println(obj);
        area.append(separator);
        scroll();
    }


    public void println(String s) {
        super.println(s);
        area.append(separator);
        scroll();
    }


    private void scroll() {
        area.setCaretPosition(area.getText().length());
    }
}
