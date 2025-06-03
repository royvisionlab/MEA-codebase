package edu.ucsc.neurobiology.vision.util;

import java.io.*;


/**
 * A table of numbers that supports printing in a table-like format.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class Table {
    private final int nRows, nCols;
    private String[][] cells;
    private String[] titles;


    public Table(final int nRows, final int nCols) {
        this.nRows = nRows;
        this.nCols = nCols;
        cells = new String[nRows][nCols];
        titles = new String[nCols];
    }


    public void setCell(int row, int col, Object data) {
        cells[row][col] = data.toString();
    }


    public void setTitle(int col, Object title) {
        titles[col] = "(" + title.toString() + ")";
    }


    public void excelDraw(PrintStream out) {

        for (int i = 0; i < titles.length; i++) {
            out.print(titles[i] + "\t");
        }
        out.print("\r\n");
        for (int i = 0; i < cells.length; i++) {
            for (int j = 0; j < cells[0].length; j++) {
                out.print(cells[i][j] + "\t");
            }
            out.print("\r\n");
        }
    }


    public void draw(PrintStream out) {
        int[] width = new int[nCols];

        for (int col = 0; col < nCols; col++) {
            int maxWidth = titles[col].length();
            for (int row = 0; row < nRows; row++) {
                if (cells[row][col] == null) {
                    cells[row][col] = "?";
                }
                if (cells[row][col].length() > maxWidth) {
                    maxWidth = cells[row][col].length();
                }
            }
            width[col] = maxWidth;
        }

        for (int col = 0; col < nCols; col++) {
            titles[col] = justify(titles[col], width[col]);
            for (int row = 0; row < nRows; row++) {
                cells[row][col] = justify(cells[row][col], width[col]);
            }
        }

        String titleString = "";
        for (int col = 0; col < nCols; col++) {
            titleString += titles[col] + " ";
        }
        out.println(titleString);
        for (int row = 0; row < nRows; row++) {
            String rowString = "";
            for (int col = 0; col < nCols; col++) {
                rowString += cells[row][col] + " ";
            }
            out.println(rowString);
        }
    }


    private String justify(String s, int newLength) {
        int length = s.length();

        if (length > newLength) {
            throw new IllegalArgumentException("s.length() > newLength");
        } else if (length == newLength) {
            return s;
        } else {
            int dLength = newLength - length;
            for (int i = 0; i < dLength; i++) {
                s = " " + s;
            }
            return s;
        }
    }

}
