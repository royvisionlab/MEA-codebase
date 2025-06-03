package edu.ucsc.neurobiology.vision.test;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.util.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class linecount {
    static LinkedList<FileEntry> list = new LinkedList<FileEntry>();


    public static void main(String args[]) throws Exception {
        if (args.length != 1) {
            System.out.println("Wrong arguments");
            System.exit(1);
        }
        File directory = new File(args[0]);
        if (!directory.exists()) {
            System.out.println(args[0] + " does not exist");
            System.exit(1);
        } else if (!directory.isDirectory()) {
            System.out.println(args[0] + " is not a directory");
            System.exit(1);
        }
        count(directory.getAbsoluteFile());
        Collections.sort(list);
        int totalLines = 0;
        Table table = new Table(list.size() + 2, 2);
        table.setTitle(1, "File Name");
        table.setTitle(0, "Lines");
        int i = 0;
        for (FileEntry e : list) {
            totalLines += e.nLines;
            //System.out.println(e.name + ":  " + e.nLines);
            table.setCell(i, 1, e.name);
            table.setCell(i, 0, "" + e.nLines);
            i++;
        }
        table.setCell(list.size(), 0, "======");
        table.setCell(list.size(), 1, "=============");
        table.setCell(list.size() + 1, 0, "" + totalLines);
        table.setCell(list.size() + 1, 1, list.size() + " files");
        table.draw(System.err);
    }


    public static void count(File directory) throws IOException {
        File[] childs = directory.listFiles();

        for (int i = 0; i < childs.length; i++) {
            String name = childs[i].getName();
            if (childs[i].isFile() && name.endsWith(".java")) {
                int nLines = countLines(childs[i]);
                list.add(new FileEntry(name, nLines));
            } else if (childs[i].isDirectory()) {
                count(childs[i]);
            }
        }
    }


    public static int countLines(File file) throws IOException {
        LineNumberReader r = new LineNumberReader(new FileReader(file));

        int n;
        for (n = 1; r.readLine() != null; n++) {
            ;
        }

        return n;
    }

}


class FileEntry
    implements Comparable<FileEntry> {
    int nLines;
    String name;

    public FileEntry(String name, int nLines) {
        this.nLines = nLines;
        this.name = name;
    }


    public int compareTo(FileEntry o) {
        FileEntry e = o;
        if (nLines > e.nLines) {
            return 1;
        } else if (nLines < e.nLines) {
            return -1;
        } else {
            return 0;
        }
    }
}
