package edu.ucsc.neurobiology.vision.test;

import java.io.*;


/**
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public final class diff {

    public static void main(String[] args) throws IOException {
        File file1 = new File("C:\\Program Files\\WildTangent\\Apps\\GameChannel\\Games\\C2D8F0E2-6978-4409-8351-BA8785DA11EE\\Fate.exe");
        File file2 = new File("C:\\fate\\Fate.exe");
        DataInputStream dis1 = new DataInputStream(new FileInputStream(file1));
        DataInputStream dis2 = new DataInputStream(new FileInputStream(file2));

        long size1 = file1.length();
        long size2 = file1.length();
        if (size1 != size2) {
            System.err.println("different file sizes.");
            System.exit(1);
        }

        long i = 0;
        while (i < size1) {
            float b1 = dis1.readByte();
            float b2 = dis2.readByte();
            if (b1 != b2) {
                System.out.println(
                    "different at byte " + i + " : " + b1 + " : " + b2);
                //   System.exit(1);
            }
            i++;
        }

        System.out.println("equal");

        dis1.close();
        dis2.close();
    }
}
