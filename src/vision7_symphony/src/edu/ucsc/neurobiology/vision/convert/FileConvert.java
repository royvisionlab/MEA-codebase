package edu.ucsc.neurobiology.vision.convert;

import java.io.*;
import java.util.*;

import java.awt.geom.*;

import edu.ucsc.neurobiology.vision.electrodemap.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class FileConvert {
    public static final float standardSigmaNoise = 20;


    public static void oldFormat2newFormat(String inputFileName, String outputFileName) throws
        IOException {

        final SynchronizationObject sync = new SynchronizationObject();
        File inputFile = new File(inputFileName);
        FileInputStream inputStream = new FileInputStream(inputFile);
        RawDataHeader65 header65 = new RawDataHeader65(inputStream);
        AsynchronousInputStream ais =
            new AsynchronousInputStream(inputStream, 100 * 1024, 1000);
        final int nElectrodes = header65.getNumberOfElectrodes();
        final long nSamples =
            (inputFile.length() - header65.getHeaderSize()) / (2 * nElectrodes);
        SampleInputStream sis = new SampleInputStream(ais, 500, nElectrodes);

        FileOutputStream outputStream = new FileOutputStream(outputFileName);
        RawDataHeader512 header512 = new RawDataHeader512(
            (int) (new Date().getTime() / 1000),
            nElectrodes,
            20000,
            (int) nSamples,
            0,
            RawDataHeader.FORMAT_12BIT_COMPRESSED,
            "", // FIXME, wrong dataset identifier
            "compressed");
        outputStream.write(header512.getBinaryRepresentation());
        final BufferedOutputStream bOutputStream =
            new BufferedOutputStream(outputStream, 1024 * 1024);

        sis.addDataListener(new SampleListener() {
            long sampleIndex, percentage, oldPercentage;

            public void processSample(short[] sample) {
                try { //BEU
                    bOutputStream.write( (byte) (sample[0] >> 8));
                    bOutputStream.write( (byte) (sample[0] & 0x00ff));

                    for (int i = 1; i < nElectrodes; ) {
                        int s1 = sample[i++] + 2048;
                        int s2 = sample[i++] + 2048;

                        bOutputStream.write(s1 >> 4);
                        bOutputStream.write( ( (s1 & 0x000f) << 4) + (s2 >> 8));
                        bOutputStream.write(s2 & 0x00ff);
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }

                sampleIndex++;
                percentage = 100 * sampleIndex / nSamples;
                if (percentage != oldPercentage) {
                    System.out.println(percentage + "%");
                    oldPercentage = percentage;
                }
            }


            public void finishSampleProcessing() throws IOException {
                bOutputStream.close();

                sync.done();
            }
        });

        sync.setWorking();

        ais.start();
        sis.start();

        sync.waitUntilDone();
    }


    /**
     * Converts a 512 format compressed file into a 64 format uncompressed file.
     * An optional old-type header can be included.
     */
    public static void newFormat2oldFormat(
        String inputFileName, String outputFileName, boolean includeHeader) throws
        IOException {

        final SynchronizationObject sync = new SynchronizationObject();
        MultipleCompressedSampleInputStream sis =
            new MultipleCompressedSampleInputStream(inputFileName);
        RawDataHeader header512 = sis.getHeader();
        final int nElectrodes = header512.getNumberOfElectrodes();
        final int nSamples = header512.getNumberOfSamples();
        FileOutputStream outputStream = new FileOutputStream(outputFileName);
        if (includeHeader) {
            RawDataHeader header65 = new RawDataHeader65(nElectrodes,
                header512.getSamplingFrequency());
            byte[] binHeader = header65.getBinaryRepresentation();
            outputStream.write(binHeader);
            System.out.println("Header Length : " + binHeader.length);
        }
        final BufferedOutputStream bOutputStream =
            new BufferedOutputStream(outputStream, 10 * 1024 * 1024);
        sis.addSampleListener(new SampleListener() {
            long sampleIndex, percentage, oldPercentage = -1;
            public void processSample(short[] sample) {
                try { //BEU
                    for (int i = 0; i < nElectrodes; i++) {
                        bOutputStream.write( (byte) (sample[i] >> 8));
                        bOutputStream.write( (byte) (sample[i] & 0x00ff));
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }
                sampleIndex++;
                percentage = 100 * sampleIndex / nSamples;
                if (percentage != oldPercentage) {
                    System.out.println(percentage + "%");
                    oldPercentage = percentage;
                }
            }


            public void finishSampleProcessing() throws IOException {
                bOutputStream.close();
                sync.done();
            }
        });

        sync.setWorking();

        sis.start();

        sync.waitUntilDone();
        System.out.println("now exit the method");
    }


    /*
        public static void filterAndCorrectGain(
            String inputFileName, String outputFileName, String sigmasFileName) throws
            IOException {
            final SynchronizationObject sync = new SynchronizationObject();
            sync.setWorking();
            MultipleCompressedSampleInputStream sis =
                new MultipleCompressedSampleInputStream(inputFileName, 500);
            RawDataHeader inputHeader = sis.getHeader();
            final int nElectrodes = inputHeader.getNumberOfElectrodes();
            final int nSamples = inputHeader.getNumberOfSamples();
            final FileOutputStream outputStream = new FileOutputStream(outputFileName);
            final BufferedOutputStream bOutputStream =
                new BufferedOutputStream(outputStream, 100 * 1024 * 1024);
            bOutputStream.write(inputHeader.getBinaryRepresentation());
            final FIRFilter[] filter = new FIRFilter[nElectrodes];
            for (int i = 0; i < nElectrodes; i++) {
                filter[i] = new FIRFilter();
            }
            final float[] factor = SpikeFinding.getSigmas(sigmasFileName, nElectrodes);
            for (int i = 0; i < nElectrodes; i++) {
                factor[i] = standardSigmaNoise / factor[i];
            }
            sis.addSampleListener(new SampleListener() {
                long sampleIndex = 0;
                long percentage, oldPercentage = -1;
                byte[] buffer = new byte[2 + (nElectrodes - 1) * 3 / 2];
                int byteIndex;
                public void processSample(short[] sample) {
                    sampleIndex++;
                    byteIndex = 0;
                    buffer[byteIndex++] = (byte) (sample[0] >> 8);
                    buffer[byteIndex++] = (byte) (sample[0] & 0x00ff);
                    for (int i = 1; i < nElectrodes; ) {
                        int s1 = (int) (filter[i].filter(sample[i]) * factor[i]) + 2048;
                        i++;
                        int s2 = (int) (filter[i].filter(sample[i]) * factor[i]) + 2048;
                        i++;
                        buffer[byteIndex++] = (byte) (s1 >> 4);
         buffer[byteIndex++] = (byte) ( ( (s1 & 0x000f) << 4) + (s2 >> 8));
                        buffer[byteIndex++] = (byte) (s2 & 0x00ff);
                    }
                    try {
                        bOutputStream.write(buffer);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    percentage = 100 * sampleIndex / nSamples;
                    if (percentage != oldPercentage) {
                        System.out.println(percentage + "%");
                        oldPercentage = percentage;
                    }
                }
                public void finishProcessing() throws IOException {
                        bOutputStream.close();
                        outputStream.close();
                    sync.done();
                }
            });
            sis.start();
            sync.waitUntilDone();
        }
     */


    /*
        public static void filter(
     String inputFileName, String outputFileName, int timeInSeconds) throws IOException {
            final SynchronizationObject sync = new SynchronizationObject();
            sync.setWorking();
            MultipleCompressedSampleInputStream sis =
                new MultipleCompressedSampleInputStream(inputFileName, 500);
            RawDataHeader inputHeader = sis.getHeader();
            final int nElectrodes = inputHeader.getNumberOfElectrodes();
            final int nSamples = inputHeader.getNumberOfSamples();
            final FileOutputStream outputStream = new FileOutputStream(outputFileName);
            final BufferedOutputStream bOutputStream =
                new BufferedOutputStream(outputStream, 100 * 1024 * 1024);
            bOutputStream.write(inputHeader.getBinaryRepresentation());
            sis.addSampleListener(new SampleListener() {
                long sampleIndex = 0;
                long percentage, oldPercentage = -1;
                byte[] buffer = new byte[2 + (nElectrodes - 1) * 3 / 2];
                int byteIndex;
                public void processSample(short[] sample) {
                    sampleIndex++;
                    byteIndex = 0;
                    buffer[byteIndex++] = (byte) (sample[0] >> 8);
                    buffer[byteIndex++] = (byte) (sample[0] & 0x00ff);
                    for (int i = 1; i < nElectrodes; ) {
                        int s1 = (int) filter[i].filter(sample[i]) + 2048;
                        i++;
                        int s2 = (int) filter[i].filter(sample[i]) + 2048;
                        i++;
                        buffer[byteIndex++] = (byte) (s1 >> 4);
         buffer[byteIndex++] = (byte) ( ( (s1 & 0x000f) << 4) + (s2 >> 8));
                        buffer[byteIndex++] = (byte) (s2 & 0x00ff);
                    }
                    try {
                        bOutputStream.write(buffer);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    percentage = 100 * sampleIndex / nSamples;
                    if (percentage != oldPercentage) {
                        System.out.println(percentage + "%");
                        oldPercentage = percentage;
                    }
                }
                public void finishProcessing() throws IOException {
                        bOutputStream.close();
                        outputStream.close();
                    sync.done();
                }
            });
            sis.start();
            sync.waitUntilDone();
        }
     */

    /*
        public static void filter1(String inputFileName, String outputFileName,
                                   double timeConstant) throws IOException {
            final double t1 = System.currentTimeMillis();
            MultipleCompressedSampleInputStream sis =
                new MultipleCompressedSampleInputStream(inputFileName, 500);
            RawDataHeader inputHeader = sis.getHeader();
            final int nElectrodes = inputHeader.getNumberOfElectrodes();
            final int nSamples = inputHeader.getNumberOfSamples();
            FileOutputStream outputStream = new FileOutputStream(outputFileName);
            outputStream.write(inputHeader.getBinaryRepresentation());
            final BufferedOutputStream bOutputStream =
                new BufferedOutputStream(outputStream, 1024 * 1024);
            final double[] mean = new double[nElectrodes];
            final double alpha = 1.0 / (timeConstant * 20.0);
            sis.addSampleListener(new SampleListener() {
                long sampleIndex, percentage, oldPercentage = -1;
                public void processSample(short[] sample) {
                    try {
                        // write the TTL electrode
                        bOutputStream.write( (byte) (sample[0] >> 8));
                        bOutputStream.write( (byte) (sample[0] & 0x00ff));
                        for (int i = 1; i < nElectrodes; ) {
                            mean[i] += alpha * (sample[i] - mean[i]);
                            int s1 = (int) (sample[i] - mean[i]) + 2048;
                            i++;
                            mean[i] += alpha * (sample[i] - mean[i]);
                            int s2 = (int) (sample[i] - mean[i]) + 2048;
                            i++;
                            bOutputStream.write(s1 >> 4);
                            bOutputStream.write( ( (s1 & 0x000f) << 4) + (s2 >> 8));
                            bOutputStream.write(s2 & 0x00ff);
                        }
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                    sampleIndex++;
                    percentage = 100 * sampleIndex / nSamples;
                    if (percentage != oldPercentage) {
                        System.out.println(percentage + "%");
                        oldPercentage = percentage;
                    }
                }
                public void finishProcessing() throws IOException {
                        bOutputStream.close();
                    long t2 = System.currentTimeMillis();
                    System.out.println( (t2 - t1) / 1000.);
                }
            });
            sis.start();
        }
     */


    /**
     * The input is a 512 format file and the output is also 512 format file
     * but contraining only the data from the given electrodes plus the data from
     * electrode 0.
     */
    public static void extractElectrodes(
        String inputFileName, String outputFolder, final int[] electrodes,
        int arrayID) throws IOException {

        String datasetName = new File(outputFolder).getName();
        String outputFileName = outputFolder + File.separator + datasetName + ".bin";
        File f = new File(outputFileName);
        if (f.exists()) {
            throw new IOException("The output file " + outputFileName +
                                  " allready exists.");
        }

        final SynchronizationObject sync = new SynchronizationObject();
        if ( (electrodes.length % 2) != 0) {
            throw new IllegalArgumentException("Even");
        }
        final double t1 = System.currentTimeMillis();
        MultipleCompressedSampleInputStream sis =
            new MultipleCompressedSampleInputStream(inputFileName);
        RawDataHeader512 inputHeader = (RawDataHeader512) sis.getHeader();
        final int nElectrodes = inputHeader.getNumberOfElectrodes();
        final int nSamples = inputHeader.getNumberOfSamples();
        FileOutputStream outputStream = new FileOutputStream(outputFileName);
        final RawDataHeader outputHeader = new RawDataHeader512(
            0,
            electrodes.length + 1,
            inputHeader.getSamplingFrequency(),
            nSamples,
            arrayID,
            RawDataHeader.FORMAT_12BIT_COMPRESSED,
            inputHeader.getDatasetIdentifier(),
            inputHeader.getComment());

        outputStream.write(outputHeader.getBinaryRepresentation());
        final BufferedOutputStream bOutputStream =
            new BufferedOutputStream(outputStream, 1024 * 1024);

        sis.addSampleListener(new SampleListener() {
            long sampleIndex, percentage, oldPercentage;
            byte[] buffer = new byte[outputHeader.getSampleSize()];

            public void processSample(short[] sample) {
                int byteIndex = 0;
                buffer[byteIndex++] = (byte) (sample[0] >> 8);
                buffer[byteIndex++] = (byte) (sample[0] & 0x00ff);

                for (int i = 0; i < electrodes.length; ) {
                    int s1 = sample[electrodes[i]] + 2048;
                    i++;

                    int s2 = sample[electrodes[i]] + 2048;
                    i++;

                    buffer[byteIndex++] = (byte) (s1 >> 4);
                    buffer[byteIndex++] = (byte) ( ( (s1 & 0x000f) << 4) + (s2 >> 8));
                    buffer[byteIndex++] = (byte) (s2 & 0x00ff);
                }

                try {
                    bOutputStream.write(buffer);
                } catch (IOException e) {
                    e.printStackTrace();
                }

                sampleIndex++;
                percentage = 100 * sampleIndex / nSamples;
                if (percentage != oldPercentage) {
                    System.out.println(percentage + "% ");
                    oldPercentage = percentage;
                }
            }


            public void finishSampleProcessing() throws IOException {
                bOutputStream.close();
                long t2 = System.currentTimeMillis();
                System.out.println( (t2 - t1) / 1000.);
                sync.done();
            }
        });

        sync.setWorking();
        sis.start();
        sync.waitUntilDone();
    }


    public static int[] get61ElectrodeSubset(int centerElectrode512) {
        int[] electrodes = new int[64];
        ElectrodeMap m64 = ElectrodeMapFactory.getElectrodeMap(0);
        ElectrodeMap m512 = ElectrodeMapFactory.getElectrodeMap(501);
        Point2D.Double center512 = new Point2D.Double(0, 0);
        m512.getPosition(centerElectrode512, center512);

        Point2D.Double p64 = new Point2D.Double(0, 0);
        Point2D.Double p512 = new Point2D.Double(0, 0);
        for (int i = 1; i <= electrodes.length; i++) {
            if (m64.isDisconnected(i)) {
                electrodes[i - 1] = 0;
                continue;
            }
            m64.getPosition(i, p64);
            for (int j = 0; j < m512.getNumberOfElectrodes(); j++) {
                m512.getPosition(j, p512);
                if (p512.x - center512.x == p64.x && p512.y - center512.y == p64.y) {
                    electrodes[i - 1] = j;
                    break;
                }
            }
        }

        for (int i = 0; i < electrodes.length; i++) {
            System.out.println(i + 1 + " - " + electrodes[i]);
        }

        return electrodes;
    }

}
