package edu.ucsc.neurobiology.vision.convert;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.anf.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.snf.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class ConvertRawDataToElectrodeMajor
    extends AbstractCalculation implements SampleListener {

    private String rawDataFileName, outputFolder, sigmasFileName;
    private Vision app;
    private int nSamples;
    private long sampleIndex = 0;
    private long percentage, oldPercentage = -1;
    private int index = 0, writeNumber = 0;
    private short[][] buffer;
    private byte[] byteBuffer;
    private int nElectrodes;
    private long headerLength;
    private final int bufferSize = 512 * 1024;
    private RandomAccessFile outputStream;
    private String outputFileName;
    private SynchronizationObject sync = new SynchronizationObject();


    public void startCalculation() throws Exception {
//        if (bufferSize % 2 != 0) {
//            throw new Exception("bufferSize should be a multiple of 2");
//        }

        app = Vision.getInstance();
        app.sendMessage("Convert Data to Electrode Major...");
        app.startProgressBar();
        String datasetName =
            StringUtil.removeExtension(new File(outputFolder).getName());
        outputFileName = outputFolder + File.separator + datasetName +
                         VisionParams.REM_FILE_EXTENSION_512;


        sync.setWorking();

        MultipleCompressedSampleInputStream sis =
            new MultipleCompressedSampleInputStream(rawDataFileName);
        RawDataHeader inputHeader = sis.getHeader();
        nElectrodes = inputHeader.getNumberOfElectrodes();
        nSamples = inputHeader.getNumberOfSamples();
        outputStream = new RandomAccessFile(outputFileName, "rw");
        inputHeader.setFormat(RawDataHeader.FORMAT_8BIT_UNCOMPRESSED);
        final byte[] binHeader = inputHeader.getBinaryRepresentation();
        outputStream.write(binHeader);
        headerLength = binHeader.length;

        buffer = new short[nElectrodes][bufferSize];
        byteBuffer = new byte[3 * bufferSize / 2];

        sis.addSampleListener(this);
        sis.start();
        sync.waitUntilDone();

        buffer = null;
        byteBuffer = null;

        calculateSpikeAmpHistograms(outputFileName, sigmasFileName);

        app.sendMessage("Convert To Electrode Major: Done.");
        app.endProgressBar();
        Vision.getInstance().getCalculationManager().calculationDone();
    }


    private static void calculateSpikeAmpHistograms(
        String outputFileName, String sigmasFileName) throws IOException {

        // calculate spike amplitude histograms
        ElectrodeMajorRawDataFile rawFile = new ElectrodeMajorRawDataFile(outputFileName);
        final int nElectrodes = rawFile.getHeader().getNumberOfElectrodes();
        float[] sigmas = SpikeFinding.getSigmas(sigmasFileName, nElectrodes);
        final int nSamples = rawFile.getHeader().getNumberOfSamples();
        int[][] hist = new int[nElectrodes][2048 + 2];
        short[] data = new short[nSamples];
        double[] spikeT = new double[10 * 1024 * 1024];
        float[] spikeA = new float[10 * 1024 * 1024];
        for (int electrode = 1; electrode < nElectrodes; electrode++) {
            rawFile.readRawData(electrode, data, nSamples);
            int n = SpikeFinder.findSpikes(
                data, 2 * sigmas[electrode], spikeT, spikeA, 0, 0);
  //          System.out.println("Spike Finding On " + electrode + ", " + n);
            for (int i = 0; i < n; i++) {
                hist[electrode][ (int) spikeA[i]]++;
            }
        }
        rawFile.close();

        SerialNeuronPanel.saveAmplitudeHistograms(
            hist, StringUtil.removeExtension(outputFileName) + ".sah");
    }


    private void write() throws IOException {
        for (int electrode = 0; electrode < nElectrodes; electrode++) {
            outputStream.seek(
                headerLength +
                3 * ( (long) electrode * nSamples + writeNumber * bufferSize) / 2);

            for (int i = 0, byteIndex = 0; i < bufferSize; ) {
                int s1 = buffer[electrode][i++] + 2048;
                int s2 = buffer[electrode][i++] + 2048;

                byteBuffer[byteIndex++] = (byte) (s1 >> 4);
                byteBuffer[byteIndex++] =
                    (byte) ( ( (s1 & 0x000f) << 4) + (s2 >> 8));
                byteBuffer[byteIndex++] = (byte) (s2 & 0x00ff);
            }

            outputStream.write(byteBuffer, 0, 3 * index / 2);
        }
    }


    public void processSample(short[] sample) {
        for (int electrode = 0; electrode < nElectrodes; electrode++) {
            if (sample[electrode] < -2048) {
                sample[electrode] = -2048;
            }
            if (sample[electrode] > 2047) {
                sample[electrode] = 2047;
            }
            buffer[electrode][index] = sample[electrode];
        }

        index++;
        if (index == bufferSize) {
            try {
                write();
            } catch (IOException e) {
                e.printStackTrace();
            }
            index = 0;
            writeNumber++;
        }

        sampleIndex++;
        percentage = 100L * sampleIndex / nSamples;
        if (percentage != oldPercentage) {
            app.setProgress((int) percentage);
        //    System.out.println(percentage + "%");
            oldPercentage = percentage;
        }
    }


    public void finishSampleProcessing() throws IOException {
        write();
        outputStream.close();

        sync.done();
    }


    public void setParameters(HashMap<String, String> parameters) {
        rawDataFileName = (String) parameters.get("Raw_Data_File");
        sigmasFileName = (String) parameters.get("Sigmas_File");
        outputFolder = (String) parameters.get("Save_To_File");
    }
}
