package edu.ucsc.neurobiology.vision.anf;

import java.io.*;
import java.util.*;

import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.calculations.*;
import edu.ucsc.neurobiology.vision.io.*;
import edu.ucsc.neurobiology.vision.math.*;
import edu.ucsc.neurobiology.vision.util.*;


/**
 * Calculation to filter the raw data using a high pass FIR filter:
 * HP FIR, Least-Squares, order 100, fStop 60, fPass 230, wStop 1, wPass 2
 * -3dB point is at 200 Hz, DC supression -70 dB, passband ripple -0.3:0.4 dB.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public class RawDataFiltering
    extends AbstractCalculation {

    // HP FIR, Least-Squares, order 100, fStop 60, fPass 230, wStop 1, wPass 2
    // -3dB point is at 200 Hz, DC supression -70 dB, passband ripple -0.3:0.4 dB.
    private final static double[] highPassFilter = {
        -0.002626261104497,
        -0.002885383430372,
        -0.003153123468302,
        -0.003429132199827,
        -0.003713029714947,
        -0.004004405924799,
        -0.004302821386413,
        -0.004607808237606,
        -0.004918871239608,
        -0.005235488924628,
        -0.005557114845115,
        -0.005883178921084,
        -0.006213088881465,
        -0.006546231795059,
        -0.006881975686323,
        -0.007219671230831,
        -0.00755865352497,
        -0.007898243924055,
        -0.008237751942834,
        -0.008576477212001,
        -0.00891371148416,
        -0.009248740682421,
        -0.009580846984623,
        -0.009909310936025,
        -0.01023341358313,
        -0.01055243862126,
        -0.01086567454828,
        -0.01117241681702,
        -0.01147196997871,
        -0.01176364980987,
        -0.01204678541506,
        -0.01232072129806,
        -0.01258481939396,
        -0.012838461055,
        -0.01308104898285,
        -0.01331200910053,
        -0.01353079235708,
        -0.01373687645847,
        -0.01392976751855,
        -0.01410900162391,
        -0.01427414630698,
        -0.01442480192214,
        -0.0145606029196,
        -0.0146812190125,
        -0.01478635623304,
        -0.01487575787357,
        -0.01494920530937,
        -0.01500651870003,
        -0.01504755756682,
        -0.01507222124403,
        0.9849195507976,
        -0.01507222124403,
        -0.01504755756682,
        -0.01500651870003,
        -0.01494920530937,
        -0.01487575787357,
        -0.01478635623304,
        -0.0146812190125,
        -0.0145606029196,
        -0.01442480192214,
        -0.01427414630698,
        -0.01410900162391,
        -0.01392976751855,
        -0.01373687645847,
        -0.01353079235708,
        -0.01331200910053,
        -0.01308104898285,
        -0.012838461055,
        -0.01258481939396,
        -0.01232072129806,
        -0.01204678541506,
        -0.01176364980987,
        -0.01147196997871,
        -0.01117241681702,
        -0.01086567454828,
        -0.01055243862126,
        -0.01023341358313,
        -0.009909310936025,
        -0.009580846984623,
        -0.009248740682421,
        -0.00891371148416,
        -0.008576477212001,
        -0.008237751942834,
        -0.007898243924055,
        -0.00755865352497,
        -0.007219671230831,
        -0.006881975686323,
        -0.006546231795059,
        -0.006213088881465,
        -0.005883178921084,
        -0.005557114845115,
        -0.005235488924628,
        -0.004918871239608,
        -0.004607808237606,
        -0.004302821386413,
        -0.004004405924799,
        -0.003713029714947,
        -0.003429132199827,
        -0.003153123468302,
        -0.002885383430372,
        -0.002626261104497
    };


    private String rawDataFileName, outputFolder;
    private Vision app;
    int nSamples;


    public void startCalculation() throws IOException {
        app = Vision.getInstance();
       
        
      File folder = new File(outputFolder);
        folder.mkdirs();
        String datasetName = StringUtil.removeExtension(folder.getName());
        String saveToFile = outputFolder + File.separator + datasetName +
                            VisionParams.BIN_FILE_EXTENSION_512;
   //     System.out.println("Save to: " + saveToFile);

        double[] taps = highPassFilter;

        app.sendMessage("Raw Data Filtering ...");
        app.startProgressBar();
        filter(rawDataFileName, saveToFile, taps);
        app.endProgressBar();

        Vision.getInstance().getCalculationManager().calculationDone();
    }


    public void setParameters(HashMap<String, String> parameters) {
        rawDataFileName = (String) parameters.get("Raw_Data_File");
        outputFolder = (String) parameters.get("Output Folder");
    }


    public void filter(
        String inputFileName, String outputFileName, double[] taps) throws IOException {
 
        final SynchronizationObject sync = new SynchronizationObject();
        sync.setWorking();

        final MultipleCompressedSampleInputStream sis =
            new MultipleCompressedSampleInputStream(inputFileName);
        RawDataHeader inputHeader = sis.getHeader();

        final int nElectrodes = inputHeader.getNumberOfElectrodes();
        final int nSamples = inputHeader.getNumberOfSamples();

        File f = new File(outputFileName);
        if (f.exists() && f.isFile()) {
            f.delete();
        }
        final FileOutputStream outputStream = new FileOutputStream(outputFileName);
        final BufferedOutputStream bOutputStream =
            new BufferedOutputStream(outputStream, 100 * 1024 * 1024);
//        final BufferedOutputStream bOutputStream =
//            new BufferedOutputStream(outputStream, 1 * 1024 * 1024);
        bOutputStream.write(inputHeader.getBinaryRepresentation());

        final FIRFilter[] filter = new FIRFilter[nElectrodes];
        for (int i = 0; i < nElectrodes; i++) {
            filter[i] = new FIRFilter(taps);
        }

        // no TTL filtering
//                buffer[byteIndex++] = (byte) (sample[0] >> 8);
//                buffer[byteIndex++] = (byte) (sample[0] & 0x00ff);


        sis.addSampleListener(new SampleListener() {
            long sampleIndex = 0;
            long percentage, oldPercentage = -1;
            byte[] buffer = new byte[2 + (nElectrodes - 1) * 3 / 2];
            int byteIndex;

            public void processSample(short[] sample) {
                sampleIndex++;
                byteIndex = 0;
                
                
                // TTL filtering, keeps timing synchronized
                int s = (int) filter[0].filter(sample[0]);
                buffer[byteIndex++] = (byte) (s >> 8);
                buffer[byteIndex++] = (byte) (s & 0x00ff);

                for (int i = 1; i < nElectrodes; ) {
                    int s1 = (int) filter[i].filter(sample[i]) + 2048;
                    i++;
                    int s2 = (int) filter[i].filter(sample[i]) + 2048;
                    i++;

                    // zero the artifact
//                    if (sampleIndex <= 50) {
//                        s1 = 2048;
//                        s2 = 2048;
//                    }

                    buffer[byteIndex++] = (byte) (s1 >> 4);
                    buffer[byteIndex++] = (byte) ( ( (s1 & 0x000f) << 4) + (s2 >> 8));
                    buffer[byteIndex++] = (byte) (s2 & 0x00ff);
                }

                try {
                    //remove offset
                    if(sampleIndex > 50) {
                        bOutputStream.write(buffer);
                    }
                } catch (IOException e) {
                    e.printStackTrace();
                }

                percentage = 100 * sampleIndex / nSamples;
                if (percentage != oldPercentage) {
                    app.setProgress((int) percentage);
                    //System.out.println(percentage + "%");
                    //oldPercentage = percentage;
                }
            }


            public void finishSampleProcessing() throws IOException {

                bOutputStream.close();
                outputStream.close();
                sync.done();
            }
        });

        sis.start();
        sync.waitUntilDone();
    }
}
