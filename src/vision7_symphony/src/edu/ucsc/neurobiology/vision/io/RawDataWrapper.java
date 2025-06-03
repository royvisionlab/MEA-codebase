package edu.ucsc.neurobiology.vision.io;

import java.io.IOException;



/**
 * Allows access to raw data, based on a a string that follows the rules of io.DataFileStringParser.
 * This function is recommended for when the data to read is small, or not much time is available 
 * for development.  io.MultipleCompressedSampleInputStream is faster otherwise.  
 * io.MultipleCompressedSampleInputStream requires you to read the data sequentially.
 * 
 * @author Matthew Grivich, The Salk Institute
 */

public class RawDataWrapper  {

    int count = 0;
    short[][] dataBuffer;
    boolean finished = false;
    long samplesToRead;

    private long[] startSamples, stopSamples;

    private RawDataHeader[] headers;
    private RawDataFile[] files;

    private RawDataHeader totalHeader;

    public RawDataWrapper(String rawDataSources) throws IOException {

        DataFileStringParser parser = new DataFileStringParser(rawDataSources);

        String[] rawDataSourcesSplit = parser.getDatasets();
        double[] startTimes = parser.getStartTimes();
        double[] stopTimes = parser.getStopTimes();

        headers = new RawDataHeader[rawDataSourcesSplit.length];
        startSamples = new long[rawDataSourcesSplit.length];
        stopSamples = new long[rawDataSourcesSplit.length];
        files = new RawDataFile[rawDataSourcesSplit.length];

        for (int i = 0; i < rawDataSourcesSplit.length; i++) {
            files[i] = new RawDataFile(rawDataSourcesSplit[i]);
            headers[i] = files[i].getHeader();

            startSamples[i] = ((long) startTimes[i]) * headers[i].getSamplingFrequency();

            stopSamples[i] = (long) stopTimes[i];
            if(stopSamples[i] == Long.MAX_VALUE) {
                stopSamples[i] = headers[i].getNumberOfSamples();	
            } else {
                stopSamples[i] = ((long) stopTimes[i]) * headers[i].getSamplingFrequency();
            }


        }

        int nSamples = 0;
        for (int i = 0; i < rawDataSourcesSplit.length; i++) {
            long nCurSamples = headers[i].getNumberOfSamples();
            //samples before start are zeroed, not removed.  samples after stop are removed.
            if(stopSamples[i] < nCurSamples) { 
                nCurSamples = stopSamples[i];
            }
            nSamples += nCurSamples;
        }
        totalHeader = new RawDataHeader512(
                0,
                headers[0].getNumberOfElectrodes(),
                headers[0].getSamplingFrequency(),
                nSamples,
                headers[0].getArrayID(),
                headers[0].getFormat(),
                headers[0].getDatasetIdentifier(), 
        "");


    }

    public RawDataHeader getHeader() {
        return totalHeader;
    }


    /**
     * 
     * @param data short[][] The data array to be filled [sample][nElectrodes]
     * 
     */
    public void getData(long startSample, short[][] data) throws IOException  {

        if(startSample + data.length > totalHeader.getNumberOfSamples()) {
            throw new IOException("Requested raw data is after the end of the file(s).");
        }

        if(startSample < 0) {
            throw new IOException("Seek index is negative.");
        }

        long currentStartSample = startSample;
        int writingSample = 0; //keeps track of the output position.

        while(writingSample < data.length) {

            int file=0;

            //determines which file is next, and where to start reading it.
            while(currentStartSample >= stopSamples[file]) {
                currentStartSample -= stopSamples[file];
                file++;
            }

            //either read till the end of the requested data or till the end of the current file.
            int currentSamplesToRead = (int) Math.min(data.length - writingSample, stopSamples[file] - currentStartSample);
            short[][] currentData = files[file].getData(currentStartSample, currentSamplesToRead);

            for(int j=0; j<currentData.length; j++) {

                //set samples before startSamples[file] to zero
                for(int electrode=0; electrode<currentData[0].length; electrode++) {
                    if(currentStartSample + j < startSamples[file] && electrode > 0) {
                        data[writingSample][electrode] = 0; 
                    } else {
                        data[writingSample][electrode] = currentData[j][electrode]; 
                    }
                }

                writingSample++;
            }

            //set current start sample to the overall start sample
            currentStartSample = startSample + writingSample;


        }
    }


    public short[][] getData(long startSample, int nSamples, int nElectrodes) throws IOException  {
        short[][] data = new short[nSamples][nElectrodes];
        getData(startSample, data);
        return data;
    }



    public static void main(String args[]) {
        try {
            RawDataWrapper rdw = new RawDataWrapper("/data/raw/2000-12-14-1/data051.bin(1-2),data051.bin(3-4)");

            System.out.println(rdw.getHeader().getNumberOfSamples());
            short[][] data = new short[20000][65];
            rdw.getData(29998, data);
            System.out.println(data[2][1]);


        } catch(IOException ex) {
            ex.printStackTrace();
        }

    }
}




//package edu.ucsc.neurobiology.vision.io;

//import java.io.IOException;

//public class RawDataWrapper implements SampleListener {

//int count = 0;
//short[][] dataBuffer;
//boolean finished = false;
//long samplesToRead;



//MultipleCompressedSampleInputStream ms;

//public RawDataHeader getHeader(String rawDataSources) throws IOException {
//ms = new MultipleCompressedSampleInputStream(rawDataSources, 500*770, 100, 0, 50);
//return ms.getHeader();
//}



//public short[][] getRawData(String rawDataSources, long relativeStart, short[][] dataBuffer) throws IOException {


//this.dataBuffer = dataBuffer;
//this.samplesToRead = dataBuffer.length;
//count = 0;
//ms = new MultipleCompressedSampleInputStream(rawDataSources,
//500*770, 100, relativeStart, samplesToRead);

//ms.addSampleListener(this);

//ms.start();

//while(!ms.isFinished()) {
////try {
////Thread.sleep(2);
////} catch(Exception e) {
////e.printStackTrace();
////}
//}

////System.out.println("finished: "  + ms.isFinished());
////System.out.println("count:" + count);


//return dataBuffer;


//}

//public void processSample(short[] sample) {
////SAMPLE GETS MODIFIED, NEED TO COPY
//if(count < samplesToRead) {
//for(int i=0; i<sample.length; i++) {
//dataBuffer[count][i] = sample[i];
//}			
//count++;
//} else {

//finishSampleProcessing();


//}

//}

//public void finishSampleProcessing()  {
////System.out.println("finished");
////System.out.println("count: " + count);
//finished = true;
////System.out.println(ms.isFinished());

//}



//public static void main(String args[]) {
//RawDataWrapper rdw = new RawDataWrapper();
//try {
//rdw.getRawData("/data/2000-12-14-1vision/data051.bin", 0, new short[100][65]);
//} catch(IOException ex) {
//ex.printStackTrace();
//}
//}
//}
