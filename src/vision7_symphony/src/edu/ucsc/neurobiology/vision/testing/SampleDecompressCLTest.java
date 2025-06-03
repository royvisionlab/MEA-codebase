package edu.ucsc.neurobiology.vision.testing;

import java.io.IOException;
import java.util.ListIterator;

import org.bridj.Pointer;

import com.nativelibs4java.opencl.CLContext;
import com.nativelibs4java.opencl.CLPlatform;
import com.nativelibs4java.opencl.JavaCL;

import edu.ucsc.neurobiology.vision.io.AsyncBufferCL;
import edu.ucsc.neurobiology.vision.io.BlockSampleDecompressorCL;
import edu.ucsc.neurobiology.vision.io.RawDataFile;
import edu.ucsc.neurobiology.vision.io.SampleBlockQueueCL;

public class SampleDecompressCLTest {

    private static long startTime;
    private static String filename;
    private static CLContext context;
    private static BlockSampleDecompressorCL decompressor;
    private static RawDataFile rdf;
    
    
public static void sampleBlockQueueTest() throws IOException, InterruptedException {
        int samples = 1000;
        int nBuffers = 5;

        SampleBlockQueueCL sampleQueue = new SampleBlockQueueCL(decompressor, nBuffers, nBuffers, samples);
        
        sampleQueue.start();		
        int numResults = 0;
//		int i = 0;
        while (true) {
            AsyncBufferCL<Short> sampleBuffer = sampleQueue.take();
            if (sampleBuffer.getBuffer() == null) break;

//			Pointer<Short> ptr = sampleBuffer.getBuffer().read(context.createDefaultQueue(), sampleBuffer.getEvts());
            sampleBuffer.close();
            numResults += samples;
            
//			rdf.getData(i++*samples, data);
//
//			ListIterator<Short> it = ptr.iterator();
//			for (int sample = 0; sample < 1; sample++) {
//				for (int electrode = 0; electrode < nElectrodes; electrode++) {
//					short gpu = it.next();
//					if (gpu != data[sample][electrode])
//						System.out.println("s" + sample + "e" + electrode + "    " + gpu + "    " + data[sample][electrode]);
//				}
//			}

        }
        System.out.println("Samples read: " + numResults);
        System.out.println("Total time: " + ((System.currentTimeMillis()-startTime) / 1000) + " s");
        
//		// Debug; compare to output from RawDataFile
//		RawDataFile rdf = new RawDataFile(args[0]);
//		int nElectrodes = rdf.getHeader().getNumberOfElectrodes();
//		short[][] data = new short[samples][nElectrodes];
//		for (int i = 0; i < 20; i++) {
//			AsyncBufferCL<Short> sampleBuffer = sampleQueue.take();
//			Pointer<Short> ptr = sampleBuffer.getBuffer().read(context.createDefaultQueue(), sampleBuffer.getEvts());
//			sampleBuffer.close();
//
//			startTime = System.currentTimeMillis();
//			rdf.getData(i*samples, data);
//			System.out.println("Finished CPU decompression: " + samples + " samples in " + (System.currentTimeMillis()-startTime) + " ms");		
//
//			ListIterator<Short> it = ptr.iterator();
//			for (int sample = 0; sample < 1; sample++) {
//				for (int electrode = 0; electrode < nElectrodes; electrode++) {
//					short gpu = it.next();
//					if (gpu != data[sample][electrode])
//						System.out.println("s" + sample + "e" + electrode + "    " + gpu + "    " + data[sample][electrode]);
//				}
//			}
//		}
//		
//		sampleQueue.close();
    }
    
    
    public static void blockSampleDecompressorTest(int samples) throws IOException {
        startTime = System.currentTimeMillis();
        AsyncBufferCL<Short> sampleBuffer = decompressor.decompress(samples);
        decompressor.close();
        sampleBuffer.waitFor();
        System.out.println("Complete GPU decompress including wait for events: " + samples + " samples in " + (System.currentTimeMillis()-startTime) + " ms");
        
        startTime = System.currentTimeMillis();
        Pointer<Short> ptr = sampleBuffer.getBuffer().read(context.createDefaultQueue(), sampleBuffer.getEvts());
        sampleBuffer.close();
        System.out.println("Retrieved output: " + (System.currentTimeMillis()-startTime) + " ms");
        
        // Debug; compare to output from RawDataFile
        startTime = System.currentTimeMillis();
        rdf = new RawDataFile(filename);
        int nElectrodes = rdf.getHeader().getNumberOfElectrodes();
        short[][] data = new short[samples][nElectrodes];
        rdf.getData(0, data);
        System.out.println("Finished CPU decompression: " + samples + " samples in " + (System.currentTimeMillis()-startTime) + " ms");		

        ListIterator<Short> it = ptr.iterator();
        for (int sample = 0; sample < samples; sample++) {
            for (int electrode = 0; electrode < nElectrodes; electrode++) {
                short gpu = it.next();
                if (gpu != data[sample][electrode])
                    System.out.println("s" + sample + "e" + electrode + "    " + gpu + "    " + data[sample][electrode]);
            }
        }
    }
    
    
    public static void main(String[] args) throws IOException, InterruptedException {
        filename = args[0];
        int samples = 1000;
        
        startTime = System.currentTimeMillis();
        context = JavaCL.createBestContext(CLPlatform.DeviceFeature.GPU);
        System.out.println("Created context: " + (System.currentTimeMillis()-startTime) + " ms");

        startTime = System.currentTimeMillis();
        decompressor = BlockSampleDecompressorCL.create(context, context.createDefaultQueue(), args[0]);
        System.out.println("Created decompressor: " + (System.currentTimeMillis()-startTime) + " ms");
        
//		blockSampleDecompressorTest(samples);
        sampleBlockQueueTest();
    }
    
}
