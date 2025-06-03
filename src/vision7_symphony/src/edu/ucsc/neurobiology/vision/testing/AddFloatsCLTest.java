package edu.ucsc.neurobiology.vision.testing;

import java.io.IOException;

import com.nativelibs4java.opencl.*;
import com.nativelibs4java.opencl.CLMem.Usage;
import com.nativelibs4java.util.IOUtils;

import org.bridj.Pointer;
import static org.bridj.Pointer.*;
import static java.lang.Math.*;


public class AddFloatsCLTest {

    public static void main(String[] args) throws IOException {
        CLContext context = JavaCL.createBestContext(CLPlatform.DeviceFeature.GPU);
        CLQueue queue = context.createDefaultQueue();

        int n = 1024;
        Pointer<Float> aPtr = allocateFloats(n), bPtr = allocateFloats(n);

        for (int i = 0; i < n; i++) {
            aPtr.set(i, (float)cos(i));
            bPtr.set(i, (float)sin(i));
        }

        // Create OpenCL input and output buffers
        CLBuffer<Float> a = context.createBuffer(Usage.Input, aPtr);
        CLBuffer<Float> b = context.createBuffer(Usage.Input, bPtr);
        CLBuffer<Float> out = context.createFloatBuffer(Usage.Output, n);

        // Read the program sources and compile them :
        String src = IOUtils.readText(AddFloatsCLTest.class.getResource("add_floats.cl"));
        CLProgram program = context.createProgram(src);

        // Get and call the kernel :
        CLKernel addFloatsKernel = program.createKernel("add_floats");
        addFloatsKernel.setArgs(a, b, out, n);
        CLEvent addEvt = addFloatsKernel.enqueueNDRange(queue, new int[] { n });

        Pointer<Float> outPtr = out.read(queue, addEvt); // blocks until add_floats finished

        // Print the first 10 output values :
        for (int i = 0; i < 10 && i < n; i++)
            System.out.println("out[" + i + "] = " + outPtr.get(i));		
    }
}