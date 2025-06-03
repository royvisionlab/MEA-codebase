/**
 * The simpler decompression used for even number of electrodes cases
 */
__kernel void decompress_samples(__global const uchar* in, __global short* out, int n)
{
	int i = get_global_id(0);
	if (i >= n) return;
	
	int startin = i*3;
	uchar c1 = in[startin];
	uchar c2 = in[startin+1];
	uchar c3 = in[startin+2];

	int startout = i*2;
	out[startout]   = (short)  (c1         << 4) + (c2 >> 4) - 2048;
	out[startout+1] = (short) ((c2 & 0x0f) << 8) +  c3       - 2048;
}


/**
 * Casting @out as an int has two possible advantages: 
 *   1) Runs on older GPUs that don't allow writing for smaller than 32-bit widths
 *   2) May allow improved memory coalescence on output step
 * 
 * I'm not seeing any speed-up on creampuff (X58 + Radeon HD 5770 + OSX 10.6)
 */
__kernel void decompress_samples_int_out(__global const uchar* in, __global int* out, int n)
{
	int i = get_global_id(0);
	if (i >= n) return;
	
	int startin = i*3;
	uchar c1 = in[startin];
	uchar c2 = in[startin+1];
	uchar c3 = in[startin+2];

	ushort s1 = (short)  (c1         << 4) + (c2 >> 4) - 2048;
	ushort s2 = (short) ((c2 & 0x0f) << 8) +  c3       - 2048;
	out[i] = s1 + (s2<<16);
}


/**
 * A more complicated scheme was cooked up for odd numbers of electrodes...
 */
__kernel void decompress_samples_odd(__global const uchar* in, __global short* out, int numelectrodes, int n) 
{
	int i = get_global_id(0);
	if (i >= n)	return;
	
	int ops_per_sample = (numelectrodes-1)/2 + 1;
	int op = i % ops_per_sample;
	int sample = i / ops_per_sample;
	int samplesize = (numelectrodes-1)/2*3 + 2;
		
	if (op == 0) {
		int startin  = sample*samplesize;
		int startout = sample*numelectrodes;
		uchar c1 = in[startin];
		uchar c2 = in[startin+1];
		out[startout] = (short) (c1 << 8) + c2;
	} else {
		int startin  = sample*samplesize    + 2 + (op-1)*3;
		int startout = sample*numelectrodes + 1 + (op-1)*2;

		uchar c1 = in[startin];
		uchar c2 = in[startin+1];
		uchar c3 = in[startin+2];
		
		out[startout]   = (short)  (c1         << 4) + (c2 >> 4) - 2048;
		out[startout+1] = (short) ((c2 & 0x0f) << 8) +  c3       - 2048;
	}
}

// Had to use uchars otherwise C was trying to keep the sign under the shift operation