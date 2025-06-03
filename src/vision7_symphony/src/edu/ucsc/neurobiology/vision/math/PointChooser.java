package edu.ucsc.neurobiology.vision.math;

import edu.ucsc.neurobiology.vision.math.RandomMersenne;

/**
 * Class for providing a mapping between a given array, and a shorter 
 * array with points randomly removed.
 * 
 * Usage:
 * double[] array = new array[nTotalPoints];
 * PointChooser pc = new PointChooser(nChosenPoints, nTotalPoints);
 * 
 * for(int i=0; i<pc.getNChosen; i++) {
 * 	System.out.println(array[pc.get(i)]);
 * }
 * 
 * 
 * 
 * @author Matthew Grivich, The Salk Institute
 *
 */

public class PointChooser {
    private int nChosenPoints;
    private int nTotalPoints;

    private int[] redirect;
    
    public PointChooser(int nChosenPoints, int nTotalPoints) {
        if (nChosenPoints < 0 || nTotalPoints < 0)
            throw new IllegalArgumentException("Requested number points cannot be negative.");		
        
        // if choose all, set up the trivial solution
        if (nChosenPoints >= nTotalPoints) {
            nChosenPoints = nTotalPoints;
            redirect = new int[nTotalPoints];
            for(int i=0; i<nTotalPoints; i++) {
                redirect[i] = i;
            }
            this.nChosenPoints = nChosenPoints;
            this.nTotalPoints = nTotalPoints;
            return;
        }
        
        RandomMersenne rng = new RandomMersenne(RandomMersenne.DEFAULT_SEED); 
        
        // find number of bits to pull from rng.
        int bits = 0;
        while(nTotalPoints > (1 << bits)) bits++;
        
        this.nChosenPoints = nChosenPoints;
        this.nTotalPoints = nTotalPoints;
        boolean[] chosen = new boolean[nTotalPoints]; // =false
        boolean trueMeansKeep;
        int nToFlip;
        
        //require that no more that 50% of the points are flipped
        if (nChosenPoints <= nTotalPoints - nChosenPoints) {
            nToFlip = nChosenPoints;
            trueMeansKeep = true;
        } else {
            nToFlip = nTotalPoints - nChosenPoints;
            trueMeansKeep = false;
        }
        nToFlip  = (nChosenPoints <= nTotalPoints - nChosenPoints) ? nChosenPoints : nTotalPoints - nChosenPoints;
        for (int i = 0; i < nToFlip; i++) {
            while(true) {
                //use most significant bits
                int candidate = (rng.nextInt() >>> (32-bits));
                //don't use the next number if it is out of the acceptable range
                if (candidate >= nTotalPoints) continue;
                
                //if point has not already been flipped, flip it.  Otherwise, reroll.
                //other algorithms lead to bias or efficiency problems
                //algorithm is guaranteed to converge, because never are more than 50% of the
                //points flipped.
                
                if (!chosen[candidate])  {
                    chosen[candidate] = true;
                    break;
                }
            }
        }
        
        redirect = new int[nChosenPoints];
        int index = 0;
        for (int i = 0; i < nChosenPoints; i++) {
            while (chosen[index++] != trueMeansKeep); //while not a keeper
            //if a keeper, add it to the redirect list.
            redirect[i] = index - 1;
        }
        

        
    }
    /**
     * Given an index in the reduced length array, the index in the original array is returned.
     * 
     * 
     * @param point next point in reduced length array
     * @return index  point in original array
     */
    public int get(int point) {
        return redirect[point];
    }	
    
    public int getNChosen() {
        return nChosenPoints;
    }
    
    public int getNTotal() {
        return nTotalPoints;
    }
    
    public static void main(String args[]) {
        PointChooser pc = new PointChooser(9, 16);
        for(int i=0; i<pc.getNChosen(); i++) {
            System.out.println(pc.get(i));
        }
        
    }

}
