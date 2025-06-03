package edu.ucsc.neurobiology.vision.util;

import java.util.List;

public final class ParallelUtil {

    public static class StartEndIndices {
        public int start;
        public int end;
        
        StartEndIndices(int start, int end) {
            this.start = start;
            this.end = end;
        }
    }
    
    public static StartEndIndices subIndex(int thisPart, int numParts, int startIndex, int endIndex) {
        float[] proportions = new float[numParts];
        for (int i = 0; i < numParts; i++) proportions[i] = 1.0f / ((float) numParts);
        return subIndex(thisPart, numParts, startIndex, endIndex, proportions);
    }

    
    public static StartEndIndices subIndex(int thisPart, int numParts, int numIndices) {
        return subIndex(thisPart, numParts, 0, numIndices - 1);
    }
    
    public static StartEndIndices subIndex(int thisPart, int numParts, int numIndices, float[] proportions) {
        return subIndex(thisPart, numParts, 0, numIndices - 1, proportions);
    }

    
    public static StartEndIndices subIndex(int thisPart, int numParts, int startIndex, int endIndex, float[] proportions) {
        float startProportion = 0;
        for (int i = 0; i < thisPart; i++) startProportion += proportions[i];
        float endProportion = startProportion + proportions[thisPart];

        int numIndices = endIndex - startIndex + 1;
        int thisSubStartIndex = Math.round(numIndices * startProportion) + startIndex;
        int thisSubEndIndex   = Math.round(numIndices * endProportion)   + startIndex - 1;
        
        if (thisPart == (numParts - 1)) thisSubEndIndex = endIndex;
        return new StartEndIndices(thisSubStartIndex, thisSubEndIndex);
    }
    

    /**
     * Was having trouble hacking this into the generic version below, so just recode for now :(
     * Return is not cloned; you must clone it yourself if you want to modify without affecting original.
     */
    private static int[] getMyPart(int[] full, StartEndIndices myInd) {
        int myLength = myInd.end - myInd.start + 1;
        int[] myPart = new int[myLength];
        for (int i = 0; i < myLength; i++) myPart[i] = full[i+myInd.start];
        return myPart;
    }

    public static int[] getMyPart(int thisPart, int numParts, int[] full) {
        StartEndIndices myInd = subIndex(thisPart, numParts, full.length);
        return getMyPart(full, myInd);
    }
    
    public static int[] getMyPart(int thisPart, int numParts, int[] full, float[] proportions) {
        StartEndIndices myInd = subIndex(thisPart, numParts, full.length, proportions);
        return getMyPart(full, myInd);
    }	

    
    /**
     * The meta version...
     */
    @SuppressWarnings("unchecked")
    private static <T> List<T> getMyPart(List<T> full, StartEndIndices myInd) throws Exception {
        Class<? extends List> fullClass = full.getClass();
        List<T> myPart = fullClass.getConstructor().newInstance();		
        for (int i = myInd.start; i <= myInd.end; i++) myPart.add(full.get(i));
        return myPart;
    }
    
    public static <T> List<T> getMyPart(int thisPart, int numParts, List<T> full) throws Exception {
        StartEndIndices myInd = subIndex(thisPart, numParts, full.size());
        return getMyPart(full, myInd);
    }
    
    public static <T> List<T> getMyPart(int thisPart, int numParts, List<T> full, float[] proportions) throws Exception {
        StartEndIndices myInd = subIndex(thisPart, numParts, full.size(), proportions);
        return getMyPart(full, myInd);
    }
    
    
}