package edu.ucsc.neurobiology.vision.test;


import java.util.regex.*;
import java.io.*;
import java.util.*;


public class RegexTestHarness {

    /**
     * Examples legal values:
     * /folder/data000, data003, data005,  concatenate data000, data003, and data005
     * /data000 - data003,  concatenate data000, data001, data002, and data003
     * /data000.bin, data003.bin, concatenate data000.bin and data003.bin (individual files, rather than folders of files).
     * /data000(5-10), data000 from second 5  - second 10
     * /folder/data000(5-) - data001(-10), from data000 at second 5 to data001 at second 10.
     * 
     * /data001.bin(2-) - data003(-1), data005
     * 
     * White space ignored.
     * 
     * 
     * 
     */




    public static void main(String[] args) throws PatternSyntaxException {


        ArrayList<String> datasList = new ArrayList<String>();
        ArrayList<Double> startList = new ArrayList<Double>();
        ArrayList<Double> stopList = new ArrayList<Double>();

        String regex = File.separator + "snld" + File.separator + "data005(5-19), data002(5-)-data005(-26)" + File.separator;


        //remove white space	
        regex = regex.replaceAll("[\\s]", "");

        //remove trailing file separator
        if(regex.endsWith(File.separator)) {
            regex.substring(0, regex.length()-1);
        }

        regex = regex.substring(0, regex.length() -1);
        int cutPoint = regex.lastIndexOf(File.separator);
        String filePath = regex.substring(0, cutPoint+1);
        regex = regex.substring(cutPoint+1);



        String[] rawDataSources = regex.split(",");

        String[] firstDatas = new String[rawDataSources.length];
        String[] secondDatas = new String[rawDataSources.length];
        String[] firstRanges = new String[rawDataSources.length];
        String[] secondRanges = new String[rawDataSources.length];
        double[] startSeconds = new double[rawDataSources.length];
        double[] stopSeconds = new double[rawDataSources.length];


        for(int i=0; i<rawDataSources.length; i++) {
            startSeconds[i] = 0;
            stopSeconds[i] = Double.POSITIVE_INFINITY;
            firstDatas[i] = "";
            secondDatas[i] = "";
            firstRanges[i] = "";
            secondRanges[i] = "";
        }


        Pattern pattern = null;
        Matcher matcher = null;

        //For each group in comma separated list.
        for(int i=0; i<rawDataSources.length; i++) {
            boolean[] parsed = new boolean[rawDataSources[i].length()];


            //Find anything of the form data*   These are folders.  
            //Will find data* before .bin as well.  Overridden in next step.
            pattern = Pattern.compile("data[\\d][\\d][\\d]");
            matcher = pattern.matcher(rawDataSources[i]);

            if(matcher.find()) {
                firstDatas[i] =  matcher.group();
                for(int j=matcher.start(); j< matcher.end(); j++) parsed[j] = true;	
            }
            if(matcher.find()) {
                secondDatas[i] = matcher.group();
                for(int j=matcher.start(); j< matcher.end(); j++) parsed[j] = true;	
            }


            //find anything of the form dataXXX.bin   These are files.
            pattern = Pattern.compile("data[\\d][\\d][\\d]\\.bin");
            matcher = pattern.matcher(rawDataSources[i]);

            if(matcher.find()) {
                firstDatas[i] = matcher.group();
                for(int j=matcher.start(); j< matcher.end(); j++) parsed[j] = true;	
            }
            if(matcher.find()) {
                secondDatas[i] = matcher.group();
                for(int j=matcher.start(); j< matcher.end(); j++) parsed[j] = true;	
            }
            
            if(firstDatas[i].endsWith(".bin") && !secondDatas[i].endsWith(".bin") ||
                    !firstDatas[i].endsWith(".bin") && secondDatas[i].endsWith(".bin")) {
                throw new IllegalStateException("Parse error: All datasets in a range must be either files or folders.");
            }




            //Find the range, if any, associated with first data.  It is of the form (3-7)
            pattern = Pattern.compile(firstDatas[i] + "\\(.*?\\)");
            matcher = pattern.matcher(rawDataSources[i]);

            if(matcher.find()) {
                firstRanges[i] = matcher.group().replaceAll(firstDatas[i], "");
                for(int j=matcher.start(); j< matcher.end(); j++) parsed[j] = true;	
            }


            //Find the range, if any, associated with second data.  It is of the form (3-7)
            if(secondDatas[i].length() > 0) {
                pattern = Pattern.compile(secondDatas[i] + "\\(.*?\\)");
                matcher = pattern.matcher(rawDataSources[i]);

                if(matcher.find()) {
                    secondRanges[i] = matcher.group().replaceAll(secondDatas[i], "");
                    for(int j=matcher.start(); j< matcher.end(); j++) parsed[j] = true;	
                }
            }



            //Check for dash, in a roundabout way.
            if(secondDatas[i].length() > 0) {
                pattern = Pattern.compile("\\Q" + firstDatas[i] + firstRanges[i] + "-" + secondDatas[i] + secondRanges[i] + "\\E");
                matcher = pattern.matcher(rawDataSources[i]);

                if(matcher.find()) {

                    matcher.group();
                    for(int j=matcher.start(); j< matcher.end(); j++) parsed[j] = true;	
                }
            }


            //Throw error for any character that was not successfully parsed.
            for(int j=0; j<parsed.length; j++) {
                if(!parsed[j]) {
                    throw new PatternSyntaxException("Unable to parse input", rawDataSources[i], j);
                }

            }



            if(firstRanges[i].length() == 0) {
                startSeconds[i] = 0;
                stopSeconds[i] = Double.POSITIVE_INFINITY;
            } else {
                double[] range = parseRange(firstRanges[i]);
                if(secondRanges[i].length() ==0 ) {
                    startSeconds[i] = range[0];
                    stopSeconds[i] = range[1];
                } else {
                    startSeconds[i] = range[0];
                    if(range[1] != Double.POSITIVE_INFINITY) {
                        throw new IllegalStateException("Parse error: You cannot have a range of the form data000(0-x) - data005(y-5).  The x and y must be removed, or you must use a comma to separate data000 and data005.");
                    }

                    range = parseRange(secondRanges[i]);
                    if(range[0] != 0) {
                        throw new IllegalStateException("You cannot have a range of the form data000(0-x) - data005(y-5).  The x and y must be removed, or you must use a comma to separate data000 and data005.");
                    }
                    stopSeconds[i] = range[1];
                }
            }


        }



        int datasetCount = 0;
        for(int i=0; i<firstDatas.length; i++) {
            if(secondDatas[i].length() == 0) {
                datasetCount++;
                datasList.add(filePath + firstDatas[i]);
                startList.add(new Double(startSeconds[i]));
                stopList.add(new Double(stopSeconds[i]));
            } else {
                int first = (new Integer(firstDatas[i].substring(4,7))).intValue();
                int last = (new Integer(secondDatas[i].substring(4,7))).intValue();
                //if a file
                if(firstDatas[i].endsWith(".bin")) {
                    
                    for(int j=first; j<=last; j++) {
                        datasetCount++;
                        if(j==first) {
                            startList.add(new Double(startSeconds[i]));
                        } else {
                            startList.add(new Double(0.0));
                        }
                        if(j==last) {
                            stopList.add(new Double(stopSeconds[i]));
                        } else {
                            stopList.add(new Double(Double.POSITIVE_INFINITY));
                        }
                        String nString = new Integer(j).toString();
                        if(nString.length() == 1) nString = "data00" + nString + ".bin";
                        else if(nString.length() == 2) nString = "data0" + nString + ".bin";
                        else if(nString.length() == 3) nString = "data" + nString + ".bin";
                        else throw new IllegalStateException(firstDatas[i] + " is not a legal dataset name.");
                        
                        datasList.add(filePath + nString);
                    }



                } else {

                    for(int j=first; j<=last; j++) {
                        datasetCount++;
                        if(j==first) {
                            startList.add(new Double(startSeconds[i]));
                        } else {
                            startList.add(new Double(0.0));
                        }
                        if(j==last) {
                            stopList.add(new Double(stopSeconds[i]));
                        } else {
                            stopList.add(new Double(Double.POSITIVE_INFINITY));
                        }
                        String nString = new Integer(j).toString();
                        if(nString.length() == 1) nString = "data00" + nString;
                        else if(nString.length() == 2) nString = "data0" + nString;
                        else if(nString.length() == 3) nString = "data" + nString;
                        else throw new IllegalStateException(firstDatas[i] + " is not a legal dataset name.");
                        datasList.add(filePath + nString);




                    }
                }
            }
        }

        for(int i=0; i<datasList.size(); i++) {
            String dataset = datasList.get(i);
            double start = startList.get(i).doubleValue();
            double stop = stopList.get(i).doubleValue();
            System.out.println(dataset + " " + start + " " + stop);
        }



    }

    public static double[] parseRange(String rangeString) {
        double[] range = new double[2];


        //parse range
        if(rangeString.length() > 0) {
            Pattern pattern = Pattern.compile("\\([\\d]*-[\\d]*\\)");
            Matcher matcher = pattern.matcher(rangeString);

            if(matcher.find()) {

                String[] temp = matcher.group().split("-");
                temp[0] = temp[0].substring(1);
                temp[1] = temp[1].substring(0, temp[1].length()-1);



                if(temp[0].length() > 0) {
                    range[0] = (new Double(temp[0])).doubleValue();
                } else {
                    range[0] = 0;
                }

                if(temp[1].length() > 0) {
                    range[1] = (new Double(temp[1])).doubleValue();
                } else {
                    range[1] = Double.POSITIVE_INFINITY;
                }
            } else {
                throw new IllegalStateException("Unable to Parse Range: " + rangeString);
            }

            //pattern = Pattern.compile("\\(-[\\d]*");
            //patcher = 
        } else {
            range[0] = 0;
            range[1] = Double.POSITIVE_INFINITY;
        }


        return range;
    }
}