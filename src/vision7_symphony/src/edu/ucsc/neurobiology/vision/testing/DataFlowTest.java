package edu.ucsc.neurobiology.vision.testing;

import java.io.File;
import java.util.*;
import edu.ucsc.neurobiology.vision.tasks.*;
import edu.ucsc.neurobiology.vision.*;
import edu.ucsc.neurobiology.vision.anf.*;

/*
 * Data flow test for standard 8-way streaming 512 electrode experiment. 
 * 
 * 
 * @author Matthew Grivich, The Salk Institute
 * 
 */



public class DataFlowTest {

    private static void split(String rawData, String splitPath) throws Exception {


        HashMap<String, String> p = new HashMap<String, String>();
        p.put("Raw_Data_Source", rawData);
        p.put("Buffer Size (Kb)", "111");
        p.put("Buffers Count", "100");
        p.put("Spike Threshold", "3.0");
        p.put("TTL Threshold", "1000.0");
        p.put("Sigma", "30");
        p.put("Mean Time Constant", ".01");

        p.put("Diagnostic Plots", "false");
        p.put("saveRawData", "true");
        p.put("saveRawData.Common Path", splitPath);
        p.put("saveRawData.outputFolders", "0;1;2;3;4;5;6;7");
        p.put("saveRawData.Stop Streaming After (sec)", "-1");

        Vision.getInstance().getCalculationManager().runCalculation("Spike Finding", p);

    }

    private static void standard512() {
        String rawData = "/data/raw/2007-03-02-1/data001/";
        String splitPath = "/data/rawsplit/";
        String movieXML = "/data/raw/2007-03-02-1/data001.xml";
        String out = "/data/2007-03-02-1/data001";

        try {

            split(rawData, splitPath);

            for(int i=0; i<8; i++) {
                WhiteNoiseDatasetAnalysis.main(new String[]{splitPath + i + "/2007-03-02-1/data001.bin",
                        "/data/" + i + "/", movieXML, "-e", "-g", "-c", "defaults.xml"});
            }

            String in = "";
            for(int i=0; i<8; i++) {
                in = in + "/data/" + i + "/2007-03-02-1/data001" + ";";
            }

            Join.main(new String[]{in,out});
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private static void standard512b() {
        
        try {
    /*		File outFile = new File("/Volumes/Viola/Data/Grivich/2000-12-14-1/data051.bin");
            
              File[] files =  outFile.listFiles();
               for(int i=0; i<files.length; i++) {
                   files[i].delete();
               }
*/
            WhiteNoiseDatasetAnalysis.main(new String[]{
                    "/Volumes/Leopard Intel A/temp/Data/Grivich/raw/2009-02-28-2/data007",
                    "/Volumes/Viola/Data/Grivich/",
                    "/Volumes/Viola/Data/Grivich/xml/RGB-10-32-0.48-11111.xml",
                    "-g", 
                    "-c", "defaults.xml"} );
        /*	WhiteNoiseDatasetAnalysis.main(new String[]{
                "/Data/Machado/2009-04-13-0/data000/data000.bin",
                "/Data/Machado/2009-04-13-0/data000/test2/",
                "/snle/lab/Development/RRS/xml-library/versions/2008-11-17/primate.xml",
                "-eg","-u", //unwhitened
                "-c", "defaults.xml"} );
                */

        } catch(Exception e) {
            e.printStackTrace();
        }
    }

    private static void standard61() {
        String root = "/data/2000-12-14-1/";
        try {
    /*		File outFile = new File("/Volumes/Viola/Data/Grivich/2000-12-14-1/data051.bin");
            
              File[] files =  outFile.listFiles();
               for(int i=0; i<files.length; i++) {
                   files[i].delete();
               }
*/
            WhiteNoiseDatasetAnalysis.main(new String[]{
                    "/Volumes/Viola/Data/Grivich/raw/2000-12-14-1/data051.bin",
                    "/Volumes/Viola/Data/Grivich/",
                    "/Volumes/Viola/Data/Grivich/xml/data051.xml",
                    "-egu", 
                    "-c", "defaults.xml"} );
        /*	WhiteNoiseDatasetAnalysis.main(new String[]{
                "/Data/Machado/2009-04-13-0/data000/data000.bin",
                "/Data/Machado/2009-04-13-0/data000/test2/",
                "/snle/lab/Development/RRS/xml-library/versions/2008-11-17/primate.xml",
                "-eg","-u", //unwhitened
                "-c", "defaults.xml"} );
                */

        } catch(Exception e) {
            e.printStackTrace();
        }
    }
    
    private static void simple512() {

        try {
            File outFile = new File("/Volumes/Wharf/Data/Grivich/2000-12-14-1/data051.bin");
            
              File[] files =  outFile.listFiles();
               for(int i=0; i<files.length; i++) {
                   files[i].delete();
               }

            WhiteNoiseDatasetAnalysis.main(new String[]{
                    "/Volumes/Windows A/data/raw/2007-03-02-1/data001/",
                    "/Volumes/Wharf/Data/Grivich/",
                    "/Volumes/Windows A/data/raw/2007-03-02-1/data001.xml",
                    "-eg", 
                    "-c", "defaults.xml"} );
        /*	WhiteNoiseDatasetAnalysis.main(new String[]{
                "/Data/Machado/2009-04-13-0/data000/data000.bin",
                "/Data/Machado/2009-04-13-0/data000/test2/",
                "/snle/lab/Development/RRS/xml-library/versions/2008-11-17/primate.xml",
                "-eg","-u", //unwhitened
                "-c", "defaults.xml"} );
                */

        } catch(Exception e) {
            e.printStackTrace();
        }
    }
    
    private static void simple512b() {

        try {
    /*		File outFile = new File("/Data/Grivich/2010-03-05-2/data000");
            
              File[] files =  outFile.listFiles();
               for(int i=0; i<files.length; i++) {
                   files[i].delete();
               }
*/
            WhiteNoiseDatasetAnalysis.main(new String[]{
                    "/Data/Grivich/raw/2010-03-05-2/data000/",
                    "/Data/Grivich/",
                    "/Data/Grivich/xml/RGB-8-16-0.48-11111-519elect.xml",
                    "-eg", 
                    "-c", "defaults.xml"} );
        /*	WhiteNoiseDatasetAnalysis.main(new String[]{
                "/Data/Machado/2009-04-13-0/data000/data000.bin",
                "/Data/Machado/2009-04-13-0/data000/test2/",
                "/snle/lab/Development/RRS/xml-library/versions/2008-11-17/primate.xml",
                "-eg","-u", //unwhitened
                "-c", "defaults.xml"} );
                */

        } catch(Exception e) {
            e.printStackTrace();
        }
    }
    
    
    //assumes that standard 61 has been run first
    private static void splitAndMap61() {
        try {
        MappingAnalysis.main(new String[]{
                "/Volumes/Wharf/Data/Grivich/2000-12-14-1/data051.bin",
                "/Volumes/Wharf/Data/Grivich/",
                "/Volumes/Rat/Data/Grivich/raw/2000-12-14-1/data051.bin(400-800)",	
                "-c", "defaults.xml",
                "-g"
                } );
        
        
        } catch(Exception e) {
            e.printStackTrace();
        }
    
    }
    
    private static void concatenatedMappingAnalysis() {
        try {
        ConcatenatedMappingAnalysis.main(new String[] {
                "/Volumes/Twist/Data/Greschner/2009-12-03-2/data000,data001",
                "/Volumes/Twist/Data/Greschner/9999-99-99-9/data000_001-refit/data000_001",
                "100",
                ".1",
                "-c", "/Volumes/Twist/Vision.app/Contents/Resources/Java/config.xml"
                
        });
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private static void concatenatedMappingAnalysis2() {
        try {
        ConcatenatedMappingAnalysis.main(new String[] {
                "/Volumes/Scratch/Data/Max/2010-06-30-0/data006,data007,data010,data011,data012,data013,data014",
                "/Volumes/Scratch/Data/Max/2010-06-30-0/white+allLED/",
                "100",
                ".1",
                "-n",
                "-c",
                "/snle/lab/Applications/java-volatile/vision/Vision.app/Contents/Resources/Java/config.xml"
                
        });
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private static void concatenatedMappingAnalysis3() {
        try {
        ConcatenatedMappingAnalysis.main(new String[] {
                "/Volumes/Scratch/Data/Max/2010-08-03-1/data002,data003,data004,data005,data006,data007,data008,data009,data010,data011,data012,data013,data014,data015",
                "/Volumes/Scratch/Data/Max/2010-08-03-1/white_allLED/",
                "100",
                ".1",
                "-n",
                "-c",
                "/snle/lab/Applications/java-volatile/vision/Vision.app/Contents/Resources/Java/config.xml"
                
        });
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        concatenatedMappingAnalysis3();

    }
}
