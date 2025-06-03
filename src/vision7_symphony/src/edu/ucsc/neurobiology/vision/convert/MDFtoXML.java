package edu.ucsc.neurobiology.vision.convert;

import java.io.*;
import java.util.*;
import edu.ucsc.neurobiology.vision.util.*;
import edu.ucsc.neurobiology.vision.*;

/**
 *  Tool to convert legacy movie definition files (mdf) to
 *  xml files.  The new xml files serve the same purpose as the old
 *  mdf files, but are easier to work with and extend. 
 * 
 * 
 * 
 * 
 */

public class MDFtoXML {


    public MDFtoXML(String mdfPath, String xmlPath) throws Exception {
        String xmlTemplate = xmlPath + "template.xml";
        File mdfDir = new File(mdfPath);

        File[] mdfFiles = mdfDir.listFiles();

        for(int i=0; i<mdfFiles.length; i++) {

            String mdfName = mdfFiles[i].getName();

            if(mdfName.startsWith(".")) continue;
            if(mdfName.startsWith("2007")) continue;
            if(mdfName.startsWith("template")) continue;

            LinkedHashMap<String, String> mdfHash = mdfToHashMap(new File(mdfPath+mdfName));

            String xmlName = StringUtil.removeExtension(mdfName) + ".xml";

            BZip2Compress.copyFile(new File(xmlTemplate), new File(xmlPath + xmlName));

            Config movieConfig = new Config(xmlPath + xmlName);
            LinkedHashMap<String, String> makeHash = movieConfig.getParameterList("Make White Noise Movie");
            LinkedHashMap<String, String> auxHash = movieConfig.getParameterList("Calculate Auxiliary Parameters.Set Movie");

            Set<String> mdfKeys = mdfHash.keySet();
            for(String key: mdfKeys) {
                String value = mdfHash.get(key);



                if (key.matches("MovieType")) {
                    continue;
                }

                if (key.matches("NoiseType")) {
                    if(value.matches("bin")) {
                        makeHash.put("NoiseType", "0.0");
                    } else if(value.matches("gauss")) {
                        makeHash.put("NoiseType", "1.0");
                    } else {
                        throw new IllegalStateException("Key:Value " + key + ":" + value + " not recognized.");
                    }
                    continue;
                }

                if (key.matches("ColorType")) {
                    if(value.matches("rgb")) {
                        makeHash.put("ColorType", "0.0");
                    } else if(value.matches("gray")) {
                        makeHash.put("ColorType", "1.0");
                    } else if(value.matches("separated")) {
                        makeHash.put("ColorType", "2.0");
                    }
                    else {
                        throw new IllegalStateException("Key:Value " + key + ":" + value + " not recognized.");
                    }
                    continue;
                }

                if (key.matches("RandomNumberGenerator")) {
                    if (value.matches("mac")) {
                        makeHash.put("RandomNumberGenerator", "0.0");
                    } else if(value.matches("java")) {
                        makeHash.put("RandomNumberGenerator", "1.0");
                    } else if(value.matches("javaV2")) {
                        makeHash.put("RandomNumberGenerator", "2.0");
                    }
                    else {
                        throw new IllegalStateException("Key:Value " + key + ":" + value + " not recognized.");
                    }
                    continue;
                }

                if (key.matches("Seed")) {
                    makeHash.put(key, value);
                    continue;
                }

                if (key.matches("Width")) {
                    makeHash.put(key, value);
                    continue;
                }

                if (key.matches("Height")) {
                    makeHash.put(key, value);
                    continue;
                }

                if (key.matches("RefreshInterval")) {
                    auxHash.put("refreshInterval", value);
                    continue;
                }

                if(key.matches("ContrastSigma")) {
                    makeHash.put(key, value);
                    continue;
                }

                if(key.matches("PixelSize")) {
                    auxHash.put("pixelsPerStixelX", new Double(value).intValue() + "");
                    auxHash.put("pixelsPerStixelY", new Double(value).intValue() + "");
                    continue;
                }



                throw new IllegalStateException("Key " + key + " not found.");

            }


            //if pixel size is not in file, use the line from the file name.  It is the first int.
            if(!mdfHash.containsKey("PixelSize")) {
                String[] split = mdfName.split("-");
                for(int j=0; j<split.length; j++) {
                    boolean failed = false;
                    try {
                        new Integer(split[j]);
                    } catch (NumberFormatException ex) {
                        failed = true;
                    }
                    if(!failed) {
                        //zero means full screen
                        if(!split[j].matches("0")) {
                            auxHash.put("pixelsPerStixelX", split[j]);
                            auxHash.put("pixelsPerStixelY", split[j]);
                        } else {
                            auxHash.put("pixelsPerStixelX", "640");
                            auxHash.put("pixelsPerStixelY", "320");
                        }
                        break;
                    }
                }

            }

            //determine X and Y offset, assuming that the stimulus is centered on a 640x480 display

            int pixelsPerStixelX = new Double(auxHash.get("pixelsPerStixelX")).intValue();
            int pixelsPerStixelY = new Double(auxHash.get("pixelsPerStixelY")).intValue();
            int width = new Double(makeHash.get("Width")).intValue();
            int height = new Double(makeHash.get("Height")).intValue();
            auxHash.put("xOffset", (640 - pixelsPerStixelX*width)/2 + "");
            auxHash.put("yOffset", (480 - pixelsPerStixelY*height)/2 + "");


            movieConfig.setParameterGroup("Make White Noise Movie", makeHash);
            movieConfig.setParameterGroup("Calculate Auxiliary Parameters.Set Movie", auxHash);



        }





    }

    private LinkedHashMap<String,String> mdfToHashMap(File mdfFile) throws Exception {
        BufferedReader input = new BufferedReader(new FileReader(mdfFile));
        String line = null;
        LinkedHashMap<String, String> mdfHash = new LinkedHashMap<String, String>();
        while((line = input.readLine()) != null) {
            if(line.contains("=")) {
                String[] parsed = line.split("=");
                mdfHash.put(parsed[0].trim(), parsed[1].trim());

            }
        }

        return mdfHash;
    }


    public static void main(String[] args) {
        try {
            MDFtoXML mdfToXml = new MDFtoXML("//Volumes//snleacquisition//mdf//", "//Volumes//snleacquisition//movie-xml//");
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

}
