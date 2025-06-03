package edu.ucsc.neurobiology.vision.io;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.InvocationTargetException;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import edu.ucsc.neurobiology.vision.util.StringUtil;
import edu.ucsc.neurobiology.vision.util.VisionParams;

public class VisionFileInfo {
    
    public static void main(String[] args) throws IOException, NoSuchMethodException, InvocationTargetException, IllegalArgumentException, InstantiationException, IllegalAccessException {
        String filename = args[0];
        File f = new File(filename);
        if (!f.exists()) throw new Error("File " + filename + " not found!");
            
        if (f.isDirectory()) {
            String[] fileMatches = f.list(new visInfoFilenameFilter());
            for (String match : fileMatches) {
                System.out.println(match + ":");
                printVisInfo(filename + "/" + match);
                System.out.println("");
            }
            return;
        }

        printVisInfo(filename);
    }
    
    
    static private class visInfoFilenameFilter implements FilenameFilter {
        final private static Set<String> EXTENSIONS = VisionParams.EXTENSION_MAP.keySet();
        public boolean accept(File dir, String name) {
            String ext = StringUtil.getExtension(name);
            if (!EXTENSIONS.contains(ext)) return false;
            
            // Special handling of .bin; only include the first bin file
            if (ext.equals(VisionParams.BIN_FILE_EXTENSION_512)) {
                try { return isFirstBinFile(name); }
                catch (IllegalArgumentException iae) { 
                    // If something unexpected happened, let it through
                    System.err.println(iae.getMessage()); 
                    return true;
                }
            }
            
            return true;
        };
        
        static private final Pattern BIN_FILE_PAT = Pattern.compile("data\\d{3}(\\d{3})?.bin");
        private static boolean isFirstBinFile(String filename) {
            Matcher m = BIN_FILE_PAT.matcher(filename);
            if (!m.find()) throw new IllegalArgumentException("Unrecognized .bin file naming convention");
            String g = m.group(1);
            if (g == null || g.equals("000")) return true;
            return false;
        }
    }
    
    
    private static void printVisInfo(String filename) throws IllegalArgumentException, SecurityException, InstantiationException, IllegalAccessException, InvocationTargetException, NoSuchMethodException {
        Object visFile = buildVisFile(filename);
        if (visFile != null) System.out.println(visFile);
        else                 System.err.println("Unrecognized file type: " + filename);
    }
    
    private static Object buildVisFile(String filename) throws IllegalArgumentException, InstantiationException, IllegalAccessException, InvocationTargetException, SecurityException, NoSuchMethodException {
        Object visFile = null;
        String ext = StringUtil.getExtension(filename);
        for (String k : VisionParams.EXTENSION_MAP.keySet()) {
            if (ext.equals(k)) {
                Constructor<?> ctor = VisionParams.EXTENSION_MAP.get(k).getDeclaredConstructor(String.class);
                visFile = ctor.newInstance(filename);
                break;
            }
        }
        
        return visFile;
    }
}
