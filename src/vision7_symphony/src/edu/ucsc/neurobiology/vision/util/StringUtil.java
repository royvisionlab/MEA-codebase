package edu.ucsc.neurobiology.vision.util;

import java.text.*;
import java.util.*;


/**
 * This class is intended to contain miscelaneous string utility
 * methods used allover the program.
 *
 * @author Dumitru Petrusca, University of California, Santa Cruz
 */
public final class StringUtil {
    private static DecimalFormat format = new DecimalFormat();

    private StringUtil() {}


    /**
     * Extracts the extension (with the dot) from a given filename.
     * @return the extension (with the dot)
     */
    public static String getExtension(String fileName) {
        int lastIndex = fileName.lastIndexOf('.');
        if (lastIndex == -1) {
            return ".";
        } else {
            return fileName.substring(fileName.lastIndexOf('.'));
        }
    }


    /**
     * Returns a string representation of the array. The values are separated by commas.
     */
    public static String toString(int[] x) {
        String s = "";

        for (int i = 0; i < x.length; i++) {
            s += x[i] + ", ";
        }

        return s;
    }


    /**
     * Takes a filename and returns the name without the extension
     */
    public static String removeExtension(String fileName) {
        StringBuffer s = new StringBuffer(fileName);
        for (int i = s.length()-1; i > -1; i--) {
            char c = s.charAt(i);
            if (c == '.') {
                s.setLength(i);
                break;
            }
            else if (c == '\\') break;
        }
        return s.toString();
    }



    /**
     * Replaces the extension of the given filename with the given string.
     */
    public static String setExtension(String fileName, String newExtension) {
        return removeExtension(fileName) + newExtension;
    }


    // Should probably be deprecated in favor of Java 1.4 String#split
    public static String[] decomposeString(String string, String delimiter) {
        StringTokenizer st = new StringTokenizer(string, delimiter, false);
        ArrayList<String> names = new ArrayList<String>();

        for (int i = 0; st.hasMoreTokens(); i++) {
            String s = st.nextToken().trim();
            if (s.length() > 0) {
                names.add(s);
            }
        }

        return (String[]) names.toArray(new String[names.size()]);
    }
    
    // Could be replaced by Apache StringUtils.join
    public static String joinStringCollection(Collection<String> strings, String delim) {
        if (strings.isEmpty()) return "";
        Iterator<String> it = strings.iterator();
        String retstr = it.next();
        while (it.hasNext()) retstr += delim + it.next();
        return retstr;
    }


    /**
     * Converts a double value to a string keeping a minimum of 0 and a maximum of 3
     * fractionar digits.
     */
    public static final String format(double v) {
        format.setMaximumFractionDigits(3);
        format.setMinimumFractionDigits(0);
        format.setMinimumIntegerDigits(0);
        format.setMaximumIntegerDigits(Integer.MAX_VALUE);

        return format.format(v);
    }


    /**
     * Converts a double value to a string keeping the required number of
     * fractionar digits.
     */
    public static final String format(double v, int fractionarDigits) {
        format.setMaximumFractionDigits(fractionarDigits);
        format.setMinimumFractionDigits(fractionarDigits);
        format.setMinimumIntegerDigits(0);
        format.setMaximumIntegerDigits(Integer.MAX_VALUE);
        return format.format(v);
    }


    /**
     * Converts a double value to a string keeping the required number of
     * fractionar digits.
     */
    public static final String format(double v, int fractionarDigits, int intgerDigits) {
        format.setMaximumFractionDigits(fractionarDigits);
        format.setMinimumFractionDigits(fractionarDigits);
        format.setMinimumIntegerDigits(intgerDigits);
//        format.setMaximumIntegerDigits(intgerDigits);
        format.setDecimalSeparatorAlwaysShown(false);
        format.setGroupingSize(100);
        return format.format(v);
    }

    
    public static final String DEFAULT_SEP = ", ";

    public static final String listString(List<?> lst) {
        return listString(lst, DEFAULT_SEP);
    }
    
    public static final String listString(List<?> lst, String sep) {
        String retstr = "";
        for (int i = 0; i < lst.size()-1; i++) {
            retstr += lst.get(i).toString() + sep;
        }
        retstr += lst.get(lst.size()-1).toString();
        return retstr;
    }

    public static final String listString(Object[] arr) {
        return listString(arr, DEFAULT_SEP);
    }
    
    public static final String listString(Object[] arr, String sep) {
        return listString(Arrays.asList(arr), sep);
    }
}
