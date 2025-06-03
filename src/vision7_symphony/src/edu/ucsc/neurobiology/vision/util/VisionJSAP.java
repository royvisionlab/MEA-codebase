/*
 * Copyright (c) 2002-2004, Martian Software, Inc.
 * This file is made available under the LGPL as described in the accompanying
 * LICENSE.TXT file.
 */
package edu.ucsc.neurobiology.vision.util;

import com.martiansoftware.jsap.*;

import java.util.Iterator;
import java.util.List;

import com.martiansoftware.util.StringUtils;

/** A simple interface to {@link com.martiansoftware.jsap.JSAP} that handles directly help,
 * explanation and an array of parameters.
 * 
 * <P>More precisely, instances of this class behave exactly like those of
 * {@link com.martiansoftware.jsap.JSAP}, but additionally require a command name and an array 
 * of parameters (which will be registered automatically). 
 * A switch activated by <samp>--help</samp> is always
 * registered under the ID <samp>help</samp>. 
 * 
 * <p>A message will be automatically printed upon invocation
 * of the <code>parse()</code> methods if an error occurs, or if the help switch is detected. In this
 * case, the code will exit.
 * 
 * See http://martiansoftware.com/jsap/ for full documentation.  Mutated from SimpleJSAP by Matthew Grivich.
 * 
 * 
 * @author Sebastiano Vigna
 * @author Matthew Grivich, The Salk Institute
 */
public class VisionJSAP extends JSAP {

    /** The screen width used for formatting. */
    private int screenWidth = 100;


    /** The name of the command that will appear in the help message. */
    final private String name;

    /** Creates a new simple JSAP with default screen width. 
     * 
     * @param name the name of the command for which help will be printed.
     * @param parameter an array of parameters, which will be registered for you, or <code>null</code>.
     */
    public VisionJSAP( final String name, final Parameter[] parameter ) throws JSAPException {
        super();

        this.name = "java -Xmx1564m -Xss1m -classpath Vision.jar " + name;


        final Switch help = new Switch( "help", 'h', "help" );
        help.setHelp( "Prints this help message." );
        this.registerParameter( help );

        if ( parameter != null ) 
            for( int i = 0; i < parameter.length; i++ ) this.registerParameter( parameter[ i ] );
    }

    public JSAPResult parse( String arg ) {
        JSAPResult jsapResult = super.parse( arg );
        printMessageIfUnsuccessfulOrHelpRequired( jsapResult );
        return jsapResult;
    }

    public JSAPResult parse( String[] arg ) {
        JSAPResult jsapResult = super.parse( arg );
        printMessageIfUnsuccessfulOrHelpRequired( jsapResult );
        return jsapResult;
    }

    // PHLI: Abstracted to here from original inline in CalculationManager.main()
    //Hack around because JSAP cannot handle unflagged negative numbers.
    //This has been on the bug list for JSAP since 9/06.  I'm (mgrivich) writing this on 5/07.	
    public JSAPResult parseWithMangledNegs(String[] args, String mangle) {
        // Hack around part I
        for (int i = 0; i < args.length; i++) {
            boolean isNeg = false;
            if (args[i].startsWith("-")) {
                isNeg = true;
                try {
                    new Double(args[i]);
                } catch (Exception ex) {
                    isNeg = false;
                }
                if (isNeg) {
                    args[i] = args[i].replaceFirst("-", mangle);
                }
            }
        }
        return parse(args);
    }
    
    public JSAPResult parseWithMangledNegs(String[] args) {
        return parseWithMangledNegs(args, "negative");
    }

    // PHLI: Abstracted to here from original inline in CalculationManager.main()
    //Hack around part II.
    public static String[] demangle(String[] mangled, String mangle) {
        for (int i = 0; i < mangled.length; i++) {
            if (mangled[i].startsWith("negative")) {
                mangled[i] = mangled[i].replaceFirst(mangle, "-");
            }
        }		
        return mangled;
    }

    public static String[] demangleNegs(String[] mangled) {
        return demangle(mangled, "negative");
    }
    
    /** Checks the given JSAP result for errors or help requests and acts accordingly.
     * 
     * @param jsapResult the result of a JSAP parsing.
     */

    private void printMessageIfUnsuccessfulOrHelpRequired( final JSAPResult jsapResult ) {
        if (  !jsapResult.success() || jsapResult.getBoolean( "help" ) ) {

            System.err.println();
            for ( Iterator err = jsapResult.getErrorMessageIterator(); err.hasNext(); ) 
                System.err.println( "Error: " + err.next() );

            System.err.println();


            System.err.println( "Usage:" );
            List l = StringUtils.wrapToList( name + " " + getUsage(), screenWidth );


            for( Iterator i = l.iterator(); i.hasNext(); ) 
                System.err.println( "  " + i.next().toString() );
            
            System.err.println();
            System.err.println();

            System.err.println("Details:");
            System.err.println( getHelp( screenWidth ) );

            System.exit(1);

        }

    }

}