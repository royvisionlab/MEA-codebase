@echo off

set CP=.;..\classes;
set CP=%CP%;..\..\java-library\colt.jar;
set CP=%CP%;..\..\java-library\gif.jar;
set CP=%CP%;..\..\java-library\java_cup.jar;
set CP=%CP%;..\..\java-library\jh.jar;
set CP=%CP%;..\..\java-library\jmf.jar;
set CP=%CP%;..\..\java-library\jx.jar;
set CP=%CP%;..\..\java-library\junit-4.0.jar;

javadoc -source 1.5 -classpath %CP% -sourcepath ..\src -d doc -author -link http://xml.apache.org/xerces-j/apiDocs/ -link http://java.sun.com/j2se/1.3/docs/api/ -link http://java.freehep.org/lib/freehep/api/ -windowtitle "Vision API Documentation" @packages.txt
