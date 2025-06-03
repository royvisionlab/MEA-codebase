#!/usr/bin/env bash
# Creates the Vision.app file.
# Will only work when ran from OSX.

# Check for JAVA_HOME variable - used to bundle jdk
javahomest=$JAVA_HOME
if [[ "$javahomest" == "" ]]; then
    echo "[-] JAVA_HOME is not properly set on this machine."
    echo "    JAVA_HOME should be set to something like '/usr/libexec/java_home' on OSX."
    echo "    Please update your .bash_profile (export JAVA_HOME = [...])."
    echo "    Bundling will fail, aborting now."
    exit
else
    echo "[x] Found JAVA_HOME:" $javahomest
fi

unamestr=`uname`
if [[ "$unamestr" == 'Darwin' ]]; then
    echo "    Bundling app... "
    ant -f bundle-app.xml
    echo "[x] Bundling done."
else
    echo "[-] App bundling will only work on OSX."
    echo "    Aborting bundling now."
    exit
fi