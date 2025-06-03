#!/bin/bash

# USAGE: bash online_vision.sh 20221117C data015 1278606431 FastNoise

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: $0 <EXPERIMENT_DATE> <FILE_NAME> <SEED> <PROTCOL_ID>"
    echo "Example w/ default protocol: $0 20240926C data015 1278606431"
    echo "Example: $0 20240926C data015 1278606431 504"
    exit 0
fi


DEFAULT_PROTOCOL="FastNoise"

EXPERIMENT_NAME=$1
shift
FILE_NAME=$1
shift
START_SEED=$1
shift
PROTOCOL_ID="${1:-$DEFAULT_PROTOCOL}"

# Get the start time.
START_TIME=$(date +%s)


# Set the paths by the computer.
if [[ "$HOSTNAME" == *"hyak"* || "$HOSTNAME" == *"klone"* ]]; then
    YASS_SPIKE_PATH='/gscratch/retina/data/sorted'; #"/home/mike/ftp/files/ksort_out/"; #'/home/mike/ftp/files/20220420C/data021/kilosort2/';
    RAW_DATA_PATH='/gscratch/scrubbed/retina/data/raw';
    VISIONPATH="/gscratch/retina/GitRepos/manookin-lab/MEA/src/Vision7_for_2015DAQ/";
    TMP_PATH="/gscratch/scrubbed/retina/data/sorted/${EXPERIMENT_NAME}/"
else
    YASS_SPIKE_PATH='/data/data/sorted'; #"/home/mike/ftp/files/ksort_out/"; #'/home/mike/ftp/files/20220420C/data021/kilosort2/';
    KILOSORT_TTL_PATH='/data/data/sorted'; #'/home/mike/ftp/files/ksort_out/';
    RAW_DATA_PATH='/data/data/raw';
    VISIONPATH="/home/mike/Documents/git_repos/manookin-lab/MEA/src/Vision7_for_2015DAQ/";
    TMP_PATH="/data/data/sorted/${EXPERIMENT_NAME}/"
fi

# YASS_SPIKE_PATH='/usr/share/pool/data/sorted'; #"/home/mike/ftp/files/ksort_out/"; #'/home/mike/ftp/files/20220420C/data021/kilosort2/';
    # KILOSORT_TTL_PATH='/usr/share/pool/data/sorted'; #'/home/mike/ftp/files/ksort_out/';
    # RAW_DATA_PATH='/usr/share/pool/data/raw';
    # VISIONPATH="/home/mike/Documents/git_repos/manookin-lab/MEA/src/Vision7_for_2015DAQ/";
    # TMP_PATH="/usr/share/pool/data/sorted/${EXPERIMENT_NAME}/"

SORT_PATH="${YASS_SPIKE_PATH}/${EXPERIMENT_NAME}/${FILE_NAME}/vision/"
DATA_PATH="${RAW_DATA_PATH}/${EXPERIMENT_NAME}/"


data_string="${DATA_PATH}$FILE_NAME"

if [[ ! -e $TMP_PATH ]]; then
   mkdir -p ${TMP_PATH}
fi

if [[ ! -e $SORT_PATH ]]; then
   mkdir -p ${SORT_PATH}
fi

# Run spike sorting in Vision.
java -d64 -Xmx8000m -Xss2m -cp ${VISIONPATH}Vision.jar edu.ucsc.neurobiology.vision.tasks.NeuronIdentification ${DATA_PATH}${FILE_NAME} ${SORT_PATH} -c ${VISIONPATH}config.xml

# Convert to Vision format.
# python ../utilities/kilosort_to_vision.py ${SORT_PATH}${FILE_NAME}/kilosort2 ${SORT_PATH}${FILE_NAME} ${SORT_PATH}${FILE_NAME} ${FILE_NAME} -l -d ${RAW_DATA_PATH}/${EXPERIMENT_NAME}/${FILE_NAME}/

# Get the EI file.
java -Xmx8G -cp ${VISIONPATH}Vision.jar edu.ucsc.neurobiology.vision.calculations.CalculationManager "Electrophysiological Imaging Fast" ${SORT_PATH}/ ${RAW_DATA_PATH}/${EXPERIMENT_NAME}/${FILE_NAME}/ 0.01 20 40 1000000 8

# Run the simple STA analysis.
# python sta_single_file.py ${SORT_PATH} ${FILE_NAME} ${START_SEED} -p ${PROTOCOL_ID}
python sta_single_file.py ${SORT_PATH} vision ${START_SEED} -p ${PROTOCOL_ID}

END_TIME=$(date +%s)
RUN_TIME=$((END_TIME-START_TIME))
MINUTES_TIME=$((RUN_TIME/60))
SECONDS_TIME=$((RUN_TIME-MINUTES_TIME*60))
echo "Completed running after ${MINUTES_TIME} minutes and ${SECONDS_TIME} seconds."