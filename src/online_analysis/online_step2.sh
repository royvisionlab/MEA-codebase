#!/bin/bash

# USAGE: bash online_step1.sh 20221117C data015 1278606431 FastNoise

DEFAULT_PROTOCOL="FastNoise"

# litke_bin_path=$1
# shift
experimentName=$1
shift
out_name=$1
shift
start_seed=$1
shift
PROTOCOL_ID="${1:-$DEFAULT_PROTOCOL}"

# Get the start time.
START_TIME=$(date +%s)


# Set the paths by the computer.
if [[ "$HOSTNAME" == *"hyak"* || "$HOSTNAME" == *"klone"* ]]; then
    YASS_SPIKE_PATH='/gscratch/retina/data/sorted'; #"/home/mike/ftp/files/ksort_out/"; #'/home/mike/ftp/files/20220420C/data021/kilosort2/';
    RAW_DATA_PATH='/gscratch/scrubbed/retina/data/raw';
    VISIONPATH="/gscratch/retina/GitRepos/manookin-lab/MEA/src/Vision7_for_2015DAQ/Vision.jar";
    TMP_PATH="/gscratch/scrubbed/retina/data/sorted/${experimentName}/"
else
    YASS_SPIKE_PATH='/usr/share/pool/data/sorted'; #"/home/mike/ftp/files/ksort_out/"; #'/home/mike/ftp/files/20220420C/data021/kilosort2/';
    KILOSORT_TTL_PATH='/usr/share/pool/data/sorted'; #'/home/mike/ftp/files/ksort_out/';
    RAW_DATA_PATH='/usr/share/pool/data/raw';
    VISIONPATH="/home/mike/Documents/git_repos/manookin-lab/MEA/src/Vision7_for_2015DAQ/Vision.jar";
    TMP_PATH="/usr/share/pool/data/sorted/${experimentName}/"
fi

SORT_PATH="${YASS_SPIKE_PATH}/${experimentName}/"
DATA_PATH="${RAW_DATA_PATH}/${experimentName}/"


data_string="${DATA_PATH}$out_name"

if [[ ! -e $TMP_PATH ]]; then
   mkdir ${TMP_PATH}
fi

# Convert to Vision format.
python ../utilities/kilosort_to_vision.py ${SORT_PATH}${out_name}/kilosort2 ${SORT_PATH}${out_name} ${SORT_PATH}${out_name} ${out_name} -l -d ${RAW_DATA_PATH}/${experimentName}/${out_name}/

# Get the EI file.
java -Xmx8G -cp $VISIONPATH edu.ucsc.neurobiology.vision.calculations.CalculationManager "Electrophysiological Imaging Fast" ${SORT_PATH}${out_name}/ ${RAW_DATA_PATH}/${experimentName}/${out_name}/ 0.01 20 40 1000000 8

# Run the simple STA analysis.
python sta_single_file.py ${SORT_PATH}${out_name} ${out_name} ${start_seed} -p ${PROTOCOL_ID}

END_TIME=$(date +%s)
RUN_TIME=$((END_TIME-START_TIME))
MINUTES_TIME=$((RUN_TIME/60))
SECONDS_TIME=$((RUN_TIME-MINUTES_TIME*60))
echo "Completed running after ${MINUTES_TIME} minutes and ${SECONDS_TIME} seconds."