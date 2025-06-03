#!/bin/bash

# USAGE: bash online_sta.sh 20221117C data015 1278606431 FastNoise

DEFAULT_PROTOCOL="FastNoise"

# litke_bin_path=$1
# shift
EXP=$1
shift
FILE_NAME=$1
shift
SEED=$1
shift
PROTOCOL_ID="${1:-$DEFAULT_PROTOCOL}"

# Get the start time.
START_TIME=$(date +%s)


# Set the paths by the computer.
# if [[ "$HOSTNAME" == *"hyak"* || "$HOSTNAME" == *"klone"* ]]; then
#     YASS_SPIKE_PATH='/gscratch/scrubbed/retina/data/sorted'; #"/home/mike/ftp/files/ksort_out/"; #'/home/mike/ftp/files/20220420C/data021/kilosort2/';
#     RAW_DATA_PATH='/gscratch/scrubbed/retina/data/raw';
#     VISIONPATH="/gscratch/retina/GitRepos/manookin-lab/MEA/src/Vision7_for_2015DAQ/Vision.jar";
#     TMP_PATH="/gscratch/scrubbed/retina/data/sorted/${EXP}/"
# else
#     YASS_SPIKE_PATH='/usr/share/pool/data/sorted'; #"/home/mike/ftp/files/ksort_out/"; #'/home/mike/ftp/files/20220420C/data021/kilosort2/';
#     KILOSORT_TTL_PATH='/usr/share/pool/data/sorted'; #'/home/mike/ftp/files/ksort_out/';
#     RAW_DATA_PATH='/usr/share/pool/data/raw';
#     VISIONPATH="/home/mike/Documents/git_repos/manookin-lab/MEA/src/Vision7_for_2015DAQ/Vision.jar";
#     TMP_PATH="/usr/share/pool/data/sorted/${EXP}/"
# fi
EXP='20230511C'
FILE_NAME='data002'
RAW_DATA_PATH='/data/data/raw';
TMP_PATH="/data/data/sorted/${EXP}/"
source ../utilities/data_paths.sh

SORT_PATH="${YASS_SPIKE_PATH}/${EXP}/"
DATA_PATH="${RAW_DATA_PATH}/${EXP}/"


data_string="${DATA_PATH}$FILE_NAME"

if [[ ! -e $TMP_PATH ]]; then
   mkdir ${TMP_PATH}
fi

python ../utilities/convert_litke_to_yass.py ${data_string} ${TMP_PATH} ${FILE_NAME} -w

# Step 2: Pre-process the data (subtract the average across all channels)
# now run MATLAB script to actually run Kilosort
/usr/local/bin/matlab -nodesktop -nosplash -r "run_kilosort_single('${EXP}','${FILE_NAME}'); quit"
rc=$?; if [[ $rc != 0 ]]; then echo "Pre-processing failed, exiting..."; exit $rc; fi

# Convert to Vision format.
python ../utilities/kilosort_to_vision.py ${TMP_PATH}${FILE_NAME}/kilosort2 ${TMP_PATH}${FILE_NAME} ${RAW_DATA_PATH}/${EXP}/${FILE_NAME}/ ${TMP_PATH}${FILE_NAME} ${FILE_NAME} -l

# Get the EI file.
java -Xmx8G -cp $VISIONPATH edu.ucsc.neurobiology.vision.calculations.CalculationManager "Electrophysiological Imaging Fast" ${TMP_PATH}${FILE_NAME}/ ${RAW_DATA_PATH}/${EXP}/${FILE_NAME}/ 0.01 20 40 1000000 8

# Run the simple STA analysis.
python sta_single_file.py ${TMP_PATH} ${FILE_NAME} ${SEED} -p ${PROTOCOL_ID}

END_TIME=$(date +%s)
RUN_TIME=$((END_TIME-START_TIME))
MINUTES_TIME=$((RUN_TIME/60))
SECONDS_TIME=$((RUN_TIME-MINUTES_TIME*60))
echo "Completed running after ${MINUTES_TIME} minutes and ${SECONDS_TIME} seconds."