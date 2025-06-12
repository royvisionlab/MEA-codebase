#!/bin/bash

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: $0 <EXPERIMENT_DATE> <CHUNK_NAME> <SORT_ALGORITHM> <NUM_CPU>"
    echo "Example: bash $0 20240926C chunk1 kilosort2.5 8"
    exit 0
fi

EXP=$1
shift
CHUNK=$1
shift
ALG=$1
shift
NUM_CPU=$1


# Set the default sorting algorithm to kilosort2.
ALG=${ALG:-'kilosort2'}

echo "Using ${ALG} for sorting."

# Set the default number of CPUs to use if no input present.
NUM_CPU=${NUM_CPU:-6}

# Kick out the number of CPUs.
echo "Using ${NUM_CPU} CPUs"

# Set the paths by the computer.
# source data_paths.sh
# VISIONPATH=${VISIONPATH}
# if [[ "$HOSTNAME" == *"hyak"* || "$HOSTNAME" == *"klone"* ]]; then
#     SORTED_SPIKE_PATH='/gscratch/retina/data/sorted'; 
#     KILOSORT_TTL_PATH='/gscratch/retina/data/sorted';
#     RAW_DATA_PATH='/gscratch/scrubbed/retina/data/raw';
#     VISIONPATH="/gscratch/retina/GitRepos/manookin-lab/MEA/src/Vision7_for_2015DAQ/Vision.jar";    
# else
#     SORTED_SPIKE_PATH='/usr/share/pool/data/sorted'; #"/home/mike/ftp/files/ksort_out/"; #'/home/mike/ftp/files/20220420C/data021/kilosort2/';
#     KILOSORT_TTL_PATH='/usr/share/pool/data/sorted'; #'/home/mike/ftp/files/ksort_out/';
#     RAW_DATA_PATH='/usr/share/pool/data/raw';
#     VISIONPATH="/home/mike/Documents/git_repos/manookin-lab/MEA/src/Vision7_for_2015DAQ/Vision.jar";
# fi

# source data_paths.sh

TMP_PATH="${SORTED_SPIKE_PATH}/${EXP}/"

data_files=($(head -n 1 "${TMP_PATH}${EXP}_${CHUNK}.txt"))

data_string=""
for arg in "${data_files[@]}"; do
   data_string=$data_string"${DATA_PATH}$arg "
done

# Concatenate the experiment name.
EXPERIMENT_SPIKE_PATH="${SORTED_SPIKE_PATH}/${EXP}";
KILOSORT_TTL_PATH="${KILOSORT_TTL_PATH}/${EXP}/TTLTriggers";
EXPERIMENT_DATA_PATH="${RAW_DATA_PATH}/${EXP}";

CHUNK_SPIKE_PATH="${EXPERIMENT_SPIKE_PATH}/${CHUNK}/${ALG}/"
VISION_OUT="${EXPERIMENT_SPIKE_PATH}/${CHUNK}/${ALG}/"
SORT_PATH=${EXPERIMENT_SPIKE_PATH} 

# Generate the neurons file for the chunk.
if [ "${ALG}" = "yass" ]; then
    echo "Collecting spikes from yass"
    python yass_output_to_neurons.py ${CHUNK_SPIKE_PATH}tmp/output/spike_train.npy ${CHUNK_SPIKE_PATH}tmp/ ${CHUNK_SPIKE_PATH} yass
    # Split the spikes into the appropriate data files.
    python split_yass_npy.py ${CHUNK_SPIKE_PATH}tmp/output/spike_train.npy ${EXPERIMENT_SPIKE_PATH}/${CHUNK}.csv ${EXPERIMENT_SPIKE_PATH}
else
    echo "Collecting spikes from ${ALG}"
    echo ${data_files[0]}
    FILE_DATA_PATH="${EXPERIMENT_DATA_PATH}/${data_files[0]}/"
    # python kilosort_to_vision.py $CHUNK_SPIKE_PATH $KILOSORT_TTL_PATH $VISION_OUT kilosort2 -l -q -d $FILE_DATA_PATH
    python kilosort_to_vision.py $CHUNK_SPIKE_PATH $KILOSORT_TTL_PATH $VISION_OUT $ALG -l -q -d $FILE_DATA_PATH
fi

# Check whether this is an old file and the EI needs to be corrected.
EXP_NUM=$(echo $EXP | cut -c1-8)

for arg in "${data_files[@]}"; do
    FNAME=$arg
    echo "Converting file ${FNAME} to Vision format."

    FILE_SPIKE_PATH="${EXPERIMENT_SPIKE_PATH}/${FNAME}/${ALG}/"
    FILE_DATA_PATH="${EXPERIMENT_DATA_PATH}/${FNAME}/"
    VISION_OUT="${EXPERIMENT_SPIKE_PATH}/${FNAME}/"

    # Create the directories if they don't exist.
    if [[ ! -e $VISION_OUT ]]; then
        mkdir -p $VISION_OUT
    fi

    if [[ ! -e $FILE_SPIKE_PATH ]]; then
        mkdir -p $FILE_SPIKE_PATH
    fi

    if [ "${ALG}" = "yass" ]; then
        # python yass_output_to_neurons.py ${SORT_PATH}/${FNAME}.npy ${SPIKE_PATH}tmp/ ${master_dir} ${FNAME} -d ${DATA_PATH}${FNAME}
        python yass_output_to_neurons.py ${EXPERIMENT_SPIKE_PATH}/${FNAME}.npy ${EXPERIMENT_SPIKE_PATH}tmp/ ${VISION_OUT} ${FNAME} -d ${FILE_DATA_PATH}
    else
        python kilosort_to_vision.py $FILE_SPIKE_PATH $KILOSORT_TTL_PATH $VISION_OUT $FNAME -l -q -d $FILE_DATA_PATH
        # python kilosort_to_vision.py $SPIKE_PATH $KILOSORT_TTL_PATH $VISION_OUT $FNAME -l -d $DATA_PATH
    fi

    java -Xmx8G -cp $VISIONPATH edu.ucsc.neurobiology.vision.calculations.CalculationManager "Electrophysiological Imaging Fast" $VISION_OUT $FILE_DATA_PATH 0.01 67 133 1000000 ${NUM_CPU}
    java -Xmx8G -cp $VISIONPATH edu.ucsc.neurobiology.vision.calculations.CalculationManager "Raw Data Noise Evaluation" $FILE_DATA_PATH $VISION_OUT
    # Now move the Vision files into the sort directory.
    mv "${VISION_OUT}${FNAME}.ei" ${FILE_SPIKE_PATH}
    mv "${VISION_OUT}${FNAME}.globals" ${FILE_SPIKE_PATH}
    mv "${VISION_OUT}${FNAME}.neurons" ${FILE_SPIKE_PATH}
    mv "${VISION_OUT}${FNAME}.noise" ${FILE_SPIKE_PATH}
    if [ "${ALG}" = "yass" ]; then
        mv "${EXPERIMENT_SPIKE_PATH}/${FNAME}.npy" ${FILE_SPIKE_PATH}
    fi

    # Fix the EI map (Manookin/Rieke lab specific).
    if typeset -p LAB_NAME 2> /dev/null | grep -q '^'; then
        if [ "${LAB_NAME}" == "RiekeManookin" ]; then
            if (( EXP_NUM < 20230228 )); then
                python fix_electrode_map.py ${EXPERIMENT_SPIKE_PATH} -a ${ALG} -f ${FNAME}
            fi
        fi
    fi
done

echo "Done!"
