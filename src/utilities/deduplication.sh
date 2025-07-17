#!/bin/bash

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: $0 <EXPERIMENT_DATE> <CHUNK_NAME> <SORT_ALGORITHM> <NUM_CPU> <EI_FILES>"
    echo "Example: bash deduplication.sh 20240926C chunk1 kilosort2.5 8 data000"
    exit 0
fi

EXP=$1
shift
CHUNK=$1
shift
ALG=$1
shift
NUM_CPU=$1
shift 
EI_FILES="$*"

# Set the default sorting algorithm to kilosort2.
ALG=${ALG:-'kilosort2'}

# Set the default number of CPUs to use if no input present.
NUM_CPU=${NUM_CPU:-6}


# Set the paths by the computer.
# source data_paths.sh
JAR_PATH=${VISIONPATH}

TMP_PATH="${SORTED_SPIKE_PATH}${EXP}/"

data_files=($(head -n 1 "${TMP_PATH}${EXP}_${CHUNK}.txt"))

# Get the data files to combine for the EI file.
EI_FILES=${EI_FILES:-${data_files[*]}}
# Get the number of data files to combine for the EI.
files=($EI_FILES)
NUM_FILES=${#files[@]}

data_string=""
for arg in "${data_files[@]}"; do
   data_string=$data_string"${DATA_PATH}$arg "
done

# Concatenate the experiment name.
EXPERIMENT_SPIKE_PATH="${SORTED_SPIKE_PATH}${EXP}";
KILOSORT_TTL_PATH="${KILOSORT_TTL_PATH}${EXP}/TTLTriggers";
EXPERIMENT_DATA_PATH="${RAW_DATA_PATH}${EXP}";

CHUNK_SPIKE_PATH="${EXPERIMENT_SPIKE_PATH}/${ALG}/${CHUNK}/"
VISION_OUT="${EXPERIMENT_SPIKE_PATH}/${ALG}/${CHUNK}/"
SORT_PATH=${EXPERIMENT_SPIKE_PATH} 

# Check whether this is an old file and the EI needs to be corrected.
EXP_NUM=$(echo $EXP | cut -c1-8)

echo "Combining EI files for the chunk."
if (($NUM_FILES > 1)); then
    # Merge the EI files for the noise runs.
    python ../analysis/protocol/typing/merge_ei.py ${SORTED_SPIKE_PATH}${EXP}/ -c ${CHUNK} -a ${ALG} -f ${EI_FILES[*]}
else
    cp ${SORTED_SPIKE_PATH}${EXP}/${EI_FILES}/${ALG}/${EI_FILES}.ei ${SORTED_SPIKE_PATH}${EXP}/${ALG}/${CHUNK}/${CHUNK}.ei
fi

# Now run deduplication on the EI files.
echo "Running deduplication on the EI files."
python ../analysis/deduplication/deduplication.py ${EXP} -a ${ALG} -c ${CHUNK}

echo "Recomputing individual EI files."
for arg in "${data_files[@]}"; do
    FNAME=$arg
    echo "Converting file ${FNAME} to Vision format."

    FILE_SPIKE_PATH="${EXPERIMENT_SPIKE_PATH}/${ALG}/${FNAME}/"
    FILE_DATA_PATH="${EXPERIMENT_DATA_PATH}/${FNAME}/"
    #VISION_OUT="${EXPERIMENT_SPIKE_PATH}/${ALG}/${FNAME}/"

    # Move out the neurons file so that it can be used for the Vision calculation.
    #mv "${FILE_SPIKE_PATH}${FNAME}.neurons" ${VISION_OUT}

    java -Xmx8G -cp $JAR_PATH edu.ucsc.neurobiology.vision.calculations.CalculationManager "Electrophysiological Imaging Fast" $FILE_SPIKE_PATH $FILE_DATA_PATH 0.01 67 133 1000000 ${NUM_CPU}

echo "Recombining EI files for the chunk."
# Now, you need to regenerate the EI files for the chunk after deduplication.
if (($NUM_FILES > 1)); then
    # Merge the EI files for the noise runs.
    python ../analysis/protocol/typing/merge_ei.py ${SORTED_SPIKE_PATH}${EXP}/ -c ${CHUNK} -a ${ALG} -f ${EI_FILES[*]}
else
    cp ${SORTED_SPIKE_PATH}${EXP}/${EI_FILES}/${ALG}/${EI_FILES}.ei ${SORTED_SPIKE_PATH}${EXP}/${ALG}/${CHUNK}/${CHUNK}.ei
fi

echo "Done!"



