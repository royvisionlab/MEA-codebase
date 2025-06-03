#!/bin/bash

# USAGE: bash analyze_chunk.sh 20220531C noise kilosort2.5 0.7 data005 data006 data010
# To use a different protocol than the default: export PROT='AdaptNoiseColorSteps'

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: $0 <EXPERIMENT_DATE> <CHUNK_NAME> <SORT_ALGORITHM> <CROP_FRACTION> <DATA_FILES>"
    echo "Example: bash analyze_chunk.sh 20240926C chunk1 kilosort2.5 0.6 data000"
    echo "To use a different protocol than the default (SpatialNoise): export PROT='AdaptNoiseColorSteps'"
    exit 0
fi

# Parse the input parameters.
EXPERIMENT_DATE=$1
shift
CHUNK_NAME=$1
shift
SORT_ALGORITHM=$1
shift
CROP_FRACTION=$1
shift 
DATA_FILES="$*"

DEFAULT_PROTOCOL="SpatialNoise" # "FastNoise" "AdaptNoiseColorSteps"
PROTOCOL_ID="${PROT:-$DEFAULT_PROTOCOL}"

# Get the number of data files.
files=($DATA_FILES)
NUM_FILES=${#files[@]}

# Check whether the directory exists.
# source ../../../utilities/data_paths.sh

# Check whether the directory exists.
DATA_PATH="${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${CHUNK_NAME}/${SORT_ALGORITHM}/"
if [[ ! -e $DATA_PATH ]]; then
    mkdir -p $DATA_PATH
elif [[ ! -d $DATA_PATH ]]; then
    echo "$DATA_PATH already exists but is not a directory" 1>&2
fi

python sta_analysis.py ${EXPERIMENT_DATE} -c ${CHUNK_NAME} -a ${SORT_ALGORITHM} -f ${DATA_FILES} -p ${PROTOCOL_ID} -x ${CROP_FRACTION}

if (($NUM_FILES > 1)); then
    # Merge the EI files for the noise runs.
    python merge_ei.py ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/ -c ${CHUNK_NAME} -a ${SORT_ALGORITHM} -f ${DATA_FILES}
else
    cp ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${DATA_FILES}/${SORT_ALGORITHM}/${DATA_FILES}.ei ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${CHUNK_NAME}/${SORT_ALGORITHM}/${SORT_ALGORITHM}.ei
fi

 # Instructions for running mapping analysis in Vision.
# Concatenated mapping analysis command:
# java -Xmx20G -Xss1G -classpath /Volumes/Lab/Development/vision7/Vision.jar edu.ucsc.neurobiology.vision.tasks.ConcatenatedMappingAnalysis '/Volumes/Data/2007-01-23-5/data000-data010' '/Volumes/Analysis/2007-01-23-5/d00-10-norefit' -n -c /Volumes/Lab/Development/vision-xml/current/primate.xml
# Will concatenate several data runs (data000 through data010 in this example), sort them together, assign the same unit ids (no refitting each run), and put them into /Volumes/Analysis/2007-01-23-5/d00-10-norefit folder. Specify config file if desired.

# Mapping analysis command:
# java -Xmx20G -Xss1G -classpath /Volumes/Lab/Development/vision7/Vision.jar edu.ucsc.neurobiology.vision.tasks.MappingAnalysis '/Volumes/Analysis/2016-02-17-4/data001' '/Volumes/Analysis/2016-02-17-4/data002_from_data001' '/Volumes/Data/2016-02-17-4/data002' -c /Volumes/Lab/Development/vision-xml/current/primate.xml
# Will take already sorted data001 (path to sorted data001 results - first argument), sort data002 and map it onto data001 (path to where to put the sorted data002 - second argument, path to the raw data002 to be sorted - third argument), so the unit ids will be the same. Specify config if desired
# for the mapping analysis, you'll need .model and .prj files produced by Vision for the template run