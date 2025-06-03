#!/bin/bash

# USAGE: bash analyze_chunk.sh 20220531C noise kilosort2.5 0.7 data005 data006 data010
# To use a different protocol than the default: export PROT='AdaptNoiseColorSteps'

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: $0 <EXPERIMENT_DATE>"
    echo "Example: bash $0 20240926C"
    exit 0
fi

EXP=$1

RAW_PATH="/big_data/data/raw/"
JSON_PATH="/data/data/metadata/json/"
H5_PATH="/data/data/h5/"

if [ ! -d "${RAW_PATH}${EXP}/" ]; then
    echo "${RAW_PATH}${EXP}/ does not exist. Exiting."
    exit 1
fi

if [ ! -d "${JSON_PATH}" ]; then
    echo "${JSON_PATH} does not exist. Exiting."
    exit 1
fi

if [ ! -f "${H5_PATH}${EXP}.h5" ]; then
    echo "${H5_PATH}${EXP}.h5 does not exist. Exiting."
    exit 1
fi


# EXP='20240813H'; 
echo "Parsing Symphony H5 data for ${EXP}."
python parse_data.py ${H5_PATH}${EXP}.h5 ${JSON_PATH} -r ${RAW_PATH}

echo "Moving summary TXT file to text folder."
mv ${JSON_PATH}${EXP}.txt /${JSON_PATH}txt/

echo "Moving JSON file to datajoint folder."
mv ${JSON_PATH}${EXP}_dj.json ${JSON_PATH}datajoint/${EXP}.json

