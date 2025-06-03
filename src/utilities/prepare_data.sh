#!/bin/bash

# USAGE: bash prepare_data.sh 20220531C chunk2 data005 data006 data010

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: $0 <EXPERIMENT_DATE> <CHUNK_NAME> <DATA_FILES>"
    echo "Example: bash $0 20240926C chunk1 data000 data001 data002"
    exit 0
fi

data_files=""
data_string=""

EXP=$1
shift
CHUNK=$1

if [ "$#" -gt 1 ]; then
    shift
    data_files=( "$@" )
else
   data_files=()
fi

# Set the paths by the computer.
# source data_paths.sh

TMP_PATH="${SORTED_SPIKE_PATH}${EXP}/"

SORT_PATH="${SORTED_SPIKE_PATH}${EXP}/"
DATA_PATH="${RAW_DATA_PATH}${EXP}/"

# If the last argument is 'all', run all of the files.
if ((${#data_files[@]})); then
   echo "array is not empty"
else
   echo "array is empty"
   for d in $DATA_PATH/*; do
   if [ -d "$d" ]; then
      FNAME=$(basename $d)
      if [[ $FNAME == *"data"* ]]; then
         echo ${FNAME}
         # Append to the array.
         data_files+=(${FNAME})
      fi
   fi
   done
fi

data_string=""
for arg in "${data_files[@]}"; do
   data_string=$data_string"${DATA_PATH}$arg "
done
# echo $data_string

if [[ ! -e $TMP_PATH ]]; then
   mkdir ${TMP_PATH}
fi

# Write the file names to a text file.
rm -f "${TMP_PATH}${EXP}_${CHUNK}.txt" # Remove the file if it exists.
echo ${data_files[@]} >> "${TMP_PATH}${EXP}_${CHUNK}.txt"

# Combine the raw Litke bin files into a yass file. (Also works for a single file)
python convert_join_litke_datasets_for_yass.py ${data_string} ${TMP_PATH} ${CHUNK} -b
