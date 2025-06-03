#!/bin/bash

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: $0 <EXPERIMENT_DATE> <FILE_NAME> <SEED> <ARRAY_ID> <PROTCOL_ID>"
    echo "Example w/ default protocol: $0 20240926C data015 1278606431 504"
    echo "Example: $0 20240926C data015 1278606431 504 FastNoise"
    exit 0
fi

DEFAULT_PROTOCOL="FastNoise"
DEFAULT_ARRAY_ID="504"

EXP=$1
shift
FILE_NAME=$1
shift
SEED=$1
shift
ARRAY_ID=$1 #"${$1:-$DEFAULT_ARRAY_ID}"
shift 
PROTOCOL_ID="${1:-$DEFAULT_PROTOCOL}"


format_time() {
  ((h=${1}/3600))
  ((m=(${1}%3600)/60))
  ((s=${1}%60))
  printf "%02d:%02d:%02d\n" $h $m $s
 }

source ../utilities/data_paths.sh
VISION_PATH=$VISIONPATH

SORT_PATH="${SORTED_SPIKE_PATH}/${EXP}/"
DATA_PATH="${RAW_DATA_PATH}/${EXP}/"

data_string="${DATA_PATH}$FILE_NAME"

if [[ ! -e $SORT_PATH ]]; then
   mkdir ${SORT_PATH}
fi

# Get the current date and time.
NOW=`date +%Y-%m-%d_%H:%M:%S`

# Log the job.
echo "Started online sorting for ${EXP}, ${FILE_NAME} on ${NOW}"

# Prepare the data.
python ../utilities/convert_litke_to_yass.py ${data_string} ${SORT_PATH} ${FILE_NAME} -w

matlab -nodesktop -nosplash -r "run_kilosort_single('${EXP}','${FILE_NAME}', ${ARRAY_ID}); quit"

# Convert to Vision format.
python ../utilities/kilosort_to_vision.py ${SORT_PATH}${FILE_NAME}/kilosort2 ${SORT_PATH}${FILE_NAME} ${SORT_PATH}${FILE_NAME} ${FILE_NAME} -l -d ${RAW_DATA_PATH}/${EXP}/${FILE_NAME}/

# Get the EI file.
echo "Computing the EI file for ${EXP}, ${FILE_NAME}."
java -Xmx8G -cp $VISION_PATH edu.ucsc.neurobiology.vision.calculations.CalculationManager "Electrophysiological Imaging Fast" ${SORT_PATH}${FILE_NAME}/ ${RAW_DATA_PATH}/${EXP}/${FILE_NAME}/ 0.01 20 40 1000000 8

# Now run deduplication on the EI files.
echo "Running deduplication on the EI files."
python ../analysis/deduplication/deduplication.py ${EXP} -a kilosort2 -c ${FILE_NAME} -o

# Run the simple STA analysis.
python sta_single_file.py ${SORT_PATH}${FILE_NAME} ${FILE_NAME} ${SEED} -p ${PROTOCOL_ID} -e ${EXP} -i ${ARRAY_ID}


echo "Online kilosort finished for ${EXP}, ${FILE_NAME} in $(format_time $SECONDS)" 