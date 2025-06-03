#!/bin/bash -l


file_range_to_file_strings() {
  result='' # Initialize the output string
  # Split the string into parts
  IFS=',' read -ra parts <<< "$1"

  for part in "${parts[@]}"; do
      # Check if the part contains a range
      if [[ $part == *-* ]]; then
          # Extract the start and end of the range
          IFS='-' read start end <<< "$part"
          # Loop through the range and print each number
          for (( i=$start; i<=$end; i++ )); do
              if [ "$i" -gt 99 ]; then
                result="${result} data${i}"
              elif [ "$i" -gt 9 ]; then
                result="${result} data0${i}"
              else
                result="${result} data00${i}"
              fi
          done
      else
        if [ "$part" -gt 99 ]; then
          result="${result} data${part}"
        elif [ "$part" -gt 9 ]; then
          result="${result} data0${part}"
        else
          result="${result} data00${part}"
        fi
      fi
  done
  echo ${result}
}

# Define defaults.
DEFAULT_SPACING="60"
ARRAY_SPACING="${ARRAY_SPACING:-$DEFAULT_SPACING}"

USE_CAR="${USE_CAR:-"false"}"

# Check if the chunk files are defined.
if [ -z "${CHUNK_FILES}" ]; then
  echo "The chunk files must be defined in order to process the chunk! Exiting..."
  exit 1
else
  # Check inputs for the file types.
  if [[ ! "$CHUNK_FILES" == *"data"* ]]; then
    CHUNK_FILES=$(file_range_to_file_strings $CHUNK_FILES)
  else
    # Replace commas with spaces.
    CHUNK_FILES=$(echo $CHUNK_FILES | tr "," " ")
  fi
  CHUNK_FILES=($CHUNK_FILES)
  NUM_FILES=${#CHUNK_FILES[@]}
fi

# Check if the noise files are defined.
if [ -z "${NOISE_FILES}" ]; then
  echo "The noise files are not defined. Skipping RF calculation."
  declare -a NOISE_FILES=()
  NUM_NOISE_FILES=0
else
  # Check if the noise files are defined by 'data0*' or by file numbers.
  if [[ ! "$NOISE_FILES" == *"data"* ]]; then
    NOISE_FILES=$(file_range_to_file_strings $NOISE_FILES)
  else
    # Replace commas with spaces.
    NOISE_FILES=$(echo $NOISE_FILES | tr "," " ")
  fi
  # Convert the noise files to an array.
  NOISE_FILES=($NOISE_FILES)
  # Get the number of noise files.
  NUM_NOISE_FILES=${#NOISE_FILES[@]}
fi

# Check if the EI files are defined.
RUN_DEDUPLICATION=true
if [ -z "${EI_FILES}" ]; then
  # Check if the number of noise files is greater than 0. If so, set the EI files to the noise files.
  if (($NUM_NOISE_FILES > 0)); then
    EI_FILES="${NOISE_FILES}"
  else
    echo "EI files not defined. Skipping EI deduplication."
    RUN_DEDUPLICATION=false
    EI_FILES="${CHUNK_FILES}"
  fi
else
  # Check if the EI files are defined by 'data0*' or by file numbers.
  if [[ ! "$EI_FILES" == *"data"* ]]; then
    EI_FILES=$(file_range_to_file_strings $EI_FILES)
  else
    # Replace commas with spaces.
    EI_FILES=$(echo $EI_FILES | tr "," " ")
  fi
fi

# Convert the EI files to an array.
EI_FILES=($EI_FILES)
# Get the number of EI files.
NUM_EI_FILES=${#EI_FILES[@]}

echo "EI FILES  = ${EI_FILES[*]}"
echo "CHUNK_FILES  = ${CHUNK_FILES[*]}"
echo "NOISE FILES  = ${NOISE_FILES[*]}"
echo "ARRAY SPACING = ${ARRAY_SPACING}"
echo "USE CAR = ${USE_CAR}"

echo "Positional Args: ${POSITIONAL_ARGS[*]}"
echo "EXP = ${EXPERIMENT_DATE}"
echo "CHUNK = ${CHUNK}"

# Make sure chunk files are defined.
if [ -z "${CHUNK_FILES}" ]; then
  echo "The chunk files must be defined in order to process the chunk! Exiting..."
  exit 1
fi
if [[ "$ARRAY_SPACING" == "30" ]]; then
    ARRAY_ID='1501'
    CROP_FRACTION='0.5'
elif [[ "$ARRAY_SPACING" == "120" ]]; then
    ARRAY_ID='3501'
    CROP_FRACTION='1.0'
else
    ARRAY_ID='504'
    CROP_FRACTION='0.7'
fi
echo "Processing data from array ID ${ARRAY_ID}."
echo "Processing data from array with ${ARRAY_SPACING} um spacing."

# Define the sorting algorithms to run.
declare -a ALGORITHMS=("2.5" "4") #ALGORITHMS=("2" "2.5" "3" "4")


ARRAY='504'
ARRAY_ID="${ARRAY_ID:-$ARRAY}"


DEFAULT_PROTOCOL="SpatialNoise"
PROTOCOL_ID="${PROT:-$DEFAULT_PROTOCOL}"

LITKE_PATH='/gscratch/scrubbed/retina/data/'
SORTED_SPIKE_PATH='/gscratch/retina/data/sorted'; 
RAW_DATA_PATH='/gscratch/scrubbed/retina/data/raw'; 
JAR_PATH="/gscratch/retina/GitRepos/manookin-lab/MEA/src/Vision7_for_2015DAQ/Vision.jar"; 

TMP_PATH="${LITKE_PATH}sorted/${EXPERIMENT_DATE}/"

SORT_PATH="${LITKE_PATH}sorted/${EXPERIMENT_DATE}/"
DATA_PATH="${RAW_DATA_PATH}/${EXPERIMENT_DATE}/"

SORTED_SPIKE_PATH='/gscratch/retina/data/sorted';

data_string=""
for arg in "${CHUNK_FILES[@]}"; do
data_string=$data_string"${DATA_PATH}$arg "
done

echo "Data string: ${data_string}"

# Write the file names to a text file.
# rm -f "${TMP_PATH}${EXPERIMENT_DATE}_${CHUNK}.txt" # Remove the file if it exists.
rm -f "${CHUNK}.txt" # Remove the file if it exists.
echo ${CHUNK_FILES[@]} >> "${CHUNK}.txt"

echo "Chunk files: ${CHUNK_FILES[@]}"

for arg in "${CHUNK_FILES[@]}"; do
  FNAME=$arg
  echo "Converting file ${FNAME} to Vision format."
done

if [ ! -z "${NOISE_FILES}" ]; then
    echo "Noise files: ${NOISE_FILES[@]}"
fi




