#!/bin/bash

# USAGE: bash pipeline.sh <experiment> <chunk> -e <ei_files> -f <chunk_files> -n <noise_files> -a <array_spacing> 
# bash pipeline.sh 20240424H chunk1 -f "data000 data002 data003 data004" -n "data000" -e "data000" -a 120

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: $0 <EXPERIMENT_DATE> <CHUNK_NAME> -f <DATA_FILES> -e <EI_FILES> -n <NOISE_FILES> -a <ARRAY_SPACING> -p <PROTOCOL> -s <ALGORITHMS> -t <NUM_CPU>" 
    echo "Example: bash $0 20240926C chunk1 -f 'data000 data001 data002' -e 'data000' -n 'data000' -a 120 -p 'SpatialNoise' -s '2.5 4' -t 8"
    echo "Note: files can also be specified as a range of numbers, e.g., '0-10,34-36'; comma and dash delimiters only"
    echo "Options: "
    echo "  -e, --ei_files: The EI files to use for deduplication."
    echo "  -f, --chunk_files: The chunk files to process."
    echo "  -n, --noise_files: The noise files to use for RF calculation."
    echo "  -a, --array_spacing: The array spacing to use for RF calculation."
    echo "  -c, --use_car: Use common average reference."
    echo "  -p, --protocol: The protocol to use for RF calculation."
    echo "  -s, --sort_algorithms: The sorting algorithms to use."
    echo "  -t, --threads: The number of threads to use for processing."
    exit 0
fi

# export 'PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True'

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

# set -e # Force exit if error occurs.
# # Catch any errors and report the line number where the error occurred.
# trap 'catch $? $LINENO' EXIT
# catch() {
#   echo "Catching Error"
#   if [ "$1" != "0" ]; then
#     # error handling goes here
#     echo "Error $1 occurred on $2"
#   fi
# }

PROT="SpatialNoise"
DEFAULT_ALGORITHM="2.5"

POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    -e|--ei_files)
      EI_FILES="$2"
      shift # past argument
      shift # past value
      ;;
    -f|--chunk_files)
      CHUNK_FILES="$2"
      shift # past argument
      shift # past value
      ;;
    -n|--noise_files)
      NOISE_FILES="$2"
      shift # past argument
      shift # past value
      ;;
    -a|--array_spacing)
      ARRAY_SPACING="$2"
      shift # past argument
      shift # past value
      ;;
    -c|--use_car)
      USE_CAR="true"
      shift # past argument
      ;;
    -p|--protocol)
      PROT="$2"
      shift # past argument
      shift # past value
      ;;
    -s|--sort_algorithms)
      ALGORITHMS="$2"
      shift # past argument
      shift # past value
      ;;
    -t|--threads)
    NUM_CPU="$2"
      shift # past argument
      shift # past value
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

# Define defaults.
DEFAULT_SPACING="60"
ARRAY_SPACING="${ARRAY_SPACING:-$DEFAULT_SPACING}"
ALGORITHMS="${ALGORITHMS:-$DEFAULT_ALGORITHM}"
USE_CAR="${USE_CAR:-"false"}"
NUM_CPU="${NUM_CPU:-"8"}"

# Convert the algorithms to an array.
ALGORITHMS=($ALGORITHMS)

# If EI files is empty set it to the noise/chunk files.
EI_FILES="${EI_FILES:-$NOISE_FILES}"
EI_FILES="${EI_FILES:-$CHUNK_FILES}"

# Check inputs for the file types.
if [[ ! "$CHUNK_FILES" == *"data"* ]]; then
  CHUNK_FILES=$(file_range_to_file_strings $CHUNK_FILES)
fi

if [[ ! "$EI_FILES" == *"data"* ]]; then
  EI_FILES=$(file_range_to_file_strings $EI_FILES)
fi

if [[ ! "$NOISE_FILES" == *"data"* ]]; then
  NOISE_FILES=$(file_range_to_file_strings $NOISE_FILES)
fi

echo "EI FILES  = ${EI_FILES}"
echo "CHUNK_FILES  = ${CHUNK_FILES}"
echo "NOISE FILES  = ${NOISE_FILES}"
echo "ARRAY SPACING = ${ARRAY_SPACING}"
echo "USE CAR = ${USE_CAR}"
echo "PROTOCOL = ${PROT}"
echo "ALGORITHMS = ${ALGORITHMS[@]}"

set -- "${POSITIONAL_ARGS[@]}" # restore positional args

EXP=${POSITIONAL_ARGS[0]}
CHUNK=${POSITIONAL_ARGS[1]}

# echo "Positional Args: ${POSITIONAL_ARGS[*]}"
echo "EXP = ${EXP}"
echo "CHUNK = ${CHUNK}"

# Make sure chunk files are defined.
if [ -z "${CHUNK_FILES}" ]; then
  echo "The chunk files must be defined in order to process the chunk! Exiting..."
  exit 1
fi

if [ -z "${NOISE_FILES}" ]; then
  echo "Noise files are not defined. Skipping RF calculation."
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

# Set the paths by the computer.
TMP_PATH="${SORTED_SPIKE_PATH}${EXP}/"

# Parse the H5 file if the JSON hasn't been generated yet.
if [ ! -z "${NOISE_FILES}" ]; then
  if [ ! -f ${LITKE_PATH}metadata/json/${EXP}.json ]; then
    # Make sure the H5 file exists.
    if [ ! -f ${LITKE_PATH}h5/${EXP}.h5 ]; then
        echo "H5 file ${LITKE_PATH}h5/${EXP}.h5 not found. Cannot continue."
        exit 1
    fi
    {
      echo "Parsing H5 file for experiment ${EXP}."
      python ../../database/parse_data.py ${LITKE_PATH}h5/${EXP}.h5 ${LITKE_PATH}metadata/json/${EXP}.json -r ${LITKE_PATH}
    } || {
      echo "An error occurred in parsing the H5 file for experiment ${EXP}."
      exit 1
    }
  fi
fi

# Prepare the data for the sorting.
if [ ! -f ${TMP_PATH}${CHUNK}.bin ]; then
    echo "Bin file not found. Running prepare_data"
    {
      bash prepare_data.sh ${EXP} ${CHUNK} ${CHUNK_FILES}
    } || {
      echo "An error occurred in preparing the data for the spike sorter."
      exit 1
    }
fi

# Check if the chunk directory exists.
if [ ! -d "${TMP_PATH}${CHUNK}" ]; then
  echo "Creating directory ${TMP_PATH}${CHUNK}"
  mkdir ${TMP_PATH}${CHUNK}
fi

## now loop through the above array
for ALG in "${ALGORITHMS[@]}"
do
    echo "Running spike sorting using kilosort version ${ALG}"
    # Check if the directory exists.
    if [ ! -d "${TMP_PATH}${CHUNK}/kilosort${ALG}" ]; then
      echo "Creating directory ${TMP_PATH}${CHUNK}/kilosort${ALG}"
      mkdir ${TMP_PATH}${CHUNK}/kilosort${ALG}
    fi
    # Run the sorting algorithm.
    {
    if [[ "$ALG" == "4" ]]; then
        # conda activate kilosort
        if [[ "${USE_CAR}" == "true" ]]; then
          python run_kilosort4.py ${EXP} ${CHUNK} -e ${ARRAY_SPACING} -c --clear_cache
        else
          python run_kilosort4.py ${EXP} ${CHUNK} -e ${ARRAY_SPACING} --clear_cache
        fi
    else
        matlab -nodesktop -nosplash -nodisplay -r "run_kilosort('${EXP}','${CHUNK}', ${ARRAY_ID}, 'version', ${ALG}, 'threshold', -4.0, 'useCAR', ${USE_CAR}); quit"
    fi
    } || {
      echo "An error occurred while running spike sorting with kilosort version ${ALG}."
      exit 1
    }

    # Process the spikes.
    {
      bash process_spikes.sh ${EXP} ${CHUNK} kilosort${ALG} ${NUM_CPU} 
    } || {
      echo "An error occurred while running spike processing for kilosort version ${ALG}."
      exit 1
    }

    # Get the RFs; this needs to be done before deduplication.
    if [ ! -z "${NOISE_FILES}" ]; then
    {
      export PROT=${PROT}
      cd ../analysis/protocol/typing/
      bash analyze_chunk.sh ${EXP} ${CHUNK} kilosort${ALG} ${CROP_FRACTION} ${NOISE_FILES} 
      cd ../../../utilities 
      echo "Completed pipeline for kilosort version ${ALG}."
    } || {
      echo "An error occurred while running RF calculation for kilosort version ${ALG}."
      exit 1
    }
    fi

    # Deduplicate the spikes.
    if [ ! -z "${EI_FILES}" ]; then
    {
      bash deduplication.sh ${EXP} ${CHUNK} kilosort${ALG} ${NUM_CPU} ${EI_FILES}
    } || {
      echo "An error occurred while running deduplication for kilosort version ${ALG}."
      exit 1
    }
    fi

    # Get the RFs.
    if [ ! -z "${NOISE_FILES}" ]; then
    {
      export PROT=${PROT}
      cd ../analysis/protocol/typing/
      bash analyze_chunk.sh ${EXP} ${CHUNK} kilosort${ALG} ${CROP_FRACTION} ${NOISE_FILES} 
      cd ../../../utilities 
      echo "Completed pipeline for kilosort version ${ALG}."
    } || {
      echo "An error occurred while running RF calculation for kilosort version ${ALG}."
      exit 1
    }
    fi
done