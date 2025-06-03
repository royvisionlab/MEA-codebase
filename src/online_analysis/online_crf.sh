#!/bin/bash

# USAGE: bash online_crf.sh 20221117C data015 

# Put in help message.
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: $0 <EXPERIMENT_NAME> <FILE_NAME>"
    echo "Example: bash $0 20221117C data015 "
    exit 0
fi

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: $0 <EXPERIMENT_DATE> <CHUNK_NAME> -f <DATA_FILES> -e <EI_FILES> -n <NOISE_FILES> -a <ARRAY_SPACING> -p <PROTOCOL>" 
    echo "Example: bash $0 20240926C chunk1 -f 'data000 data001 data002' -e 'data000' -n 'data000' -a 120 -p 'SpatialNoise'"
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

# Defaults.
DEFAULT_CONTRASTS="5 10 20"
DEFAULT_PRE_TIME="250"
DEFAULT_STIM_TIME="20000"
DEFAULT_FREQUENCY="4"

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

POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    -c|--contrasts)
      CONTRASTS="$2"
      shift # past argument
      shift # past value
      ;;
    -p|--pre_time)
      PRE_TIME="$2"
      shift # past argument
      shift # past value
      ;;
    -s|--stim_time)
      STIM_TIME="$2"
      shift # past argument
      shift # past value
      ;;
    -t|--frequency)
      FREQUENCY="$2"
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

# Set the positional args.
CONTRASTS="${CONTRASTS:-$DEFAULT_CONTRASTS}"
PRE_TIME="${PRE_TIME:-$DEFAULT_PRE_TIME}"
STIM_TIME="${STIM_TIME:-$DEFAULT_STIM_TIME}"
FREQUENCY="${FREQUENCY:-$DEFAULT_FREQUENCY}"
NUM_CPU="${NUM_CPU:-"8"}"

echo "EI FILES  = ${EI_FILES}"
echo "CHUNK_FILES  = ${CHUNK_FILES}"
echo "NOISE FILES  = ${NOISE_FILES}"
echo "ARRAY SPACING = ${ARRAY_SPACING}"
echo "USE CAR = ${USE_CAR}"
echo "PROTOCOL = ${PROT}"
echo "ALGORITHMS = ${ALGORITHMS[@]}"

set -- "${POSITIONAL_ARGS[@]}" # restore positional args

EXPERIMENT_NAME=${POSITIONAL_ARGS[0]}
FILE_NAME=${POSITIONAL_ARGS[1]}

# Get the start time.
START_TIME=$(date +%s)


# Set the paths by the computer.
# if [[ "$HOSTNAME" == *"hyak"* || "$HOSTNAME" == *"klone"* ]]; then
#     SORTED_SPIKE_PATH='/gscratch/retina/data/sorted'; #"/home/mike/ftp/files/ksort_out/"; #'/home/mike/ftp/files/20220420C/data021/kilosort2/';
#     RAW_DATA_PATH='/gscratch/scrubbed/retina/data/raw';
#     VISIONPATH="/gscratch/retina/GitRepos/manookin-lab/MEA/src/Vision7_for_2015DAQ/";
# else
#     SORTED_SPIKE_PATH='/data/data/sorted'; #"/home/mike/ftp/files/ksort_out/"; #'/home/mike/ftp/files/20220420C/data021/kilosort2/';
#     KILOSORT_TTL_PATH='/data/data/sorted'; #'/home/mike/ftp/files/ksort_out/';
#     RAW_DATA_PATH='/data/data/raw';
#     VISIONPATH="/home/mike/Documents/git_repos/manookin-lab/MEA/src/Vision7_for_2015DAQ/";
# fi

# Set the paths.
SORT_PATH="${SORTED_SPIKE_PATH}/${EXPERIMENT_NAME}/${FILE_NAME}/vision/"
DATA_PATH="${RAW_DATA_PATH}/${EXPERIMENT_NAME}/"

if [[ ! -e $SORT_PATH ]]; then
   mkdir -p ${SORT_PATH}
fi

# Run spike sorting in Vision.
java -d64 -Xmx8000m -Xss2m -cp ${VISIONPATH}Vision.jar edu.ucsc.neurobiology.vision.tasks.NeuronIdentification ${DATA_PATH}${FILE_NAME} ${SORT_PATH} -c ${VISIONPATH}config.xml

# Get the EI file.
# java -Xmx8G -cp ${VISIONPATH}Vision.jar edu.ucsc.neurobiology.vision.calculations.CalculationManager "Electrophysiological Imaging Fast" ${SORT_PATH}/ ${RAW_DATA_PATH}/${EXPERIMENT_NAME}/${FILE_NAME}/ 0.01 20 40 1000000 ${NUM_CPU}

# Run the basic CRF analysis
python crf_analysis.py ${SORT_PATH} ${FILE_NAME} -a vision -c ${CONTRASTS} -p ${PRE_TIME} -s ${STIM_TIME} -t ${FREQUENCY}

END_TIME=$(date +%s)
RUN_TIME=$((END_TIME-START_TIME))
MINUTES_TIME=$((RUN_TIME/60))
SECONDS_TIME=$((RUN_TIME-MINUTES_TIME*60))
echo "Completed running after ${MINUTES_TIME} minutes and ${SECONDS_TIME} seconds."