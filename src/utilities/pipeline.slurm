#!/bin/bash -l

# Usage: sbatch --export=EXPERIMENT_DATE='20240611H',CHUNK='chunk1',CHUNK_FILES='data001',NOISE_FILES='data001',ARRAY_SPACING='30',PROT='SpatialNoise' pipeline.slurm
# sbatch --export=EXPERIMENT_DATE='20240718H',CHUNK='chunk2',CHUNK_FILES='data000 data007 data008 data009',EI_FILES='data007',NOISE_FILES='data007',ARRAY_SPACING='120',PROT='SpatialNoise' pipeline.slurm

#SBATCH --job-name=pipeline             # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=manookin@uw.edu     # Where to send mail	
#SBATCH --account=retina
#SBATCH --partition=gpu-a40 #ckpt-all gpu-a40
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6 
#SBATCH --mem=119G 
#SBATCH --gpus=1
#SBATCH --time=2-00:00:00 # Max runtime in DD-HH:MM:SS format.
# #SBATCH --constraint="a100"
#SBATCH --chdir=/gscratch/retina/GitRepos/manookin-lab/MEA/src/utilities/
#SBATCH --export=all
#SBATCH --output=/gscratch/retina/logs/log_%j.log # where STDOUT goes
#SBATCH --error=/gscratch/retina/logs/error_%j.log # where STDERR goes


# if [[ "$1" == "-h" || "$1" == "--help" ]]; then
#     echo "Usage: sbatch --export=EXPERIMENT_DATE=<EXPERIMENT_DATE>,CHUNK=<CHUNK_NAME>,CHUNK_FILES=<DATA_FILES>,EI_FILES=<EI_FILES>,NOISE_FILES=<NOISE_FILES>,ARRAY_SPACING=<ARRAY_SPACING>,USE_CAR=<USE_CAR>,PROT=<PROTOCOL>,ALGORITHMS=<ALGORITHMS> $0"
#     echo "Example: sbatch --export=EXPERIMENT_DATE='20240611H',CHUNK='chunk1',CHUNK_FILES='data000 data007 data008 data009',EI_FILES='data007',NOISE_FILES='data007',ARRAY_SPACING='120',USE_CAR=false,PROT='SpatialNoise',ALGORITHMS='2.5 4' $0" 
#     echo "Note: files can also be specified as a range of numbers, e.g., '0-10,34-36'; comma and dash delimiters only. For example:"
#     echo "sbatch --export=EXPERIMENT_DATE='20240611H',CHUNK='chunk1',CHUNK_FILES='0,7-9',EI_FILES='7',NOISE_FILES='7',ARRAY_SPACING='120',USE_CAR=false,PROT='SpatialNoise',ALGORITHMS='2.5 4' $0" 
#     echo "Arguments: "
#     echo "  EXPERIMENT_DATE: The EI files to use for deduplication. Example: '20240611H' "
#     echo "  CHUNK: The chunk name. Example: 'chunk1' "
#     echo "  CHUNK_FILES: The chunk files to process. Example: 'data001 data002 data003' or '0-3' "
#     echo "  EI_FILES: The EI files to use for deduplication. Example: 'data001 data002 data003' or '0-3' "
#     echo "  NOISE_FILES: The noise files to use for RF calculation. Example: 'data001 data002 data003' or '0-3' "
#     echo "  ARRAY_SPACING: The array spacing to use for RF calculation. Default: 60"
#     echo "  USE_CAR: Use common average reference. Default: false"
#     echo "  PROT: The protocol to use for RF calculation. Default: 'SpatialNoise' "
#     echo "  ALGORITHMS: The sorting algorithms to use. Default: '2.5 4' "
#     exit 0
# fi

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

# Define the sorting algorithms to run.
if [ -z "${ALGORITHMS}" ]; then
  echo "The sorting algorithms are not defined. Using default algorithms."
  declare -a ALGORITHMS=("2.5" "4") #ALGORITHMS=("2" "2.5" "3" "4")
else
  # Replace commas with spaces.
  ALGORITHMS=$(echo $ALGORITHMS | tr "," " ")
  # Convert the algorithms to an array.
  ALGORITHMS=($ALGORITHMS)
fi

# Check if the array spacing is defined.
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

DEFAULT_PROTOCOL="SpatialNoise"
PROTOCOL_ID="${PROT:-$DEFAULT_PROTOCOL}"

echo "EXP = ${EXPERIMENT_DATE}"
echo "CHUNK = ${CHUNK}"
echo "NOISE PROTOCOL = ${PROTOCOL_ID}"
echo "ALGORITHMS = ${ALGORITHMS[*]}"
echo "CHUNK_FILES  = ${CHUNK_FILES[*]}"
echo "EI FILES  = ${EI_FILES[*]}"
echo "NOISE FILES  = ${NOISE_FILES[*]}"
echo "Processing data from array ID ${ARRAY_ID}."
echo "Processing data from array with ${ARRAY_SPACING} um spacing."
echo "USE CAR = ${USE_CAR}"

# Define the paths.
LITKE_PATH='/gscratch/scrubbed/retina/data/'
SORTED_SPIKE_PATH='/gscratch/retina/data/sorted'; 
RAW_DATA_PATH='/gscratch/scrubbed/retina/data/raw'; 
JAR_PATH="/gscratch/retina/GitRepos/manookin-lab/MEA/src/Vision7_for_2015DAQ/Vision.jar"; 

TMP_PATH="${LITKE_PATH}sorted/${EXPERIMENT_DATE}/"

SORT_PATH="${LITKE_PATH}sorted/${EXPERIMENT_DATE}/"
DATA_PATH="${RAW_DATA_PATH}/${EXPERIMENT_DATE}/"

SORTED_SPIKE_PATH='/gscratch/retina/data/sorted';

# Create the directories if they don't exist.
# echo "$TMP_PATH"
if [ ! -d "$TMP_PATH" ]; then
  mkdir $TMP_PATH
fi

if [ ! -d "${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}" ]; then
  mkdir ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}
fi
if [ ! -d "${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${CHUNK}" ]; then
  mkdir ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${CHUNK}
fi


# Check whether this is an old file and the EI needs to be corrected.
EXP_NUM=$(echo $EXPERIMENT_DATE | cut -c1-8)

if (( EXP_NUM < 20230228 )); then
    OLD_EXP='True'
else
    OLD_EXP='False'
fi

format_time() {
  ((h=${1}/3600))
  ((m=(${1}%3600)/60))
  ((s=${1}%60))
  printf "%02d:%02d:%02d\n" $h $m $s
 }

# Get the current date and time.
NOW=`date +%Y-%m-%d_%H:%M:%S`

# Log the job.
# echo "Job ${SLURM_JOB_ID} started for ${EXPERIMENT_DATE}, ${CHUNK} on ${NOW}" >> /gscratch/retina/logs/sort_kilosort${ALG}.txt

# Preprocess the data.
# Load singularity.
module load singularity

# Parse the H5 file if the JSON hasn't been generated yet.
if [ ! -f /gscratch/retina/data/metadata/json/${EXPERIMENT_DATE}.json ]; then
  # Make sure the H5 file exists.
  if [ ! -f /gscratch/retina/data/h5/${EXPERIMENT_DATE}.h5 ]; then
      echo "H5 file ${LITKE_PATH}h5/${EXPERIMENT_DATE}.h5 not found. Cannot continue."
      exit 1
  fi
  {
    echo "Parsing H5 file for experiment ${EXPERIMENT_DATE}."
    python ../../database/parse_data.py /gscratch/retina/data/h5/${EXPERIMENT_DATEP}.h5 /gscratch/retina/data/metadata/json/$EXPERIMENT_DATE}.json -r ${LITKE_PATH}raw/
  } || {
    echo "An error occurred in parsing the H5 file for experiment ${EXPERIMENT_DATE}."
    exit 1
  }
fi


# Prepare the data for the sorting.
if [ ! -f ${TMP_PATH}${CHUNK}.bin ]; then
    echo "Bin file not found. Running prepare_data"
    {
  data_string=""
  for arg in "${CHUNK_FILES[@]}"; do
    data_string=$data_string"${DATA_PATH}$arg "
  done
  # Write the file names to a text file.
  rm -f "${TMP_PATH}${EXPERIMENT_DATE}_${CHUNK}.txt" # Remove the file if it exists.
  echo ${CHUNK_FILES[@]} >> "${TMP_PATH}${EXPERIMENT_DATE}_${CHUNK}.txt"

  # Combine the raw Litke bin files into a binary file. (Also works for a single file)
singularity shell --bind /gscratch /gscratch/retina/containers/a40.sif <<EOF
source /opt/conda/etc/profile.d/conda.sh
conda activate yass

echo "Converting Litke data to binary format: ${data_string}"

python convert_join_litke_datasets_for_yass.py ${data_string} ${TMP_PATH} ${CHUNK} -b
EOF

  # Copy the TXT and CSV files to the sorted directory.
  cp "${TMP_PATH}${EXPERIMENT_DATE}_${CHUNK}.txt" "${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/"
  cp "${TMP_PATH}${CHUNK}.csv" "${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/"
    } || {
      echo "An error occurred in preparing the data for the spike sorter."
      exit 1
    }
fi

# Define the sorting algorithms to run.
# declare -a ALGORITHMS=("4") #ALGORITHMS=("2" "2.5" "3" "4") ALGORITHMS=("2.5" "4")

# Modules to use (optional).
module load matlab

# Loop through the algorithms.
# ALG="2.5"
for ALG in "${ALGORITHMS[@]}"
do
SORT_ALGORITHM="kilosort${ALG}"
# Check whether the directory exists.
DATA_PATH="${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${CHUNK}/${SORT_ALGORITHM}/"

if [ ! -d "${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${CHUNK}/${SORT_ALGORITHM}/" ]; then
  mkdir ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${CHUNK}/${SORT_ALGORITHM}/
fi

if [[ ! -e ${DATA_PATH} ]]; then
    mkdir ${DATA_PATH}
elif [[ ! -d ${DATA_PATH} ]]; then
    echo "${DATA_PATH} already exists but is not a directory" 1>&2
else
    echo "$DATA_PATH already exists"
fi

echo "Running spike sorting using kilosort version ${ALG}"
# Run the sorting algorithm.
{
if [[ "$ALG" == "4" ]]; then
singularity shell --nv --bind /gscratch /gscratch/retina/containers/a40.sif <<EOF
source /opt/conda/etc/profile.d/conda.sh
conda activate kilosort

python run_kilosort4.py ${EXPERIMENT_DATE} ${CHUNK} -e ${ARRAY_SPACING} --clear_cache
EOF
else
echo "Running kilosort ${ALG}"
  matlab -nodesktop -nosplash -nodisplay -r "run_kilosort('${EXPERIMENT_DATE}','${CHUNK}', ${ARRAY_ID}, 'version', ${ALG}, 'threshold', -4.0); quit"
fi
} || {
  echo "An error occurred while running spike sorting with kilosort version ${ALG}."
  exit 1
}


# Process the spikes.
{
singularity shell --bind /gscratch /gscratch/retina/containers/a40.sif <<EOF
source /opt/conda/etc/profile.d/conda.sh
conda activate yass38

# Generate the neurons file for the chunk.
if [ "${ALG}" = "yass" ]; then
    echo "Collecting spikes from yass"
    python yass_output_to_neurons.py ${DATA_PATH}tmp/output/spike_train.npy ${CHUNK_SPIKE_PATH}tmp/ ${CHUNK_SPIKE_PATH} yass
    # Split the spikes into the appropriate data files.
    python split_yass_npy.py ${DATA_PATH}tmp/output/spike_train.npy ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${CHUNK}.csv ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}
else
    echo "Collecting spikes from ${ALG}"
    python kilosort_to_vision.py ${DATA_PATH} ${DATA_PATH} ${DATA_PATH} ${SORT_ALGORITHM} -l -q -d /gscratch/scrubbed/retina/data/raw/${EXPERIMENT_DATE}/${CHUNK_FILES[0]}/
fi
EOF

for arg in "${CHUNK_FILES[@]}"; do
  FNAME=$arg
  echo "Converting file ${FNAME} to Vision format."

  FILE_SPIKE_PATH="${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${FNAME}/${SORT_ALGORITHM}/"
  FILE_DATA_PATH="/gscratch/scrubbed/retina/data/raw/${EXPERIMENT_DATE}/${FNAME}/"
  VISION_OUT="${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${FNAME}/"

singularity shell --bind /gscratch /gscratch/retina/containers/a40.sif <<EOF
source /opt/conda/etc/profile.d/conda.sh
conda activate yass38

if [ "${ALG}" = "yass" ]; then
    python yass_output_to_neurons.py ${EXPERIMENT_SPIKE_PATH}/${FNAME}.npy ${EXPERIMENT_SPIKE_PATH}tmp/ ${VISION_OUT} ${FNAME} -d ${FILE_DATA_PATH}
elif [ "${ALG}" != "4" ]; then
    python kilosort_to_vision.py $FILE_SPIKE_PATH ${DATA_PATH} $VISION_OUT $FNAME -l -q -d $FILE_DATA_PATH
else 
    python kilosort_to_vision.py $FILE_SPIKE_PATH ${DATA_PATH} $VISION_OUT $FNAME -l -q -d $FILE_DATA_PATH
fi

java -Xmx8G -cp $JAR_PATH edu.ucsc.neurobiology.vision.calculations.CalculationManager "Electrophysiological Imaging Fast" $VISION_OUT $FILE_DATA_PATH 0.01 67 133 1000000 6
java -Xmx8G -cp $JAR_PATH edu.ucsc.neurobiology.vision.calculations.CalculationManager "Raw Data Noise Evaluation" $FILE_DATA_PATH $VISION_OUT
# Now move the Vision files into the sort directory.
mv "${VISION_OUT}${FNAME}.ei" ${FILE_SPIKE_PATH}
mv "${VISION_OUT}${FNAME}.globals" ${FILE_SPIKE_PATH}
mv "${VISION_OUT}${FNAME}.neurons" ${FILE_SPIKE_PATH}
mv "${VISION_OUT}${FNAME}.noise" ${FILE_SPIKE_PATH}
if [ "${ALG}" = "yass" ]; then
    mv "${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${FNAME}.npy" ${FILE_SPIKE_PATH}
fi

# Fix the EI map.
if [ "${OLD_EXP}" = "True" ]; then
    python fix_electrode_map.py ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE} -a ${SORT_ALGORITHM} -f ${FNAME}
fi
EOF
done
} || {
  echo "An error occurred while running spike processing for kilosort version ${ALG}."
  exit 1
}

# Copy the EI files to the chunk.
{
singularity shell --nv --bind /gscratch /gscratch/retina/containers/a40.sif <<EOF
source /opt/conda/etc/profile.d/conda.sh
conda activate yass38

if (($NUM_EI_FILES > 1)); then
    # Merge the EI files for the noise runs.
    python ../analysis/protocol/typing/merge_ei.py ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/ -c ${CHUNK} -a ${SORT_ALGORITHM} -f ${EI_FILES[*]}
else
    cp ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${EI_FILES}/${SORT_ALGORITHM}/${EI_FILES}.ei ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${CHUNK}/${SORT_ALGORITHM}/${SORT_ALGORITHM}.ei
fi
EOF
} || {
  echo "An error occurred while running RF calculation for kilosort version ${ALG}."
  exit 1
}

# Deduplicate the spikes.
if [ "$RUN_DEDUPLICATION" = true ]; then
{
singularity shell --bind /gscratch /gscratch/retina/containers/a40.sif <<EOF
source /opt/conda/etc/profile.d/conda.sh
conda activate yass

echo "Running deduplication on the EI files."
python ../analysis/deduplication/deduplication.py ${EXPERIMENT_DATE} -a ${SORT_ALGORITHM} -c ${CHUNK}
EOF

for arg in "${CHUNK_FILES[@]}"; do
  FNAME=$arg
  echo "Converting file ${FNAME} to Vision format."

  FILE_SPIKE_PATH="${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${FNAME}/${SORT_ALGORITHM}/"
  FILE_DATA_PATH="/gscratch/scrubbed/retina/data/raw/${EXPERIMENT_DATE}/${FNAME}/"
  VISION_OUT="${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${FNAME}/"

singularity shell --bind /gscratch /gscratch/retina/containers/a40.sif <<EOF
source /opt/conda/etc/profile.d/conda.sh
conda activate yass

java -Xmx8G -cp $JAR_PATH edu.ucsc.neurobiology.vision.calculations.CalculationManager "Electrophysiological Imaging Fast" $VISION_OUT $FILE_DATA_PATH 0.01 67 133 1000000 6
# Now move the Vision files into the sort directory.
mv "${VISION_OUT}${FNAME}.ei" ${FILE_SPIKE_PATH}
mv "${VISION_OUT}${FNAME}.neurons" ${FILE_SPIKE_PATH}

# Fix the EI map.
if [ "${OLD_EXP}" = "True" ]; then
    python fix_electrode_map.py ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE} -a ${SORT_ALGORITHM} -f ${FNAME}
fi
EOF
done
} || {
  echo "An error occurred while running spike processing for kilosort version ${ALG}."
  exit 1
}
fi

# Get the RFs.
if (($NUM_NOISE_FILES > 0)); then
{
if [ ! -z "${NOISE_FILES}" ]; then
singularity shell --nv --bind /gscratch /gscratch/retina/containers/a40.sif <<EOF
source /opt/conda/etc/profile.d/conda.sh
conda activate yass38

# cd ../analysis/protocol/typing/

python ../analysis/protocol/typing/sta_analysis.py ${EXPERIMENT_DATE} -c ${CHUNK} -a ${ALG} -f ${NOISE_FILES[*]} -p ${PROTOCOL_ID} -x ${CROP_FRACTION}

if (($NUM_NOISE_FILES > 1)); then
    # Merge the EI files for the noise runs.
    python ../analysis/protocol/typing/merge_ei.py ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/ -c ${CHUNK} -a ${SORT_ALGORITHM} -f ${NOISE_FILES[*]}
else
    # cp ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${NOISE_FILES}/${SORT_ALGORITHM}/${NOISE_FILES}.ei ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${CHUNK}/${SORT_ALGORITHM}/${SORT_ALGORITHM}.ei
    cp ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${NOISE_FILES}/${SORT_ALGORITHM}/${NOISE_FILES}.ei ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${CHUNK}/${SORT_ALGORITHM}/
    mv ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${CHUNK}/${SORT_ALGORITHM}/${NOISE_FILES}.ei ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${CHUNK}/${SORT_ALGORITHM}/${SORT_ALGORITHM}.ei
fi
EOF
# cd ../../../utilities 
echo "Completed pipeline for kilosort version ${ALG}."
fi
} || {
  echo "An error occurred while running RF calculation for kilosort version ${ALG}."
  exit 1
}
fi
done