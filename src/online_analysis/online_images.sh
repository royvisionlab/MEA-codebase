#!/bin/bash


EXP=$1
shift
CHUNK=$1
shift
SEED=$1
shift
IMAGE_SEED=$1
shift
NUM_IMAGES=$1

if [ "$#" -gt 1 ]; then
    shift
    data_files=( "$@" )
else
   data_files=()
fi



# The noise file is the first file in the list.
NOISE_FILE=${data_files[0]}
IMAGE_FILE=${data_files[1]}

PROTOCOL_ID="FastNoise"
ALG='kilosort2'


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

data_string=""
for arg in "${data_files[@]}"; do
   data_string=$data_string"${DATA_PATH}$arg "
done
echo $data_string

if [[ ! -e $SORT_PATH ]]; then
   mkdir ${SORT_PATH}
fi

# Write the file names to a text file.
echo ${data_files[@]} >> "${SORT_PATH}${EXP}_${CHUNK}.txt"

# Get the current date and time.
NOW=`date +%Y-%m-%d_%H:%M:%S`

# Log the job.
echo "Started online sorting for ${EXP}, ${CHUNK} on ${NOW}"

# Combine the raw Litke bin files into a yass file. (Also works for a single file)
python ../utilities/convert_join_litke_datasets_for_yass.py ${data_string} ${SORT_PATH} ${CHUNK} -b

cd ../utilities
matlab -nodesktop -nosplash -nodisplay -r "run_kilosort('${EXP}','${CHUNK}'); quit"
cd ../online_analysis

for arg in "${data_files[@]}"; do
    FNAME=$arg
    echo "Converting file ${FNAME} to Vision format."

    FILE_SPIKE_PATH="${SORT_PATH}/${FNAME}/${ALG}/"
    FILE_DATA_PATH="${DATA_PATH}/${FNAME}/"
    VISION_OUT="${SORT_PATH}/${FNAME}/"

    # Create the directories if they don't exist.
    if [[ ! -e $VISION_OUT ]]; then
        mkdir -p $VISION_OUT
    fi

    if [[ ! -e $FILE_SPIKE_PATH ]]; then
        mkdir -p $FILE_SPIKE_PATH
    fi

    python ../utilities/kilosort_to_vision.py $FILE_SPIKE_PATH $KILOSORT_TTL_PATH $VISION_OUT $FNAME -l -q -d $FILE_DATA_PATH

    # java -Xmx8G -cp $JAR_PATH edu.ucsc.neurobiology.vision.calculations.CalculationManager "Electrophysiological Imaging Fast" $VISION_OUT $FILE_DATA_PATH 0.01 20 40 1000000 8
    # java -Xmx8G -cp $JAR_PATH edu.ucsc.neurobiology.vision.calculations.CalculationManager "Raw Data Noise Evaluation" $FILE_DATA_PATH $VISION_OUT
    # Now move the Vision files into the sort directory.
    # mv "${VISION_OUT}${FNAME}.ei" ${FILE_SPIKE_PATH}
    # mv "${VISION_OUT}${FNAME}.globals" ${FILE_SPIKE_PATH}
    # mv "${VISION_OUT}${FNAME}.neurons" ${FILE_SPIKE_PATH}
    # mv "${VISION_OUT}${FNAME}.noise" ${FILE_SPIKE_PATH}
done

# Convert to Vision format.
# python ../utilities/kilosort_to_vision.py ${SORT_PATH}${CHUNK}/kilosort2 ${SORT_PATH}${CHUNK} ${SORT_PATH}${CHUNK} ${CHUNK} -l -d ${RAW_DATA_PATH}/${EXP}/${CHUNK}/

# Run the simple STA analysis.
python sta_single_file.py ${SORT_PATH}${NOISE_FILE} ${NOISE_FILE} ${SEED} -p ${PROTOCOL_ID}

python images_single_file.py ${SORT_PATH}${IMAGE_FILE} -f ${IMAGE_FILE} -s ${IMAGE_SEED} -n ${NUM_IMAGES}

# Get the EI file.
java -Xmx8G -cp $VISION_PATH edu.ucsc.neurobiology.vision.calculations.CalculationManager "Electrophysiological Imaging Fast" ${SORT_PATH}${FILE_NAME}/ ${RAW_DATA_PATH}/${EXP}/${FILE_NAME}/ 0.01 20 40 1000000 8

echo "Online kilosort finished for ${EXP}, ${CHUNK} in $(format_time $SECONDS)" 

