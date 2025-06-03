#!/bin/bash

EXP=$1
shift
CHUNK=$1
shift
NUM_CPU=$1

# Set the default if no input present
NUM_CPU=${NUM_CPU:-6}

# Kick out the number of CPUs.
echo "Using ${NUM_CPU} CPUs"

# Set the paths by the computer.
if [[ "$HOSTNAME" == *"hyak"* || "$HOSTNAME" == *"klone"* ]]; then
    KILOSORT_SPIKE_PATH='/gscratch/retina/data/sorted'; #"/home/mike/ftp/files/ksort_out/"; #'/home/mike/ftp/files/20220420C/data021/kilosort2/'; 
    KILOSORT_TTL_PATH='/gscratch/retina/data/sorted';
    RAW_DATA_PATH='/gscratch/scrubbed/retina/data/raw';
    VISIONPATH="/gscratch/retina/GitRepos/manookin-lab/MEA/src/Vision7_for_2015DAQ/Vision.jar";
    TMP_PATH="/gscratch/retina/data/sorted/${EXP}/"
else
    KILOSORT_SPIKE_PATH='/usr/share/pool/data/sorted'; #"/home/mike/ftp/files/ksort_out/"; #'/home/mike/ftp/files/20220420C/data021/kilosort2/';
    KILOSORT_TTL_PATH='/usr/share/pool/data/sorted'; #'/home/mike/ftp/files/ksort_out/';
    RAW_DATA_PATH='/usr/share/pool/data/raw';
    VISIONPATH="/home/mike/Documents/git_repos/manookin-lab/MEA/src/Vision7_for_2015DAQ/Vision.jar";
    TMP_PATH="/usr/share/pool/data/sorted/${EXP}/" #'/media/mike/NVME/sort/'
fi

DATA_PATH="${RAW_DATA_PATH}/${EXP}/"
data_files=($(head -n 1 "${TMP_PATH}${EXP}_${CHUNK}.txt"))

data_string=""
for arg in "${data_files[@]}"; do
   data_string=$data_string"${DATA_PATH}$arg "
done

# Concatenate the experiment name.
KILOSORT_SPIKE_PATH="${KILOSORT_SPIKE_PATH}/${EXP}";
KILOSORT_TTL_PATH="${KILOSORT_TTL_PATH}/${EXP}/TTLTriggers";
RAW_DATA_PATH="${RAW_DATA_PATH}/${EXP}";


DATA_PATH="${RAW_DATA_PATH}/"
SPIKE_PATH="${KILOSORT_SPIKE_PATH}/${CHUNK}/kilosort2/"
VISION_OUT="${KILOSORT_SPIKE_PATH}/${CHUNK}/kilosort2/"
SORT_PATH=${KILOSORT_SPIKE_PATH} #"${KILOSORT_SPIKE_PATH}/${EXP}/"
# python kilosort_to_vision.py $SPIKE_PATH $KILOSORT_TTL_PATH $VISION_OUT $CHUNK -l
python kilosort_to_vision.py $SPIKE_PATH $KILOSORT_TTL_PATH $VISION_OUT kilosort2 -l

for arg in "${data_files[@]}"; do
  FNAME=$arg
  echo "Analyzing file ${FNAME}"
  master_dir="${SORT_PATH}/${FNAME}"
  if [[ ! -e $master_dir ]]; then
    mkdir -p $master_dir
  fi

  if [[ ! -e ${master_dir}/kilosort2 ]]; then
    mkdir -p ${master_dir}/kilosort2
  fi
  SPIKE_PATH="${KILOSORT_SPIKE_PATH}/${FNAME}/kilosort2/"
  DATA_PATH="${RAW_DATA_PATH}/${FNAME}/"
  VISION_OUT="${KILOSORT_SPIKE_PATH}/${FNAME}/"
  python kilosort_to_vision.py $SPIKE_PATH $KILOSORT_TTL_PATH $VISION_OUT $FNAME -l -d $DATA_PATH

  java -Xmx8G -cp $VISIONPATH edu.ucsc.neurobiology.vision.calculations.CalculationManager "Electrophysiological Imaging Fast" $VISION_OUT $DATA_PATH 0.01 67 133 1000000 ${NUM_CPU}
  java -Xmx8G -cp $VISIONPATH edu.ucsc.neurobiology.vision.calculations.CalculationManager "Raw Data Noise Evaluation" $DATA_PATH $VISION_OUT
  # Now move the Vision files into the sort directory.
  mv "${VISION_OUT}${FNAME}.ei" ${SPIKE_PATH}
  mv "${VISION_OUT}${FNAME}.globals" ${SPIKE_PATH}
  mv "${VISION_OUT}${FNAME}.neurons" ${SPIKE_PATH}
  mv "${VISION_OUT}${FNAME}.noise" ${SPIKE_PATH}

  # Fix the EI map.
  python fix_electrode_map.py ${KILOSORT_SPIKE_PATH} -a kilosort2 -f ${FNAME}
done

