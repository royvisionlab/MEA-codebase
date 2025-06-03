#!/bin/bash

# USAGE: bash litke_to_yass.sh 20220531C noise data005 data006 data010

data_files=""
litke_bin_path=""
dest_analysis_path=""
data_string=""

# litke_bin_path=$1
# shiftls
EXP=$1
shift
CHUNK=$1
shift

# data_files=( "$@" )

# Set the paths by the computer.
# if [[ "$HOSTNAME" == *"hyak"* || "$HOSTNAME" == *"klone"* ]]; then
#     YASS_SPIKE_PATH='/gscratch/retina/data/sorted'; #'/gscratch/scrubbed/retina/data/sorted';#"/home/mike/ftp/files/ksort_out/"; #'/home/mike/ftp/files/20220420C/data021/kilosort2/';
#     RAW_DATA_PATH='/gscratch/scrubbed/retina/data/raw'; #'/gscratch/scrubbed/retina/data/raw';
#     VISIONPATH="/gscratch/retina/GitRepos/manookin-lab/MEA/src/Vision7_for_2015DAQ/Vision.jar";
#     TMP_PATH="/gscratch/retina/data/sorted/${EXP}/" #"/gscratch/scrubbed/retina/data/sorted/${EXP}/"
# else
#     YASS_SPIKE_PATH='/usr/share/pool/SortedData'; #"/home/mike/ftp/files/ksort_out/"; #'/home/mike/ftp/files/20220420C/data021/kilosort2/';
#     KILOSORT_TTL_PATH='/usr/share/pool/SortedData'; #'/home/mike/ftp/files/ksort_out/';
#     RAW_DATA_PATH='/usr/share/pool/rawdata/MEA_Data';
#     VISIONPATH="/home/mike/Documents/Vision7_for_2015DAQ/Vision.jar" #"/home/mike/Documents/GitRepos/Manookin-Lab/MEA/src/Vision7_for_2015DAQ/Vision.jar";
#     TMP_PATH="/usr/share/pool/SortedData/${EXP}/" #'/media/mike/NVME/sort/'
# fi


source data_paths.sh

YASS_SPIKE_PATH="${SORTED_SPIKE_PATH}"

echo $YASS_SPIKE_PATH

TMP_PATH="${YASS_SPIKE_PATH}/${EXP}/"

SORT_PATH="${YASS_SPIKE_PATH}/${EXP}/"
DATA_PATH="${RAW_DATA_PATH}/${EXP}/"

data_files=($(head -n 1 "${TMP_PATH}${EXP}_${CHUNK}.txt"))

data_string=""
for arg in "${data_files[@]}"; do
   data_string=$data_string"${DATA_PATH}$arg "
done

# Step 4: Convert the YASS spikes file to Vision
master_dir="${SORT_PATH}${CHUNK}"
if [[ ! -e $master_dir ]]; then
    mkdir $master_dir
elif [[ ! -d $master_dir ]]; then
    echo "$master_dir already exists but is not a directory" 1>&2
fi

python yass_output_to_neurons.py ${TMP_PATH}${CHUNK}/yass/tmp/output/spike_train.npy ${TMP_PATH}${CHUNK}/yass/tmp/ ${TMP_PATH}${CHUNK}/yass yass
# python yass_output_to_neurons.py ${TMP_PATH}tmp/output/spike_train.npy ${TMP_PATH}tmp/ ${TMP_PATH} yass
# python yass_output_to_neurons.py ${TMP_PATH}tmp/output/spike_train.npy ${TMP_PATH}tmp/ ${master_dir} ${CHUNK} 

# Split the spikes into the appropriate data files.
python split_yass_npy.py ${TMP_PATH}${CHUNK}/yass/tmp/output/spike_train.npy ${TMP_PATH}${CHUNK}.csv ${TMP_PATH}

for arg in "${data_files[@]}"; do
   FNAME=$arg
   master_dir="${SORT_PATH}${FNAME}"
   if [[ ! -e $master_dir ]]; then
      mkdir $master_dir
   fi

   if [[ ! -e ${master_dir}/yass ]]; then
      mkdir ${master_dir}/yass
   fi
   echo $master_dir
   # echo ${master_dir}/yass 

   python yass_output_to_neurons.py ${TMP_PATH}${FNAME}.npy ${TMP_PATH}tmp/ ${master_dir} ${FNAME} -d ${DATA_PATH}${FNAME}

   java -Xmx8G -cp $VISIONPATH edu.ucsc.neurobiology.vision.calculations.CalculationManager "Electrophysiological Imaging Fast" $master_dir ${DATA_PATH}${FNAME} 0.01 67 133 1000000 6
   java -Xmx8G -cp $VISIONPATH edu.ucsc.neurobiology.vision.calculations.CalculationManager "Raw Data Noise Evaluation" ${DATA_PATH}${FNAME} $master_dir
   mv ${master_dir}/${FNAME}.ei ${master_dir}/yass
   mv ${master_dir}/${FNAME}.globals ${master_dir}/yass
   mv ${master_dir}/${FNAME}.neurons ${master_dir}/yass
   mv ${master_dir}/${FNAME}.noise ${master_dir}/yass
   mv ${TMP_PATH}${FNAME}.npy ${master_dir}/yass

   # Fix the EI map.
   python fix_electrode_map.py ${SORT_PATH} -a yass -f ${FNAME}
done
