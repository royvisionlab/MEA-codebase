#!/bin/bash -l

# Usage: sbatch --export=EXP='20230105C',CHUNK='chunk1' sort_yass.slurm

#SBATCH --job-name=sort_yass        # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=manookin@uw.edu     # Where to send mail	

#SBATCH --account=retina
#SBATCH --partition=gpu-a40 # ckpt-all or gpu-a40
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=119G
#SBATCH --gpus=1
#SBATCH --constraint="a40"
#SBATCH --time=6-00:00:00 # Max runtime in DD-HH:MM:SS format.

#SBATCH --chdir=/gscratch/retina/GitRepos/manookin-lab/MEA/src/utilities/
#SBATCH --export=all
#SBATCH --output=/gscratch/retina/logs/log_%j.log # where STDOUT goes
#SBATCH --error=/gscratch/retina/logs/error_%j.log # where STDERR goes

format_time() {
  ((h=${1}/3600))
  ((m=(${1}%3600)/60))
  ((s=${1}%60))
  printf "%02d:%02d:%02d\n" $h $m $s
 }

# Get the current date and time.
 NOW=`date +%Y-%m-%d_%H:%M:%S`

# Define the paths.
LITKE_PATH='/gscratch/scrubbed/retina/data/sorted'
SORTED_SPIKE_PATH='/gscratch/retina/data/sorted'; 
RAW_DATA_PATH='/gscratch/scrubbed/retina/data/raw'; 

# Check whether the directory exists.
DATA_PATH="${LITKE_PATH}/${EXP}/${CHUNK}"
if [[ ! -e $DATA_PATH ]]; then
    mkdir -p $DATA_PATH
elif [[ ! -d $DATA_PATH ]]; then
    echo "$DATA_PATH already exists but is not a directory" 1>&2
else
    echo "$DATA_PATH already exists"
fi

# Check whether the bin file exists.
if [[ ! -f "${DATA_PATH}/${CHUNK}.bin" ]]; then
    if [ -f "${LITKE_PATH}/${EXP}/${CHUNK}.bin" ]; then
        cp "${LITKE_PATH}/${EXP}/${CHUNK}.bin" "${DATA_PATH}/${CHUNK}.bin"
    else
        { echo "Error: File ${DATA_PATH}/${CHUNK}.bin does not exist."; exit 1; }
    fi
fi

# Change directory to the data path.
cd ${DATA_PATH}

# Log the start time.
echo "Job ${SLURM_JOB_ID} started for ${EXP}, ${CHUNK} on ${NOW}" >> /gscratch/retina/logs/sort_yass.txt

# Load singularity.
module load singularity

# Run the analysis in your container.
singularity shell --nv --bind /gscratch /gscratch/retina/containers/a40.sif <<EOF
source /opt/conda/etc/profile.d/conda.sh
conda activate yass38

yass sort config.yaml
EOF

# If the run completed successfully, copy the sorted spike files to the sorted spike directory.
# if [[ -f "$DATA_PATH/tmp/output/spike_train.npy" ]]; then
#     # Remove the temporary sorting files.
#     rm -rf ${DATA_PATH}/tmp/pre*
#     rm -rf ${DATA_PATH}/tmp/final*
#     rm -rf ${DATA_PATH}/tmp/block*
#     rm ${DATA_PATH}/tmp/yass.log

#     # Check whether the sorted spike directory exists.
#     SORTED_PATH="${SORTED_SPIKE_PATH}/${EXP}/${CHUNK}/yass"
#     if [[ ! -e $SORTED_PATH ]]; then
#         mkdir -p $SORTED_PATH
#     elif [[ ! -d $SORTED_PATH ]]; then
#         echo "$SORTED_PATH already exists but is not a directory" 1>&2
#     else
#         echo "$SORTED_PATH already exists"
#     fi

#     # Copy the sorted spike files.
#     cp -r ${DATA_PATH}/tmp/ ${SORTED_PATH}
# else
#     echo "Warning: The file ${DATA_PATH}/tmp/output/spike_train.npy does not exist! The run failed." 
# fi

# Log the end time.
echo "Job ${SLURM_JOB_ID} finished for ${EXP}, ${CHUNK} in $(format_time $SECONDS)" >> /gscratch/retina/logs/sort_yass.txt