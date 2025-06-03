#!/bin/bash -l

# Usage: sbatch --export=EXPERIMENT_DATE='20230109C',CHUNK_NAME='chunk1',SORT_ALGORITHM='kilosort2.5',CROP='0.6',DATA_FILES='data001',PROT='SpatialNoise' analyze_chunk.slurm
#   sbatch --export=EXPERIMENT_DATE='20230109C',CHUNK_NAME='chunk1',DATA_FILES='data001' analyze_chunk.slurm

#SBATCH --job-name=analyze_chunk        # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=manookin@uw.edu     # Where to send mail	

#SBATCH --account=retina
#SBATCH --partition=ckpt-all #gpu-a40 #ckpt-all
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6 #8
#SBATCH --mem=128G #119G #128G
#SBATCH --gpus=1
#SBATCH --constraint="a100" # #SBATCH --constraint="a40|a100"
#SBATCH --time=0-10:00:00 # Max runtime in DD-HH:MM:SS format.

#SBATCH --chdir=/gscratch/retina/GitRepos/manookin-lab/MEA/src/analysis/protocol/typing/
#SBATCH --export=all
#SBATCH --output=/gscratch/retina/logs/log_%j.log # where STDOUT goes
#SBATCH --error=/gscratch/retina/logs/error_%j.log # where STDERR goes

DEFAULT_ALGORITHM="kilosort2"
SORT_ALGORITHM="${SORT_ALGORITHM:-$DEFAULT_ALGORITHM}"

DEFAULT_PROTOCOL="SpatialNoise"
PROTOCOL_ID="${PROT:-$DEFAULT_PROTOCOL}"

CROP_FRACTION="${CROP:-1.0}"

SORTED_SPIKE_PATH='/gscratch/retina/data/sorted';
# Check whether the directory exists.
# srun ../../../utilities/data_paths.sh

# Check whether the directory exists.
DATA_PATH="${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${CHUNK_NAME}/${SORT_ALGORITHM}/"
if [[ ! -e $DATA_PATH ]]; then
    mkdir -p $DATA_PATH
elif [[ ! -d $DATA_PATH ]]; then
    echo "$DATA_PATH already exists but is not a directory" 1>&2
else
    echo "$DATA_PATH already exists"
fi

files=($DATA_FILES)
NUM_FILES=${#files[@]}

# Load singularity.
module load singularity

# Run the analysis in your container.
singularity shell --nv --bind /gscratch /gscratch/retina/containers/a40.sif <<EOF
source /opt/conda/etc/profile.d/conda.sh
conda activate kilosort

python sta_analysis.py ${EXPERIMENT_DATE} -c ${CHUNK_NAME} -a ${SORT_ALGORITHM} -f ${DATA_FILES} -p ${PROTOCOL_ID} -x ${CROP_FRACTION}

if (($NUM_FILES > 1)); then
    # Merge the EI files for the noise runs.
    python merge_ei.py ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/ -c ${CHUNK_NAME} -a ${SORT_ALGORITHM} -f ${DATA_FILES}
else
    cp ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${DATA_FILES}/${SORT_ALGORITHM}/${DATA_FILES}.ei ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${CHUNK_NAME}/${SORT_ALGORITHM}/${SORT_ALGORITHM}.ei
fi
EOF
