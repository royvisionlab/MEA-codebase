#!/bin/bash -l

# Usage: sbatch --export=EXP='20230105C',CHUNK='chunk1',FILES='data000 data001' prepare_data.slurm

#SBATCH --job-name=prepare_data        # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=manookin@uw.edu     # Where to send mail	

#SBATCH --account=retina
#SBATCH --partition=gpu-a40 #ckpt-all
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=119G
#SBATCH --gpus=0
# #SBATCH --constraint="mem"
#SBATCH --time=0-12:00:00 # Max runtime in DD-HH:MM:SS format.

#SBATCH --chdir=/gscratch/retina/GitRepos/manookin-lab/MEA/src/utilities/
#SBATCH --export=all
#SBATCH --output=/gscratch/retina/logs/log_%j.log # where STDOUT goes
#SBATCH --error=/gscratch/retina/logs/error_%j.log # where STDERR goes


DATA_FILES=($FILES)

LITKE_PATH='/gscratch/scrubbed/retina/data/'
SORTED_SPIKE_PATH='/gscratch/retina/data/sorted'; 
RAW_DATA_PATH='/gscratch/scrubbed/retina/data/raw'; 

TMP_PATH="${LITKE_PATH}sorted/${EXP}/"

SORT_PATH="${LITKE_PATH}sorted/${EXP}/"
DATA_PATH="${RAW_DATA_PATH}/${EXP}/"

data_string=""
for arg in "${DATA_FILES[@]}"; do
   data_string=$data_string"${DATA_PATH}$arg "
done
# echo $data_string

if [[ ! -e $TMP_PATH ]]; then
   mkdir ${TMP_PATH}
fi

# Write the file names to a text file.
rm -f "${TMP_PATH}${EXP}_${CHUNK}.txt" # Remove the file if it exists.
echo ${DATA_FILES[@]} >> "${TMP_PATH}${EXP}_${CHUNK}.txt"

# Load singularity.
module load singularity

# Combine the raw Litke bin files into a binary file. (Also works for a single file)
singularity shell --bind /gscratch /gscratch/retina/containers/a40.sif <<EOF
source /opt/conda/etc/profile.d/conda.sh
conda activate yass

python convert_join_litke_datasets_for_yass.py ${data_string} ${TMP_PATH} ${CHUNK} -b
EOF

if [[ ! -e "${SORTED_SPIKE_PATH}/${EXP}/" ]]; then
   mkdir "${SORTED_SPIKE_PATH}/${EXP}/"
fi

# Copy the TXT and CSV files to the sorted directory.
cp "${TMP_PATH}${EXP}_${CHUNK}.txt" "${SORTED_SPIKE_PATH}/${EXP}/"
cp "${TMP_PATH}${CHUNK}.csv" "${SORTED_SPIKE_PATH}/${EXP}/"
