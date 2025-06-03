#!/bin/bash -l

# Usage: sbatch --export=EXP='20230105C',CHUNK='chunk1',ARRAY='504',ALG='2.5' sort_kilosort.slurm

#SBATCH --job-name=sort_kilosort        # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=manookin@uw.edu     # Where to send mail	

#SBATCH --account=retina
#SBATCH --partition=ckpt-all #gpu-a40 #ckpt gpu-a40 ckpt-all ckpt-g2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6 # 8
#SBATCH --mem=119G #256G
#SBATCH --gpus=1
#SBATCH --constraint="a100" #--constraint="a100|a40"
#SBATCH --time=1-12:00:00 # Max runtime in DD-HH:MM:SS format.

#SBATCH --chdir=/gscratch/retina/GitRepos/manookin-lab/MEA/src/utilities/
#SBATCH --export=all
#SBATCH --output=/gscratch/retina/logs/log_%j.log # where STDOUT goes
#SBATCH --error=/gscratch/retina/logs/error_%j.log # where STDERR goes

ARRAY='504'
ARRAY_ID="${ARRAY_ID:-$ARRAY}"
ALG='2.5'
SORT_ALGORITHM="${SORT_ALGORITHM:-$ALG}"

ARRAY_SPACING="60"

format_time() {
  ((h=${1}/3600))
  ((m=(${1}%3600)/60))
  ((s=${1}%60))
  printf "%02d:%02d:%02d\n" $h $m $s
 }

# Get the current date and time.
 NOW=`date +%Y-%m-%d_%H:%M:%S`

 # Log the job.
echo "Job ${SLURM_JOB_ID} started for ${EXP}, ${CHUNK} on ${NOW}" >> /gscratch/retina/logs/sort_kilosort${ALG}.txt

# Modules to use.
module load singularity
module load matlab

 # Define the sorting algorithms to run.
declare -a ALGORITHMS=("4")  # declare -a ALGORITHMS=("2.5" "4") #ALGORITHMS=("2" "2.5" "3" "4")

## now loop through the above array
for ALG in "${ALGORITHMS[@]}"
do
  if [[ "$ALG" == "4" ]]; then
  # Run the analysis in your container.
singularity shell --nv --bind /gscratch /gscratch/retina/containers/a40.sif <<EOF
source /opt/conda/etc/profile.d/conda.sh
conda activate kilosort
python run_kilosort4.py ${EXP} ${CHUNK} -e ${ARRAY_SPACING}
EOF
  else
  matlab -nodesktop -nosplash -nodisplay -r "run_kilosort('${EXP}','${CHUNK}', ${ARRAY_ID}, 'version', ${ALG}); quit"
  fi
done



# # Run the analysis in your container.
# singularity shell --nv --bind /gscratch /gscratch/retina/containers/a40.sif <<EOF
# source /opt/conda/etc/profile.d/conda.sh
# conda activate yass

# python sta_analysis.py ${EXPERIMENT_DATE} -n -c ${CHUNK_NAME} -a ${SORT_ALGORITHM} -f ${DATA_FILES} -p ${PROTOCOL_ID}

# if (($NUM_FILES > 1)); then
#     # Merge the EI files for the noise runs.
#     python merge_ei.py ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/ -c ${CHUNK_NAME} -a ${SORT_ALGORITHM} -f ${DATA_FILES}
# else
#     cp ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${DATA_FILES}/${SORT_ALGORITHM}/${DATA_FILES}.ei ${SORTED_SPIKE_PATH}/${EXPERIMENT_DATE}/${CHUNK_NAME}/${SORT_ALGORITHM}/${SORT_ALGORITHM}.ei
# fi
# EOF

# echo "Job ${SLURM_JOB_ID} finished for ${EXP}, ${CHUNK} in $(format_time $SECONDS)" >> /gscratch/retina/logs/sort_kilosort${ALG}.txt