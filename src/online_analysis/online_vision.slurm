#!/bin/bash -l

# Usage: sbatch --export=EXP='20230105C',FILE_NAME='data001',SEED='1113434' online_vision.slurm

#SBATCH --job-name=online_vision        # Job name
#SBATCH --mail-type=ALL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=manookin@uw.edu     # Where to send mail	

#SBATCH --account=retina
#SBATCH --partition=ckpt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=128G
#SBATCH --gpus=1
#SBATCH --constraint="a100|a40"
#SBATCH --time=0-04:00:00 # Max runtime in DD-HH:MM:SS format.

#SBATCH --chdir=/gscratch/retina/GitRepos/manookin-lab/MEA/src/online_analysis/
#SBATCH --export=all
#SBATCH --output=/gscratch/retina/logs/log_%j.log # where STDOUT goes
#SBATCH --error=/gscratch/retina/logs/error_%j.log # where STDERR goes

format_time() {
  ((h=${1}/3600))
  ((m=(${1}%3600)/60))
  ((s=${1}%60))
  printf "%02d:%02d:%02d\n" $h $m $s
 }

# Set the protocol id. This could also be a command line argument.
PROTOCOL_ID="FastNoise"

# Set the paths.
# Set the paths by the computer.
YASS_SPIKE_PATH='/gscratch/retina/data/sorted'; #"/home/mike/ftp/files/ksort_out/"; #'/home/mike/ftp/files/20220420C/data021/kilosort2/';
RAW_DATA_PATH='/gscratch/scrubbed/retina/data/raw';
JAR_PATH="/gscratch/retina/GitRepos/manookin-lab/MEA/src/Vision7_for_2015DAQ/Vision.jar";
TMP_PATH="/gscratch/scrubbed/retina/data/sorted/${EXP}/"
VISIONPATH="/gscratch/retina/GitRepos/manookin-lab/MEA/src/Vision7_for_2015DAQ/"

SORT_PATH="${YASS_SPIKE_PATH}/${EXP}/"
DATA_PATH="${RAW_DATA_PATH}/${EXP}/"

data_string="${DATA_PATH}$FILE_NAME"

if [[ ! -e $TMP_PATH ]]; then
   mkdir ${TMP_PATH}
fi

# Get the current date and time.
NOW=`date +%Y-%m-%d_%H:%M:%S`

# Log the job.
echo "Job ${SLURM_JOB_ID} started for ${EXP}, ${FILE_NAME} on ${NOW}"

# Prepare the data.
singularity shell --nv --bind /gscratch /gscratch/retina/containers/a40.sif <<EOF
source /opt/conda/etc/profile.d/conda.sh
conda activate yass

python ../utilities/convert_litke_to_yass.py ${data_string} ${TMP_PATH} ${FILE_NAME} -w
EOF

# Load singularity.
module load singularity

singularity shell --nv --bind /gscratch /gscratch/retina/containers/a40.sif <<EOF
source /opt/conda/etc/profile.d/conda.sh
conda activate yass

# Run spike sorting in Vision.
java -d64 -Xmx8000m -Xss2m -cp ${JAR_PATH} edu.ucsc.neurobiology.vision.tasks.NeuronIdentification ${DATA_PATH}${FILE_NAME} ${SORT_PATH} -c ${VISIONPATH}config.xml
# java -Xmx8G -cp $JAR_PATH edu.ucsc.neurobiology.vision.calculations.CalculationManager "Electrophysiological Imaging Fast" ${SORT_PATH}${FILE_NAME}/ ${DATA_PATH}${FILE_NAME} 0.01 20 40 1000000 8

# Run the simple STA analysis.
python sta_single_file.py ${SORT_PATH}${FILE_NAME} ${FILE_NAME} ${SEED} -p ${PROTOCOL_ID}
EOF


echo "Job ${SLURM_JOB_ID} finished for ${EXP}, ${FILE_NAME} in $(format_time $SECONDS)" 