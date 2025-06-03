#!/bin/bash -l

# Usage: sbatch --export=EXP='20240105C',FILES='data000 data001',SEEDS='1113434 11134355',NUM_EPOCHS='20 5' online_blue.slurm

#SBATCH --job-name=online_blue        # Job name
#SBATCH --mail-type=ALL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=vyomr@uw.edu     # Where to send mail	

#SBATCH --account=retina
#SBATCH --partition=gpu-a40 #ckpt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=98G
#SBATCH --gpus=1
# #SBATCH --constraint="a100|a40|rtx6k"
#SBATCH --time=0-05:00:00 # Max runtime in DD-HH:MM:SS format.

#SBATCH --chdir=/gscratch/retina/GitRepos/manookin-lab/MEA/src/online_analysis/
#SBATCH --export=all
#SBATCH --output=/gscratch/retina/logs/log_%j.log # where STDOUT goes
#SBATCH --error=/gscratch/retina/logs/error_%j.log # where STDERR goes

umask 007

format_time() {
  ((h=${1}/3600))
  ((m=(${1}%3600)/60))
  ((s=${1}%60))
  printf "%02d:%02d:%02d\n" $h $m $s
 }

DATA_FILES=($FILES)
SEEDS=($SEEDS)
NUM_EPOCHS=($NUM_EPOCHS)

FILE_NAME=${DATA_FILES[0]}

# Set the protocol id. This could also be a command line argument.
PROTOCOL_ID="BlueNoise"

echo "Using seed ${SEEDS[0]}"

# Set the paths.
# Set the paths by the computer.
SORTED_SPIKE_PATH='/gscratch/scrubbed/retina/data/sorted'; #"/home/mike/ftp/files/ksort_out/"; #'/home/mike/ftp/files/20220420C/data021/kilosort2/';
RAW_DATA_PATH='/gscratch/scrubbed/retina/data/raw';
VISION_PATH="/gscratch/retina/GitRepos/manookin-lab/MEA/src/Vision7_for_2015DAQ/Vision.jar";

SORT_PATH="${SORTED_SPIKE_PATH}/${EXP}/"
DATA_PATH="${RAW_DATA_PATH}/${EXP}/"

# data_string="${DATA_PATH}$FILE_NAME"
data_string=""
for arg in "${DATA_FILES[@]}"; do
   data_string=$data_string"${DATA_PATH}$arg "
done

if [[ ! -e $SORT_PATH ]]; then
   mkdir ${SORT_PATH}
fi

if [[ ! -e $SORT_PATH ]]; then
   mkdir ${SORT_PATH}
fi

# Get the current date and time.
NOW=`date +%Y-%m-%d_%H:%M:%S`

# Log the job.
echo "Job ${SLURM_JOB_ID} started for ${EXP}, ${FILE_NAME} on ${NOW}"

# Load singularity.
module load singularity

# Prepare the data. Converts raw data files with headers and 1st TTL pulse channel (Litke format), 
# and writes remaining 512 channels to binary file
singularity shell --nv --bind /gscratch /gscratch/retina/containers/a40.sif <<EOF
umask 007
source /opt/conda/etc/profile.d/conda.sh
conda activate yass

python convert_join_litke_datasets_for_yass.py ${data_string} ${SORT_PATH} ${FILE_NAME} -b
EOF

# Modules to use (optional).
module load matlab

# Run kilosort, which outputs .npy files, then convert to vision format, and run STA analysis using .neurons output.
matlab -nodisplay -r "run_kilosort_single('${EXP}','${FILE_NAME}'); quit"

singularity shell --nv --bind /gscratch /gscratch/retina/containers/a40.sif <<EOF
umask 007
source /opt/conda/etc/profile.d/conda.sh
conda activate yass

# Convert to Vision format.
python ../utilities/kilosort_to_vision.py ${SORT_PATH}${FILE_NAME}/kilosort2 ${SORT_PATH}${FILE_NAME} ${SORT_PATH}${FILE_NAME} ${FILE_NAME} -l -d ${RAW_DATA_PATH}/${EXP}/${FILE_NAME}/

# Get the EI file.
java -Xmx8G -cp $VISION_PATH edu.ucsc.neurobiology.vision.calculations.CalculationManager "Electrophysiological Imaging Fast" ${SORT_PATH}${FILE_NAME}/ ${RAW_DATA_PATH}/${EXP}/${FILE_NAME}/ 0.01 20 40 1000000 8

# Now run deduplication on the EI files.
echo "Running deduplication on the EI files."
python ../analysis/deduplication/deduplication.py ${EXP} -a kilosort2 -c ${FILE_NAME} -o

# Run the simple STA analysis.
python sta_single_file.py ${SORT_PATH}${FILE_NAME} ${FILE_NAME} ${SEED} -p ${PROTOCOL_ID} -e ${EXP}

EOF

echo "Job ${SLURM_JOB_ID} finished for ${EXP}, ${FILE_NAME} in $(format_time $SECONDS)" 