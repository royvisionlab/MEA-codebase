#!/bin/bash -l

# Usage: sbatch --export=EXP='20230105C',CHUNK='chunk1',ALG='kilosort2.5' process_spikes.slurm
# with defaults: sbatch --export=EXP='20230105C',CHUNK='chunk1' process_spikes.slurm

#SBATCH --job-name=process_spikes       # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=manookin@uw.edu     # Where to send mail	

#SBATCH --account=retina
#SBATCH --partition=ckpt-all #gpu-a40 #ckpt-all
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=119G #256G
#SBATCH --gpus=0
# #SBATCH --constraint="mem"
#SBATCH --time=0-12:00:00 # Max runtime in DD-HH:MM:SS format.

#SBATCH --chdir=/gscratch/retina/GitRepos/manookin-lab/MEA/src/utilities/
#SBATCH --export=all
#SBATCH --output=/gscratch/retina/logs/log_%j.log # where STDOUT goes
#SBATCH --error=/gscratch/retina/logs/error_%j.log # where STDERR goes

DEFAULT_ALGORITHM="kilosort2.5"
ALG="${ALG:-$DEFAULT_ALGORITHM}"

format_time() {
  ((h=${1}/3600))
  ((m=(${1}%3600)/60))
  ((s=${1}%60))
  printf "%02d:%02d:%02d\n" $h $m $s
 }

# Get the current date and time.
 NOW=`date +%Y-%m-%d_%H:%M:%S`

 # Set the number of CPUs to use. Hard code for now.
 NUM_CPU='6'

# Log the job.
echo "Job ${SLURM_JOB_ID} started for ${EXP}, ${CHUNK} on ${NOW}" >> /gscratch/retina/logs/process_spikes.txt

SORTED_SPIKE_PATH='/gscratch/retina/data/sorted'; 
KILOSORT_TTL_PATH='/gscratch/retina/data/sorted';
RAW_DATA_PATH='/gscratch/scrubbed/retina/data/raw';
JAR_PATH="/gscratch/retina/GitRepos/manookin-lab/MEA/src/Vision7_for_2015DAQ/Vision.jar"; 

TMP_PATH="${SORTED_SPIKE_PATH}/${EXP}/"

data_files=($(head -n 1 "${TMP_PATH}${EXP}_${CHUNK}.txt"))

data_string=""
for arg in "${data_files[@]}"; do
   data_string=$data_string"${DATA_PATH}$arg "
done

# Concatenate the experiment name.
EXPERIMENT_SPIKE_PATH="${SORTED_SPIKE_PATH}/${EXP}";
KILOSORT_TTL_PATH="${KILOSORT_TTL_PATH}/${EXP}/TTLTriggers";
EXPERIMENT_DATA_PATH="${RAW_DATA_PATH}/${EXP}";

CHUNK_SPIKE_PATH="${EXPERIMENT_SPIKE_PATH}/${CHUNK}/${ALG}/"
VISION_OUT="${EXPERIMENT_SPIKE_PATH}/${CHUNK}/${ALG}/"
SORT_PATH=${EXPERIMENT_SPIKE_PATH} 

# Check if the output directory exists. If not, create it.
echo "Checking if the output directory exists."
for arg in "${data_files[@]}"; do
    FNAME=$arg

    FILE_SPIKE_PATH="${EXPERIMENT_SPIKE_PATH}/${FNAME}/${ALG}/"
    FILE_DATA_PATH="${EXPERIMENT_DATA_PATH}/${FNAME}/"
    VISION_OUT="${EXPERIMENT_SPIKE_PATH}/${FNAME}/"

    # Create the directories if they don't exist.
    if [[ ! -e $VISION_OUT ]]; then
        mkdir -p $VISION_OUT
    fi

    if [[ ! -e $FILE_SPIKE_PATH ]]; then
        mkdir -p $FILE_SPIKE_PATH
    fi
done

# Check whether this is an old file and the EI needs to be corrected.
EXP_NUM=$(echo $EXP | cut -c1-8)

if (( EXP_NUM < 20230228 )); then
    OLD_EXP='True'
else
    OLD_EXP='False'
fi

# Load singularity.
module load singularity


singularity shell --bind /gscratch /gscratch/retina/containers/a40.sif <<EOF
source /opt/conda/etc/profile.d/conda.sh
conda activate yass

# Generate the neurons file for the chunk.
if [ "${ALG}" = "yass" ]; then
    echo "Collecting spikes from yass"
    python yass_output_to_neurons.py ${CHUNK_SPIKE_PATH}tmp/output/spike_train.npy ${CHUNK_SPIKE_PATH}tmp/ ${CHUNK_SPIKE_PATH} yass
    # Split the spikes into the appropriate data files.
    python split_yass_npy.py ${CHUNK_SPIKE_PATH}tmp/output/spike_train.npy ${EXPERIMENT_SPIKE_PATH}/${CHUNK}.csv ${EXPERIMENT_SPIKE_PATH}
else
    echo "Collecting spikes from ${ALG}"
    python kilosort_to_vision.py $CHUNK_SPIKE_PATH $KILOSORT_TTL_PATH ${EXPERIMENT_SPIKE_PATH}/${CHUNK}/${ALG}/ ${ALG} -l -q
fi
EOF

for arg in "${data_files[@]}"; do
    FNAME=$arg
    echo "Converting file ${FNAME} to Vision format."

    FILE_SPIKE_PATH="${EXPERIMENT_SPIKE_PATH}/${FNAME}/${ALG}/"
    FILE_DATA_PATH="${EXPERIMENT_DATA_PATH}/${FNAME}/"
    VISION_OUT="${EXPERIMENT_SPIKE_PATH}/${FNAME}/"

    singularity shell --bind /gscratch /gscratch/retina/containers/a40.sif <<EOF
    source /opt/conda/etc/profile.d/conda.sh
    conda activate yass

    if [ "${ALG}" = "yass" ]; then
        # python yass_output_to_neurons.py ${SORT_PATH}/${FNAME}.npy ${SPIKE_PATH}tmp/ ${master_dir} ${FNAME} -d ${DATA_PATH}${FNAME}
        python yass_output_to_neurons.py ${EXPERIMENT_SPIKE_PATH}/${FNAME}.npy ${EXPERIMENT_SPIKE_PATH}tmp/ ${VISION_OUT} ${FNAME} -d ${FILE_DATA_PATH}
    else
        python kilosort_to_vision.py $FILE_SPIKE_PATH $KILOSORT_TTL_PATH $VISION_OUT $FNAME -l -q -d $FILE_DATA_PATH
        # python kilosort_to_vision.py $SPIKE_PATH $KILOSORT_TTL_PATH $VISION_OUT $FNAME -l -d $DATA_PATH
    fi

    java -Xmx8G -cp $JAR_PATH edu.ucsc.neurobiology.vision.calculations.CalculationManager "Electrophysiological Imaging Fast" $VISION_OUT $FILE_DATA_PATH 0.01 67 133 1000000 ${NUM_CPU}
    java -Xmx8G -cp $JAR_PATH edu.ucsc.neurobiology.vision.calculations.CalculationManager "Raw Data Noise Evaluation" $FILE_DATA_PATH $VISION_OUT
    # Now move the Vision files into the sort directory.
    mv "${VISION_OUT}${FNAME}.ei" ${FILE_SPIKE_PATH}
    mv "${VISION_OUT}${FNAME}.globals" ${FILE_SPIKE_PATH}
    mv "${VISION_OUT}${FNAME}.neurons" ${FILE_SPIKE_PATH}
    mv "${VISION_OUT}${FNAME}.noise" ${FILE_SPIKE_PATH}
    if [ "${ALG}" = "yass" ]; then
        mv "${EXPERIMENT_SPIKE_PATH}/${FNAME}.npy" ${FILE_SPIKE_PATH}
    fi

    # Fix the EI map.
    if [ "${OLD_EXP}" = "True" ]; then
        python fix_electrode_map.py ${EXPERIMENT_SPIKE_PATH} -a ${ALG} -f ${FNAME}
    fi
EOF
done

echo "Job ${SLURM_JOB_ID} finished for ${EXP}, ${CHUNK} in $(format_time $SECONDS)" >> /gscratch/retina/logs/process_spikes.txt
