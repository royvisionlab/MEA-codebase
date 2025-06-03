#!/bin/bash -l

# Usage: sbatch --export=EXP='20230105C',CHUNK='chunk1',ALG='kilosort2' process_spikes.slurm

#SBATCH --job-name=process_spikes       # Job name
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=manookin@uw.edu     # Where to send mail	

#SBATCH --account=retina
#SBATCH --partition=ckpt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=256G
#SBATCH --gpus=0
# #SBATCH --constraint="mem"
#SBATCH --time=0-12:00:00 # Max runtime in DD-HH:MM:SS format.

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

# Log the job.
echo "Job ${SLURM_JOB_ID} started for ${EXP}, ${CHUNK} on ${NOW}" >> /gscratch/retina/logs/process_spikes.txt

# Set the paths.
KILOSORT_SPIKE_PATH='/gscratch/retina/data/sorted'; #"/home/mike/ftp/files/ksort_out/"; #'/home/mike/ftp/files/20220420C/data021/kilosort2/'; 
KILOSORT_TTL_PATH='/gscratch/retina/data/sorted';
RAW_DATA_PATH='/gscratch/scrubbed/retina/data/raw';
VISIONPATH="/gscratch/retina/GitRepos/manookin-lab/MEA/src/Vision7_for_2015DAQ/Vision.jar";
TMP_PATH="/gscratch/retina/data/sorted/${EXP}/"

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
SPIKE_PATH="${KILOSORT_SPIKE_PATH}/${CHUNK}/${ALG}/"
VISION_OUT="${KILOSORT_SPIKE_PATH}/${CHUNK}/${ALG}/"
SORT_PATH=${KILOSORT_SPIKE_PATH}

for arg in "${data_files[@]}"; do
    FNAME=$arg
    echo "Analyzing file ${FNAME}"
    master_dir="${SORT_PATH}/${FNAME}"
    if [[ ! -e $master_dir ]]; then
        mkdir -p $master_dir
    fi
    echo "Master directory is ready."

    if [[ ! -e ${master_dir}/${ALG} ]]; then
        mkdir -p ${master_dir}/${ALG}
    fi
done

# Load singularity.
module load singularity

# Combine the raw Litke bin files into a binary file. (Also works for a single file)
# singularity shell --bind /gscratch /gscratch/retina/containers/a40.sif <<EOF
singularity shell --bind /gscratch /gscratch/retina/containers/a40.sif <<EOF
source /opt/conda/etc/profile.d/conda.sh
conda activate yass


if [ "${ALG}" = "kilosort2" ]; then
    echo "Collecting spikes from kilosort2"
    python kilosort_to_vision.py $SPIKE_PATH $KILOSORT_TTL_PATH $VISION_OUT kilosort2 -l
else
    echo "Collecting spikes from yass"
    python yass_output_to_neurons.py ${SPIKE_PATH}tmp/output/spike_train.npy ${SPIKE_PATH}tmp/ ${SORT_PATH} yass
    # Split the spikes into the appropriate data files.
    python split_yass_npy.py ${SPIKE_PATH}tmp/output/spike_train.npy ${SORT_PATH}/${CHUNK}.csv ${SORT_PATH}
fi

for arg in "${data_files[@]}"; do
    FNAME=$arg
    echo "Analyzing file ${FNAME}"
    master_dir="${SORT_PATH}/${FNAME}"
    # if [[ ! -e $master_dir ]]; then
    #     mkdir -p $master_dir
    # fi
    # echo "Master directory is ready."

    # if [[ ! -e ${master_dir}/${ALG} ]]; then
    #     mkdir -p ${master_dir}/${ALG}
    # fi
    echo "Algorithm directory is ready."
    SPIKE_PATH="${KILOSORT_SPIKE_PATH}/${FNAME}/${ALG}/"
    DATA_PATH="${RAW_DATA_PATH}/${FNAME}/"
    VISION_OUT="${KILOSORT_SPIKE_PATH}/${FNAME}/"

    if [ "${ALG}" = "kilosort2" ]; then
        python kilosort_to_vision.py $SPIKE_PATH $KILOSORT_TTL_PATH $VISION_OUT $FNAME -l -d $DATA_PATH
    else
        python yass_output_to_neurons.py ${SORT_PATH}/${FNAME}.npy ${SPIKE_PATH}tmp/ ${master_dir} ${FNAME} -d ${DATA_PATH}${FNAME}
    fi

    # Fix the EI map.
    python fix_electrode_map.py ${KILOSORT_SPIKE_PATH} -a ${ALG} -f ${FNAME}
done
EOF

module load openjdk
for arg in "${data_files[@]}"; do
    FNAME=$arg
    echo "Analyzing file ${FNAME}"
    master_dir="${SORT_PATH}/${FNAME}"
    SPIKE_PATH="${KILOSORT_SPIKE_PATH}/${FNAME}/${ALG}/"
    DATA_PATH="${RAW_DATA_PATH}/${FNAME}/"
    VISION_OUT="${KILOSORT_SPIKE_PATH}/${FNAME}/"

    java -Xmx8G -cp $VISIONPATH edu.ucsc.neurobiology.vision.calculations.CalculationManager "Electrophysiological Imaging Fast" ${KILOSORT_SPIKE_PATH}/${FNAME}/ $DATA_PATH 0.01 67 133 1000000 16
    java -Xmx8G -cp $VISIONPATH edu.ucsc.neurobiology.vision.calculations.CalculationManager "Raw Data Noise Evaluation" $DATA_PATH ${KILOSORT_SPIKE_PATH}/${FNAME}/
    # Now move the Vision files into the sort directory.
    mv "${KILOSORT_SPIKE_PATH}/${FNAME}/${FNAME}.ei" ${SPIKE_PATH}
    mv "${KILOSORT_SPIKE_PATH}/${FNAME}/${FNAME}.globals" ${SPIKE_PATH}
    mv "${KILOSORT_SPIKE_PATH}/${FNAME}/${FNAME}.neurons" ${SPIKE_PATH}
    mv "${KILOSORT_SPIKE_PATH}/${FNAME}/${FNAME}.noise" ${SPIKE_PATH}
    if [ "${ALG}" = "yass" ]; then
        mv "${SORT_PATH}/${FNAME}.npy" ${SPIKE_PATH}
    fi
done

singularity shell --bind /gscratch /gscratch/retina/containers/a40.sif <<EOF
source /opt/conda/etc/profile.d/conda.sh
conda activate yass

for arg in "${data_files[@]}"; do
    FNAME=$arg
    echo "Analyzing file ${FNAME}"
    master_dir="${SORT_PATH}/${FNAME}"
    SPIKE_PATH="${KILOSORT_SPIKE_PATH}/${FNAME}/${ALG}/"
    DATA_PATH="${RAW_DATA_PATH}/${FNAME}/"
    VISION_OUT="${KILOSORT_SPIKE_PATH}/${FNAME}/"

    # Fix the EI map.
    python fix_electrode_map.py ${KILOSORT_SPIKE_PATH} -a ${ALG} -f ${FNAME}
done
EOF

echo "Job ${SLURM_JOB_ID} finished for ${EXP}, ${CHUNK} in $(format_time $SECONDS)" >> /gscratch/retina/logs/process_spikes.txt
