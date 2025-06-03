#!/bin/bash

# USAGE: bash litke_to_yass.sh 20220531C noise data005 data006 data010

data_files=""
litke_bin_path=""
dest_analysis_path=""
data_string=""

# litke_bin_path=$1
# shift
experimentName=$1
shift
out_name=$1
shift

data_files=( "$@" )


# while [ "$1" != "" ]; do
#   case $1 in
#   -d)
#     data_files="$1"
#     shift
#     ;;
#   *)
#     litke_bin_path=$1
#     shift
#     dest_analysis_path=$1
#     shift
#     ;;
#   esac
# done

# echo $litke_bin_path

# for arg in "${data_files[@]}"; do
#    echo "$arg"
# done


# experimentName='20220531C'
YASS_SPIKE_PATH='/usr/share/pool/SortedData'; #"/home/mike/ftp/files/ksort_out/"; #'/home/mike/ftp/files/20220420C/data021/kilosort2/';
KILOSORT_TTL_PATH='/usr/share/pool/SortedData'; #'/home/mike/ftp/files/ksort_out/';
RAW_DATA_PATH='/usr/share/pool/rawdata/MEA_Data';
VISIONPATH="/home/mike/Documents/GitRepos/Manookin-Lab/MEA/src/Vision7_for_2015DAQ/Vision.jar";

TMP_PATH='/media/mike/NVME/sort/'
SORT_PATH="${YASS_SPIKE_PATH}/${experimentName}/"
DATA_PATH="${RAW_DATA_PATH}/${experimentName}/"

data_string=""
for arg in "${data_files[@]}"; do
   data_string=$data_string"${DATA_PATH}$arg "
done
# echo $data_string

# Step 1: Combine the raw Litke bin files into a yass file. (Also works for a single file)

python convert_join_litke_datasets_for_yass.py ${data_string} ${TMP_PATH} ${out_name} -m -b
