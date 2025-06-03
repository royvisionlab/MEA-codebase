#!/bin/bash

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: $0 <PASSWORD> <EXPERIMENT_DATE> <CHUNK_NAMES>"
    echo "Example: bash $0 password 20240926C chunk1 chunk2 chunk3"
    exit 0
fi

# USAGE: bash transfer_chunks.sh <password> <experiment> chunk1 chunk2 ...

PASS=$1
shift
EXP=$1
shift 
CHUNKS="$*"

# Make sure that the first 8 characters of the password are not a number.
if [[ "${PASS:0:8}" =~ ^[0-9]{8}$ ]]; then
  echo "Invalid password: ${PASS}... exiting"
  exit 1
fi

# Make sure that the first argument is a valid experiment folder.
if [[ "${EXP:0:8}" =~ ^[0-9]{8}$ ]]; then
  echo "Experiment: ${EXP} is a valid experiment folder."
else
  echo echo "Experiment: ${EXP} is NOT a valid experiment folder... exiting"
  exit 1
fi

# ALG='kilosort2.5'

chunk_names=($CHUNKS)

SERVER="mea@128.95.10.105"

echo "Transferring experiment folder."
sshpass -p ${PASS} \scp -r /data/data/sorted/${EXP} ${SERVER}:/volume1/data/data/sorted/

echo "Transferring chunk files."
for arg in "${chunk_names[@]}"; do
    echo "Transferring files from ${arg} to NAS."
    ALGS=( /data/data/sorted/${EXP}/${arg}/*/ )
    for ALG in "${ALGS[@]}"; do
        ALG=$(basename $ALG)
        echo "Transferring ${ALG} files."
        sshpass -p ${PASS} \ssh ${SERVER} "mkdir -p /volume1/data/analysis/${EXP}/${arg}/${ALG}/"
        sshpass -p ${PASS} \scp -v /data/data/sorted/${EXP}/${arg}/${ALG}/${ALG}* ${SERVER}:/volume1/data/analysis/${EXP}/${arg}/${ALG}/
    done
    # ssh {$SERVER} "mkdir -p /volume1/data/analysis/${EXP}/${arg}/${ALG}/"
    # echo "\scp /data/data/sorted/${EXP}/${arg}/${ALG}/${ALG}*"
    # sshpass -p ${PASS} \scp /data/data/sorted/${EXP}/${arg}/${ALG}/${ALG}* mea@128.95.10.105:/volume1/data/data/sorted/${EXP}/${arg}/${ALG}/
    # sshpass -p ${PASS} \scp /data/data/sorted/${EXP}/${arg}/${ALG}/${ALG}* ${SERVER}:/volume1/data/analysis/${EXP}/${arg}/${ALG}/
done

# ssh ${SERVER} "for dir in /volume1/data/analysis/20240820C/chunk2/*/; do echo ${dir}; done"
# array=( /Volumes/*/ )

# ssh ${SERVER} "alg_dirs=(/volume1/data/analysis/20240820C/chunk2/*/); for dir in ${alg_dirs[@]}; do echo ${dir}; done"

# ls /data/data/sorted/20240820C/chunk2

# sshpass -p ${PASS} \scp /data/data/sorted/${EXP}/chunk1/${ALG}/${ALG}* mea@128.95.10.183:/volume1/data/data/sorted/${EXP}/chunk1/${ALG}/
