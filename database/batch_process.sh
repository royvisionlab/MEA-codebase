#!/bin/bash

H5_PATH="/data/data/h5/"
# H5_PATH='/Users/michaelmanookin/Documents/Data/datajoint/mea/data/'
JSON_PATH="/data/data/metadata/json/"

for file in ${H5_PATH}*.h5
do
    # echo "Processing file: $file"
    IFS='/' read -ra ADDR <<< "$file"
    # echo ${ADDR[@]}
    FNAME=${ADDR[${#ADDR[@]}-1]}
    if (( ${FNAME:0:8} > 20240522)); then
        echo "Processing experiment ${FNAME:0:9}."
        bash batch_parse.sh ${FNAME:0:9} > ${JSON_PATH}logs/${FNAME:0:9}.log
    fi
    # echo ${ADDR[@]:0:8}
    # EXP=$(basename $file .h5)
    # echo "Parsing Symphony H5 data for ${EXP}."
    # python parse_data.py ${H5_PATH}${EXP}.h5 ${H5_PATH}metadata/json/ -r ${H5_PATH}raw/
    # echo "Moving summary TXT file to text folder."
    # mv ${H5_PATH}metadata/json/${EXP}.txt ${H5_PATH}metadata/json/txt/
    # echo "Moving JSON file to datajoint folder."
    # mv ${H5_PATH}metadata/json/${EXP}_dj.json ${H5_PATH}metadata/json/datajoint/
done
