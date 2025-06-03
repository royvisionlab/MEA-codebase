#!/bin/bash -l

if [[ "$HOSTNAME" == *"hyak"* || "$HOSTNAME" == *"klone"* ]]; then
    TEMPORARY_SORT_PATH='/gscratch/scrubbed/retina/data/'
    LITKE_PATH='/gscratch/scrubbed/retina/data/'
    SORTED_SPIKE_PATH='/gscratch/retina/data/sorted'; #'/gscratch/scrubbed/retina/data/sorted';#"/home/mike/ftp/files/ksort_out/"; #'/home/mike/ftp/files/20220420C/data021/kilosort2/';
    RAW_DATA_PATH='/gscratch/scrubbed/retina/data/raw'; #'/gscratch/scrubbed/retina/data/raw';
    VISIONPATH="/gscratch/retina/GitRepos/manookin-lab/MEA/src/Vision7_for_2015DAQ/Vision.jar";
elif [[ "$HOSTNAME" == *"obsidian"* ]]; then
    TEMPORARY_SORT_PATH='/data/data/'
    LITKE_PATH='/data/data/' #'/usr/share/pool/data/'
    SORTED_SPIKE_PATH='/data/data/sorted'; #"/home/mike/ftp/files/ksort_out/"; #'/home/mike/ftp/files/20220420C/data021/kilosort2/';
    KILOSORT_TTL_PATH='/data/data/sorted'; #'/home/mike/ftp/files/ksort_out/';
    RAW_DATA_PATH='/data/data/raw';
    VISIONPATH='/home/mike/Documents/git_repos/manookin-lab/MEA/src/Vision7_for_2015DAQ/Vision.jar';
else
    TEMPORARY_SORT_PATH='/mnt/md0/data/'
    LITKE_PATH='/mnt/md0/data/' #'/usr/share/pool/data/'
    SORTED_SPIKE_PATH='/usr/share/pool/data/sorted'; #"/home/mike/ftp/files/ksort_out/"; #'/home/mike/ftp/files/20220420C/data021/kilosort2/';
    KILOSORT_TTL_PATH='/usr/share/pool/data/sorted'; #'/home/mike/ftp/files/ksort_out/';
    RAW_DATA_PATH='/usr/share/pool/data/raw';
    VISIONPATH="/home/mike/Documents/GitRepos/Manookin-Lab/MEA/src/Vision7_for_2015DAQ/Vision.jar";
fi