%% you need to change most of the paths in this block

% addpath(genpath('D:\GitHub\KiloSort')) % path to kilosort folder
addpath(genpath('/Volumes/Lab/Development/Kilosort3/'));
addpath('/Volumes/Lab/Development/Kilosort2/npy-matlab/npy-matlab');
rootZ = 'D:\20220412C_ksort'; % the raw data binary file is in this folder
rootH = 'D:\20220412C_ksort'; % path to temporary binary file (same size as data, should be on fast SSD)
pathToYourConfigFile = 'C:\Users\manoo\Documents\SpikeSorting\Python\StandardConfig.m'; % take from Github folder and put it somewhere else (together with the master_file)
chanMapFile = 'D:\20220412C_ksort\data038_map.mat';

ops.trange    = [0 Inf]; % time range to sort
ops.NchanTOT  = 512; % total number of channels in your recording

run(fullfile(pathToYourConfigFile));
ops.fproc   = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD
ops.chanMap = chanMapFile;
%% this block runs all the steps of the algorithm
fprintf('Looking for data inside %s \n', rootZ)

% main parameter changes from Kilosort2 to v2.5
ops.sig        = 20;  % spatial smoothness constant for registration
ops.fshigh     = 300; % high-pass more aggresively
ops.nblocks    = 5; % blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option. 

% main parameter changes from Kilosort2.5 to v3.0
ops.Th       = [9 9];

% is there a channel map file in this folder?
fs = dir(fullfile(rootZ, 'chan*.mat'));
if ~isempty(fs)
    ops.chanMap = fullfile(rootZ, fs(1).name);
end

% find the binary file
fs          = [dir(fullfile(rootZ, '*.bin')) dir(fullfile(rootZ, '*.dat'))];
ops.fbinary = fullfile(rootZ, fs(1).name);

rez                = preprocessDataSub(ops);
rez                = datashift2(rez, 1);

[rez, st3, tF]     = extract_spikes(rez);

rez                = template_learning(rez, tF, st3);

[rez, st3, tF]     = trackAndSort(rez);

rez                = final_clustering(rez, tF, st3);

rez                = find_merges(rez, 1);

rootZ = fullfile(rootZ, 'kilosort3');
mkdir(rootZ)
rezToPhy2(rez, rootZ);

exit;
