addpath('/Volumes/Lab/Development/Kilosort2/npy-matlab/npy-matlab');
d = gpuDevice(1);


addpath(genpath('/Volumes/Lab/Development/Kilosort2')); % path to kilosort folder
%addpath('D:\GitHub\npy-matlab') % for converting to Phy
rootZ = '{0}'; % the raw data binary file is in this folder
rootH = '{0}'; % path to temporary binary file (same size as data, should be on fast SSD)
pathToYourConfigFile = '{1}'; % take from Github folder and put it somewhere else (together with the master_file)
chanMapFile = '{2}';


ops.trange = [0 Inf]; % time range to sort
ops.NchanTOT    = {3}; % total number of channels in your recording


run(fullfile(pathToYourConfigFile));
ops.fproc       = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD
%ops.chanMap = fullfile(pathToYourConfigFile, chanMapFile);
ops.chanMap = chanMapFile;

%% this block runs all the steps of the algorithm
fprintf('Looking for data inside %s \n', rootZ)

% is there a channel map file in this folder?
%fs = dir(fullfile(rootZ, 'chan*.mat'));
%if ~isempty(fs)
%    ops.chanMap = fullfile(rootZ, fs(1).name);
%end

% find the binary file
fs          = [dir(fullfile(rootZ, '*.dat')) dir(fullfile(rootZ, '*.dat'))];
ops.fbinary = fullfile(rootZ, fs(1).name);

% preprocess data to create temp_wh.dat
rez = preprocessDataSub(ops);

% time-reordering as a function of drift
rez = clusterSingleBatches(rez);

% saving here is a good idea, because the rest can be resumed after loading rez
save(fullfile(rootZ, 'rez.mat'), 'rez', '-v7.3');

% main tracking and template matching algorithm
%load(fullfile(rootZ, 'rez.mat'));
rez = learnAndSolve8b(rez);

% final merges
rez = find_merges(rez, 1);

% final splits by SVD
rez = splitAllClusters(rez, 1);

% final splits by amplitudes
rez = splitAllClusters(rez, 0);

% decide on cutoff
rez = set_cutoff(rez);

fprintf('found %d good units \n', sum(rez.good>0))

% write to Phy
fprintf('Saving results to Phy  \n')
rezToPhy(rez, rootZ);

%% if you want to save the results to a Matlab file...

% discard features in final rez file (too slow to save)
rez.cProj = [];
rez.cProjPC = [];

% save final results as rez2
fprintf('Saving final results in rez2  \n')
fname = fullfile(rootZ, 'rez2.mat');
save(fname, 'rez', '-v7.3');

exit;
