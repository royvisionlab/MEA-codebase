% Place this file in your Matlab root directory or add it to your path through your startup.m file.

function [folderPath, dataDirectory, outputDirectory, tempDirectory, csvPath, channelMap] = kilosort_paths(experiment_name, chunk_name, map_name, kilosort_version)
% [folderPath, dataDirectory, outputDirectory, tempDirectory, csvPath, channelMap] = kilosort_paths(experiment_name, chunk_name, map_name)
%
% This function sets the paths for the kilosort pipeline.
% 
% Inputs:
%   experiment_name: The name of the experiment (e.g., '20231005C')
%   chunk_name: The name of the sorting chunk (e.g., 'chunk1')
%   map_name: The name of the channel map (e.g., 'LITKE_519_ARRAY_120UM.mat')
%
% Outputs:
%   folderPath: The path to the experiment folder.
%   dataDirectory: The path to the raw data files for the experiment.
%   outputDirectory: The path to the output directory for the spike sorting.
%   tempDirectory: The path to the temporary directory for the spike sorting.
%   csvPath: The path to the output csv file that keeps track of which files belong to the chunk.
%   channelMap: The path to the channel map.

% Path to the experiment sorting folder.
folderPath = fullfile('/data/data/sorted/',experiment_name);

% Path to the raw data files for the experiment.
dataDirectory = fullfile('/data/data/raw/',experiment_name);

% Path to the output directory for the spike sorting.
outputDirectory = fullfile('/data/data/sorted/',experiment_name);

% Path to the temporary directory for the spike sorting.
tempDirectory = fullfile('/data/data/sorted/',experiment_name);  

% Path to the output csv file that keeps track of which files belong to the chunk.
csvPath = fullfile('/data/data/sorted/',experiment_name,[chunk_name, '.csv']);

% Path to the channel map file.
channelMap = fullfile('/home/mike/Documents/git_repos/manookin-lab/MEA/src/pipeline_utilities/kilosort/',map_name);

% Add any necessary paths here if you can't do it in your startup.m file.
addpath('/home/mike/Documents/git_repos/npy-matlab/npy-matlab/');
addpath('/home/mike/Documents/git_repos/manookin-lab/MEA/src/pipeline_utilities/');

% Add path to kilosort based on version.
if kilosort_version == 2.5
    addpath(genpath('/home/mike/Documents/git_repos/manookin-lab/kilosort25/'));
elseif kilosort_version == 2
    addpath(genpath('/home/mike/Documents/git_repos/Kilosort2/'));
elseif kilosort_version == 3
    addpath(genpath('/home/mike/Documents/git_repos/Kilosort3/'));
else
    error('Matlab-based Kilosort version not recognized.');
end