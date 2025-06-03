function run_kilosort(experiment_name, chunk_name, array_id, varargin)
% run_kilosort(experiment_name, chunk_name, array_id, varargin)
%
% Parameters:
%   experiment_name: The name of the experiment (e.g., '20231005C')
%   chunk_name: The name of the sorting chunk (e.g., 'chunk1')
%   array_id: The id of the array used in the experiment.
%   varargin:
%     version: The version of kilosort to run. Default is 2.5.
%     threshold: The threshold for spike detection. Default is -4.5.
%     useCAR: Use common average referencing? Default is false.
%     max_peels: Maximum number of times to peel a template. Default is 1000.
%
% Example:
%   run_kilosort('20231005C', 'chunk1', 504, 'version', 2.5, 'threshold', -4.5, 'useCAR', false, 'max_peels', 1000)
%

    ip = inputParser();
    ip.addParameter('version', 2.5, @(x)isfloat(x)); % Version of kilosort to run.
    ip.addParameter('threshold',-3.0,@(x)isfloat(x));
    ip.addParameter('useCAR', false, @(x)islogical(x)); % Use common average referencing?
    ip.addParameter('max_peels', 1000, @(x)isfloat(x)); % Maximum number of times to peel a template.
    ip.parse(varargin{:});
    
    version = ip.Results.version;
    threshold = ip.Results.threshold;
    useCAR = ip.Results.useCAR;
    max_peels = ip.Results.max_peels;

    if nargin < 3
        array_id = 504;
    end
    
    LITKE_519_120UM_ARRAY_ID_CARVEOUT = [1505,1506,1510,1601];
    if (array_id >= 500) && (array_id < 1500)
        map_name = 'LITKE_512_ARRAY.mat';
        if str2double(experiment_name(1:8)) <= 20230221
            map_name = 'chanMap512_jacked.mat';
        end
        total_channels = 512;
        variance_distance = 70;
    elseif (array_id >= 1500) && (array_id < 2500)
        if ismember(array_id, LITKE_519_120UM_ARRAY_ID_CARVEOUT)
            map_name = 'LITKE_519_ARRAY_120UM.mat';
            variance_distance = 140;
        else
            map_name = 'LITKE_519_ARRAY_30UM.mat';
            variance_distance = 35;
        end
        total_channels = 519;
    elseif (array_id >= 3500) && (array_id < 4000)
        map_name = 'LITKE_519_ARRAY_120UM.mat';
        total_channels = 519;
        variance_distance = 140;
    else
        error(['No such Litke array id ',num2str(array_id)]);
    end
    disp(['Using electrode map ',map_name,'.']);

%     [~, hostname] = system('hostname');
    hostname = getenv('HOSTNAME');
    if isempty(hostname)
        [~, hostname] = system('hostname');
    end
    this_path = pwd;
    disp(['Computer host name is: ', hostname]);

    % Get the paths. Check if you're on Hyak.
    if contains(hostname, 'orion', 'IgnoreCase', true)
        % Add path depending on the version.
        add_kilosort_path(version);
	addpath('/home/circuit/Documents/Kilosort/Kilosort-2.5-Kais/npy-matlab/npy-matlab/');
        addpath('/home/circuit/Documents/Kilosort/Kilosort-2.5-Kais/MEA/src/pipeline_utilities')
        folderPath = ['/mnt/gemini/roylab/Analyzed-MEAdata/Array/Analysis/Pipeline-Data/',experiment_name,'/'];
        dataDirectory = ['/mnt/gemini/rawdata/MEAdata/',experiment_name,'/'];
        outputDirectory = ['/mnt/gemini/roylab/Analyzed-MEAdata/Array/Analysis/Pipeline-Data/',experiment_name,'/'];
        tempDirectory = ['/mnt/gemini/roylab/Analyzed-MEAdata/Array/Analysis/Pipeline-Data/',experiment_name,'/'];  
        csvPath = ['/mnt/gemini/roylab/Analyzed-MEAdata/Array/Analysis/Pipeline-Data/',experiment_name,'/',chunk_name, '.csv'];
        channelMap = ['/home/circuit/Documents/Kilosort/Kilosort-2.5-Kais/MEA/src/pipeline_utilities/kilosort/', map_name];
    else
        [folderPath, dataDirectory, outputDirectory, tempDirectory, csvPath, channelMap] = kilosort_paths(experiment_name, chunk_name, map_name, version);
    end
    
    [chunk_names, file_numbers, durations] = parseMetaCSV(csvPath);

    obj = LitkeToKilosort(dataDirectory, outputDirectory, tempDirectory, ...
        'version', version);
    
    obj.setFileNumbers(file_numbers);
    obj.setFileNames(chunk_names);
    obj.setTimeSamples(durations);
    
    obj.setOptions();
    kilosortOptions = obj.kilosortOptions;
    kilosortOptions.NchanTOT = total_channels;
    kilosortOptions.sigmaMask = variance_distance;
    kilosortOptions.spkTh = threshold;
    kilosortOptions.CAR = useCAR;
    kilosortOptions.max_peels = max_peels;
    kilosortOptions.chanMap = channelMap;

    % Paths to raw binary and temporary binary files.
    kilosortOptions.fbinary = fullfile(folderPath,[chunk_name,'.bin']); %[folderPath, chunk_name, '.bin'];
    kilosortOptions.fproc = fullfile(folderPath,[chunk_name,'temp_wh.dat']); %[folderPath, chunk_name, 'temp_wh.dat'];

    % Display the threshold used.
    disp(['Using spike detection threshold: ', num2str(kilosortOptions.spkTh)]);

    % Process the data.
    results = preprocessDataSub( kilosortOptions );

    if version == 2.5
        results = datashift2(results,1);
    else
        results = clusterSingleBatches(results);
    end

    if version == 2.5
        iseed=1;
        results = learnAndSolve8b(results, iseed);
    else
        results = learnAndSolve8b(results);
    end

    results = find_merges(results, 1);

    % Final splits by SVD
    results = splitAllClusters(results, 1);

    % Final splits by amplitudes
    if version < 2.5
        results = splitAllClusters(results, 0);
    end
    
    % Decide on cutoff
    results = set_cutoff(results);
    
    if version == 2.5
        % eliminate widely spread waveforms (likely noise)
        results.good = get_good_units(results);
    end
    
    % Get the x/y coordinates of the spikes.
    xy = obj.getXYSpikeCoordinates(results);
    results.xy = xy;

    fprintf('Found %d good units \n', sum(results.good>0))
    
    % Save out the results.
    obj.saveSingleData(results, chunk_name);
    obj.saveBatchData(results);
end

function chunk_names = parseMetaTxt(txt)
    fileID = fopen(txt, 'r');
    tline = fgetl(fileID);
    fclose(fileID);
    chunk_names = strsplit(tline,' ');
end

function [chunk_names, file_numbers, durations] = parseMetaCSV(csvPath)
    chunk_names = {};
    file_numbers = [];
    durations = [];
    fileID = fopen(csvPath, 'r');
    while true
      tline = fgetl(fileID);
      if ~ischar(tline) 
          break; 
      end
      
      splits = strsplit(tline,',');
      durations = [durations, str2double(splits{2})]; %#ok<AGROW>
      splits = strsplit(splits{1},'/');
      fname = splits{end};
      chunk_names = [chunk_names,fname]; %#ok<AGROW>
      file_numbers = [file_numbers,str2double(strrep(fname,'data',''))]; %#ok<AGROW>
    end
    fclose(fileID);
end

function add_kilosort_path(version)
    hostname = getenv('HOSTNAME');
    if isempty(hostname)
        [~, hostname] = system('hostname');
    end
    if version == 3
            addpath(genpath('/home/mike/Documents/git_repos/Kilosort3/'));
    elseif version == 2.5
            addpath(genpath('/home/circuit/Documents/Kilosort/Kilosort-2.5.2'));
    else
            addpath(genpath('/home/mike/Documents/git_repos/Kilosort2/'));
    end
end

%%
function kilosortOptions = setOptions(outputDirectory, varargin)
            % Parse the input options
    ip = inputParser();
    ip.addParameter('trange', [0, Inf], @(x)isfloat(x)); % Time range to sort
    ip.addParameter('NchanTOT', 512, @(x)isfloat(x)); % Total number of channels in your recording
    ip.addParameter('fproc', [outputDirectory,'temp_wh.dat'], @(x)ischar(x)); % WTF is this?
    ip.addParameter('sig', 20, @(x)isfloat(x)); % Spatial smoothness constant for registration
    ip.addParameter('nblocks', 5, @(x)isfloat(x)); % Blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option.

    ip.addParameter('chanMap', ['/mmfs1/gscratch/retina/GitRepos/manookin-lab/MEA/src/pipeline_utilities/kilosort/','chanMap.mat'], @(x)ischar(x)); % Path to channel configuration .mat file.
    ip.addParameter('fs', 20000, @(x)isfloat(x)); % Sample rate during acquisition
    ip.addParameter('fshigh', 150, @(x)isfloat(x)); % High-pass filtering frequency in Hz
    ip.addParameter('minfr_goodchannels', 0.1, @(x)isfloat(x)); % Minimum firing rate on a "good" channel (0 to skip)
    ip.addParameter('Th', [10, 4], @(x)isfloat(x)); % Threshold on projections (like in Kilosort1, can be different for last pass like [10 4]). Default was [9,9]...
    ip.addParameter('lam', 10, @(x)isfloat(x)); % How important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot)
    ip.addParameter('AUCsplit', 0.97, @(x)isfloat(x));% was 0.97 % Splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
    ip.addParameter('minFR', 1/50, @(x)isfloat(x)); % Minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
    ip.addParameter('momentum', [20, 400], @(x)isfloat(x)); % Number of samples to average over (annealed from first to second value)
    ip.addParameter('sigmaMask', 70, @(x)isfloat(x));% was 70 % Spatial constant in um for computing residual variance of spike
    ip.addParameter('ThPre', 8, @(x)isfloat(x)); % Threshold crossings for pre-clustering (in PCA projection space)

    % TODO: check if you can wavelet filter instead of normal
    % highpass filtering.

    % Don't mess with the following parameters.
    ip.addParameter('spkTh', -4.5, @(x)isfloat(x)); % was 4.5 % Spike threshold in standard deviations (-6)
    ip.addParameter('reorder', true, @(x)islogical(x)); % Whether to reorder batches for drift correction.
    ip.addParameter('nskip', 25, @(x)isfloat(x)); % How many batches to skip for determining spike PCs
    ip.addParameter('GPU', true, @(x)islogical(x)); % Is this a boolean??? % Has to be 1, no CPU version yet, sorry
    ip.addParameter('nfilt_factor', 4, @(x)isfloat(x)); % Max number of clusters per good channel (even temporary ones)
    ip.addParameter('ntbuff', 64, @(x)isfloat(x)); % Samples of symmetrical buffer for whitening and spike detection
    ip.addParameter('NT', 64*1024 - 64, @(x)isfloat(x)); % Batch size; Must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory).
    ip.addParameter('whiteningRange', 32, @(x)isfloat(x)); % Number of channels to use for whitening each channel
    ip.addParameter('nSkipCov', 25, @(x)isfloat(x)); % Compute whitening matrix from every N-th batch
    ip.addParameter('scaleproc', 200, @(x)isfloat(x)); % int16 scaling of whitened data
    ip.addParameter('nPCs', 3, @(x)isfloat(x)); % How many PCs to project the spikes into
    ip.addParameter('useRAM', false, @(x)islogical(x)); % Make sure this is boolean; pretty sure...

    % Parse the inputs.
    ip.parse(varargin{:});

    % Get the field names from the input parser.
    fnames = fieldnames(ip.Results);

    % Copy the parameters into the options structure.
    kilosortOptions = struct();
    for jj = 1 : length(fnames)
        kilosortOptions.(fnames{jj}) = ip.Results.(fnames{jj});
    end
end
