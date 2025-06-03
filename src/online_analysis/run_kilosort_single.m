function run_kilosort_single(experiment_name, file_name, array_id)

    if nargin < 3
        array_id = 504;
    end

    addpath('../pipeline_utilities');
    
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
    
    hostname = getenv('HOSTNAME');
    if isempty(hostname)
        [~, hostname] = system('hostname');
    end
    if contains(hostname,'mike','IgnoreCase',true) || contains(hostname,'maverick','IgnoreCase',true)
        folder_path = ['/usr/share/pool/data/sorted/',experiment_name,'/'];
        dataDirectory = ['/usr/share/pool/data/raw/',experiment_name,'/'];
        outputDirectory = ['/usr/share/pool/data/sorted/',experiment_name,'/'];
        tempDirectory = ['/usr/share/pool/data/sorted/',experiment_name,'/'];  
        txt = ['/usr/share/pool/data/sorted/',experiment_name,'/',file_name, '.csv'];
    elseif contains(hostname,'obsidian','IgnoreCase',true) 
        addpath(genpath('/home/mike/Documents/git_repos/Kilosort2'));
        folder_path = ['/data/data/sorted/',experiment_name,'/'];
        dataDirectory = ['/data/data/raw/',experiment_name,'/'];
        outputDirectory = ['/data/data/sorted/',experiment_name,'/'];
        tempDirectory = ['/data/data/sorted/',experiment_name,'/'];  
        txt = ['/data/data/sorted/',experiment_name,'/',file_name, '.csv'];
    else
        % Add path to the kilosort directory.
        addpath(genpath('/mmfs1/gscratch/retina/GitRepos/chichilnisky-lab/Kilosort2'));
        addpath('/mmfs1/gscratch/retina/GitRepos/npy-matlab/npy-matlab/');
        addpath('/mmfs1/gscratch/retina/GitRepos/manookin-lab/MEA/src/pipeline_utilities/')
        folder_path = ['/mmfs1/gscratch/scrubbed/retina/data/sorted/',experiment_name,'/'];
        dataDirectory = ['/mmfs1/gscratch/scrubbed/retina/data/raw/',experiment_name,'/'];
        outputDirectory = ['/mmfs1/gscratch/scrubbed/retina/data/sorted/',experiment_name,'/'];
        tempDirectory = ['/mmfs1/gscratch/scrubbed/retina/data/sorted/',experiment_name,'/'];  
        txt = ['/mmfs1/gscratch/scrubbed/retina/data/sorted/',experiment_name,'/',file_name, '.csv'];
    end

    obj = LitkeToKilosort(dataDirectory, outputDirectory, tempDirectory);
    
    obj.setOptions();
    kilosortOptions = obj.kilosortOptions;
    kilosortOptions.NchanTOT = total_channels;
    kilosortOptions.sigmaMask = variance_distance;
    
    if contains(hostname,'mike','IgnoreCase',true) || contains(hostname,'maverick','IgnoreCase',true)
        kilosortOptions.chanMap = ['~/Documents/git_repos/manookin-lab/MEA/src/pipeline_utilities/kilosort/',map_name];
    elseif contains(hostname,'obsidian','IgnoreCase',true) 
        kilosortOptions.chanMap = ['/home/mike/Documents/git_repos/manookin-lab/MEA/src/pipeline_utilities/kilosort/',map_name];
    else
        kilosortOptions.chanMap = ['/mmfs1/gscratch/retina/GitRepos/manookin-lab/MEA/src/pipeline_utilities/kilosort/',map_name];
    end

    % Get the kilosort options
%     kilosortOptions = setOptions(folder_path);
    kilosortOptions.fbinary = [folder_path, file_name, '.bin'];
    kilosortOptions.fproc = [folder_path, file_name, 'temp_wh.dat'];

    results = preprocessDataSub(kilosortOptions);
    
    results = clusterSingleBatches(results);
    
    results = learnAndSolve8b(results);
    
    results = find_merges(results, 1);

    % final splits by SVD
    results = splitAllClusters(results, 1);

    % final splits by amplitudes
    results = splitAllClusters(results, 0);
    
    % decide on cutoff
    results = set_cutoff(results);
    
    xy = obj.getXYSpikeCoordinates(results);
                
    results.xy = xy;

    fprintf('found %d good units \n', sum(results.good>0))
    
    % Save out the results.
    obj.saveSingleData(results, file_name);
end

function file_names = parseMetaTxt(txt)
    fileID = fopen(txt, 'r');
    tline = fgetl(fileID);
    fclose(fileID);
    file_names = strsplit(tline,' ');
end

function [file_names, file_numbers, durations] = parseMetaCSV(txt)
    file_names = {};
    file_numbers = [];
    durations = [];
    fileID = fopen(txt, 'r');
    while true
      tline = fgetl(fileID);
      if ~ischar(tline) 
          break; 
      end
      
      splits = strsplit(tline,',');
      durations = [durations, str2double(splits{2})]; %#ok<AGROW>
      splits = strsplit(splits{1},'/');
      fname = splits{end};
      file_names = [file_names,fname]; %#ok<AGROW>
      file_numbers = [file_numbers,str2double(strrep(fname,'data',''))]; %#ok<AGROW>
    end
    fclose(fileID);
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
    ip.addParameter('NT', 64*1024 + 64, @(x)isfloat(x)); % Batch size; Must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory).
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
