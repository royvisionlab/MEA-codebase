classdef SpikeAnalysis < handle
    
    properties
        dataDir
        spikeTimes
        spikeLocations
        spikeCenters
        spikeIdentity
    end
    
    methods
        function obj = SpikeAnalysis(dataDir)
            obj.dataDir = dataDir;
        end
        
        function readNPYResults(obj, dataFileName)
            resultsDir = fullfile(fileparts(obj.dataDir),dataFileName,'kilosort2',filesep);
            
            SPIKE_TIMES_FILENAME = 'spike_times.npy';
            SPIKE_IDENTITY_FILENAME = 'spike_clusters.npy';
            CLUSTER_QUALITY_FILENAME = 'cluster_KSLabel.tsv';
            SPIKE_LOCATIONS_FILENAME = 'spike_xy.npy';

            spike_times_vector = readNPY(fullfile(fileparts(resultsDir), SPIKE_TIMES_FILENAME));
            spike_identity_vector = readNPY(fullfile(fileparts(resultsDir), SPIKE_IDENTITY_FILENAME));
            quality_by_cluster_id = obj.buildClusterQualityDict(fullfile(fileparts(resultsDir), CLUSTER_QUALITY_FILENAME));
            spike_xy = readNPY(fullfile(fileparts(resultsDir), SPIKE_LOCATIONS_FILENAME));

            % Convert identity vector to Matlab indexing.
            spike_identity_vector = double(spike_identity_vector + 1);

            % Loop through each identified cell and get the spike times.
            goodCells = find(strcmp(quality_by_cluster_id, 'good'));
            if ~isempty(goodCells)
                obj.spikeTimes = cell(1, length(goodCells));
                obj.spikeLocations = cell(1, length(goodCells));
                obj.spikeCenters = zeros(length(goodCells),2);
                obj.spikeIdentity = zeros(1, length(goodCells));
                for jj = 1 : length(goodCells)
                    idx = (spike_identity_vector == goodCells(jj));
                    obj.spikeTimes{jj} = spike_times_vector(idx); % Grab the spike times for the cell.
                    obj.spikeLocations{jj} = spike_xy(idx,:); % Get the xy locations of the spikes on the array.
                    obj.spikeCenters(jj,:) = mean(spike_xy(idx,:), 1); % Get the xy center of the spike cluster.
                    obj.spikeIdentity(jj) = goodCells(jj); % Save the cluster number
                end
            end
        end
        
    end
    
    methods (Static)
        
        function qualityDict = buildClusterQualityDict(filepath)
            qualityDict = {};

            % Open the file for reading.
            fid = fopen(filepath, 'r');
            fgetl(fid); % Header line.
            tline = fgetl(fid);
            while ischar(tline)
                parts = strsplit(tline,'\t');
                % Add to the cell array
                qualityDict{str2double(parts{1})+1} = parts{2}; %#ok<AGROW>
                tline = fgetl(fid);
            end
            fclose(fid);
        end
        
    end
end