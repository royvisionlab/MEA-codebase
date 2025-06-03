classdef LitkeToKilosort < SpikeSorting
    properties
        
    end

    properties (SetAccess = private)
        kilosortVersion
        kilosortOptions % Spike sorting options for Kilosort
    end

    methods
        function obj = LitkeToKilosort(dataDirectory, outputDirectory, tempDirectory, varargin)
            % Call the superclass on initialization.
            obj = obj@SpikeSorting(dataDirectory, outputDirectory, tempDirectory);
            
            ip = inputParser();
            ip.addParameter('version', 2, @(x)isfloat(x)); % Version of kilosort to run.
            ip.parse(varargin{:});
            
            if (ip.Results.version == 3)
                obj.kilosortVersion = 3;
            elseif (ip.Results.version == 2.5)
                obj.kilosortVersion = 2.5;
            else
                obj.kilosortVersion = 2;
            end
            

            % Set the default options for running Kilosort3.
            obj.setOptions(); % User can change these by passing arguments to method.
        end
        
        function setVersion(obj, version)
            if (version == 3)
                obj.kilosortVersion = 3;
            elseif version == 2.5
                obj.kilosortVersion = 2.5;
            else
                obj.kilosortVersion = 2;
            end
        end
        
        
        % This method runs Kilosort on a specified .bin file.
        function results = run(obj, dirName, varargin)
            ip = inputParser;
            ip.addParameter('saveResults', false, @(x)islogical(x));
            % Parse the inputs.
            ip.parse(varargin{:});
            
            obj.saveResults = ip.Results.saveResults;
            
            % Save the metadata.
%             if obj.saveResults
%                 outdir = fullfile(fileparts(obj.outputDirectory),dirName);
%                 if ~exist(outdir,'dir')
%                     mkdir(outdir);
%                 end
%                 files = struct();
%                 files.fileNumbers = obj.fileNumbers;
%                 files.fileNames = obj.fileNames;
%                 files.timeSamples = obj.timeSamples;
%                 save([outdir,filesep,'metadata.mat'],'-v7.3','files');
%             end

            % Set the binary file.
            obj.kilosortOptions.fbinary = fullfile(fileparts(obj.tempDirectory),dirName,[dirName,'.bin']);
            
            % If the GPU device has not been selected, do it.
            if ~obj.gpuInitiated
                obj.setGPUDevice();
                obj.gpuInitiated = true;
            end

            if obj.kilosortVersion == 3
                results = obj.runVersion3(dirName);
            elseif obj.kilosortVersion == 2.5
                results = obj.runVersion25(dirName);
            else
                results = obj.runVersion2(dirName);
            end
            
        end
        
        function results = runVersion2(obj, dirName)
            % OK. Now, we'll call native kilosort functions.
            results = preprocessDataSub(obj.kilosortOptions);
            results = clusterSingleBatches(results);
            
            % They save shit here...
%             if (ip.Results.outputResults)
%                 % If the output directory does not exist, create it.
%                 outdir = fullfile(fileparts(obj.outputDirectory),dirName,'kilosort2');
%                 if ~exist(outdir,'dir')
%                     disp(['Output directory ',outdir,' does not exist, creating it...'])
%                     mkdir(outdir);
%                 end
%                 save([outdir,filesep,'results.mat'],'results');
%             end
            
            results = learnAndSolve8b(results);
            
            results = find_merges(results, 1);
            
            % final splits by SVD
            results = splitAllClusters(results, 1);
            
            % final splits by amplitudes
            results = splitAllClusters(results, 0);

            % decide on cutoff
            results = set_cutoff(results);

            fprintf('found %d good units \n', sum(results.good>0))
            
            if (obj.saveResults)
                % If the output directory does not exist, create it.
                %outdir = fullfile(fileparts(obj.outputDirectory),dirName,'kilosort2');
                %if ~exist(outdir,'dir')
                    %disp(['Output directory ',outdir,' does not exist, creating it...'])
                    %mkdir(outdir);
                %end
                
                
                % Get the xy coordinates of each spike.
                xy = obj.getXYSpikeCoordinates(results);
                
                results.xy = xy;
                
%                 results.cProj = [];
%                 results.cProjPC = [];
                
                if length(obj.fileNumbers) > 1
                    obj.saveBatchData(results);
                end

                % Save the full results no matter what.
                obj.saveSingleData(results, dirName);
% 
%                 % save final results as results2
%                 fprintf('Saving final results in results2  \n')
%                 fname = fullfile(rootZ, 'results2.mat');
%                 save(fname, 'results', '-v7.3');
                
                % Save a mat file containing the results.
                try
                    %results = gather(results);
%                     save([outdir,filesep,dirName,'.mat'],'results');
                    %save([outdir,filesep,'results2.mat'],'-v7.3','results');
                catch
                end
            end
        end
        
        function results = runVersion25(obj, dirName)
            % OK. Now, we'll call native kilosort functions.
            results = preprocessDataSub(obj.kilosortOptions);
            results = datashift2(results, 1);
            
            iseed = 1;
            results = learnAndSolve8b(results, iseed);
            
            results = find_merges(results, 1);
            
            % final splits by SVD
            results = splitAllClusters(results, 1);
            
            % final splits by amplitudes
%             results = splitAllClusters(results, 0);

            % decide on cutoff
            results = set_cutoff(results);
            % eliminate widely spread waveforms (likely noise)
            results.good = get_good_units(results);

            fprintf('found %d good units \n', sum(results.good>0))
            
            if (obj.saveResults)
                
                % Get the xy coordinates of each spike.
                xy = obj.getXYSpikeCoordinates(results);
                
                results.xy = xy;
                
                if length(obj.fileNumbers) > 1
                    obj.saveBatchData(results);
                end

                % Save the full results no matter what.
                obj.saveSingleData(results, dirName);
            end
        end
        
        function results = runVersion3(obj, dirName)
            results                = preprocessDataSub(obj.kilosortOptions);
            results                = datashift2(results, 1);
            
            [results, st3, tF]     = extract_spikes(results);
            
            results                = template_learning(results, tF, st3);
            
            [results, st3, tF]     = trackAndSort(results);
            
            results                = final_clustering(results, tF, st3);
            
            results                = find_merges(results, 1);

            % Output the results file.
            if (obj.saveResults)
                % If the output directory does not exist, create it.
                outdir = fullfile(fileparts(obj.outputDirectory),dirName,'kilosort3');
                if ~exist(outdir,'dir')
                    disp(['Output directory ',outdir,' does not exist, creating it...'])
                    mkdir(outdir);
                end
                % Save the sorting results in NPY format
                resultsToPhy2(results, outdir);
                
                % Get the xy coordinates of each spike.
                xy = gather(results.xy);
                writeNPY(xy, [outdir,filesep,'spike_xy.npy']);
                

                % save final results as results2
                fprintf('Saving final results in results2  \n')
                fname = fullfile(rootZ, 'results2.mat');
                save(fname, 'results', '-v7.3');
                
                % Save a mat file containing the results.
                try
                    results = gather(results);
%                     save([outdir,filesep,dirName,'.mat'],'results');
                    save([outdir,filesep,'results2.mat'],'-v7.3','results');
                catch
                end
                
            end
        end
        
        function saveSingleData(obj, rez, dirName)
            outdir = fullfile(fileparts(obj.outputDirectory),dirName,['kilosort',num2str(obj.kilosortVersion)]);
            % Create the output directory if it doesn't exist.
            if ~exist(outdir, 'dir')
                mkdir(outdir);
            end
            
            [~, isort]   = sort(rez.st3(:,1), 'ascend');
            rez.st3      = rez.st3(isort, :);
            rez.cProj    = rez.cProj(isort, :);
            rez.cProjPC  = rez.cProjPC(isort, :, :);
            rez.xy = rez.xy(isort,:);
            
            % Save the sorting results in NPY format.
            rezToPhy(rez, outdir);
            writeNPY(rez.xy, [outdir,filesep,'spike_xy.npy']);
        end
        
        function saveBatchData(obj, results)
            copyParams = {'good','ops','W','Wphy','U','mu','xcoords','ycoords','iNeigh','iNeighPC','Wrot','est_contam_rate','simScore'};
            copyParams = unique(copyParams);
            
            % Sort the spike times vector.
            [~, isort]      = sort(results.st3(:,1), 'ascend');
            results.st3     = results.st3(isort, :);
            results.cProj   = results.cProj(isort, :);
            results.cProjPC = results.cProjPC(isort, :, :);
            results.xy      = results.xy(isort,:);
            
            % Loop through each of the files and pull out the appropriate
            % data.
            for jj = 1 : length(obj.fileNames)
                % Get the time index.
                if (jj == 1)
                    idx = [1, obj.timeSamples(1)];
                else
                    idx = sum(obj.timeSamples(1:jj-1)) + [1, obj.timeSamples(jj)];
                end
                % Find the spike times that fall in the file window.
                spIdx = (results.st3(:,1) >= idx(1) & results.st3(:,1) <= idx(end));
                if sum(spIdx) > 0
                    
                    outdir = fullfile(fileparts(obj.outputDirectory),obj.fileNames{jj},['kilosort',num2str(obj.kilosortVersion)]);
                    % Create the output directory if it doesn't exist.
                    if ~exist(outdir, 'dir')
                        mkdir(outdir);
                    end

                    % Create a results structure for this data file.
                    rez = struct();
                    % Copy important parameters.
                    for pp = 1 : length(copyParams)
                        rez.(copyParams{pp}) = gather(results.(copyParams{pp}));
                    end
                    
                    rez.st3 = results.st3(spIdx,:);
                    rez.cProj = results.cProj(spIdx,:);
                    rez.cProjPC = results.cProjPC(spIdx,:,:);
                    rez.xy = results.xy(spIdx,:);
                    % Correct for the spike timing due to concatenation.
                    rez.st3(:,1) = rez.st3(:,1) - idx(1) + 1;

                    % Save the sorting results in NPY format.
                    rezToPhy(rez, outdir);
                    writeNPY(rez.xy, [outdir,filesep,'spike_xy.npy']);
                end
            end
        end

        function setOptions(obj, varargin)
            % Parse the input options
            ip = inputParser();
            ip.addParameter('trange', [0, Inf], @(x)isfloat(x)); % Time range to sort
            ip.addParameter('NchanTOT', 512, @(x)isfloat(x)); % Total number of channels in your recording
            ip.addParameter('fproc', fullfile(obj.outputDirectory,'temp_wh.dat'), @(x)ischar(x)); % WTF is this?
            ip.addParameter('sig', 20, @(x)isfloat(x)); % Spatial smoothness constant for registration
            ip.addParameter('nblocks', 5, @(x)isfloat(x)); % Blocks for registration. 0 turns it off, 1 does rigid registration. Replaces "datashift" option.

            ip.addParameter('chanMap', ['.',filesep,'kilosort',filesep,'chanMap.mat'], @(x)ischar(x)); % Path to channel configuration .mat file.
            ip.addParameter('fs', 20000, @(x)isfloat(x)); % Sample rate during acquisition
            ip.addParameter('fshigh', 300, @(x)isfloat(x)); % High-pass filtering frequency in Hz
            ip.addParameter('minfr_goodchannels', 0.1, @(x)isfloat(x)); % Minimum firing rate on a "good" channel (0 to skip)
            ip.addParameter('Th', [10, 4], @(x)isfloat(x)); % Threshold on projections (like in Kilosort1, can be different for last pass like [10 4]). Default was [9,9]...
            ip.addParameter('lam', 10, @(x)isfloat(x)); % How important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot)
            ip.addParameter('AUCsplit', 0.92, @(x)isfloat(x));% was 0.97 % Splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
            ip.addParameter('minFR', 1/50, @(x)isfloat(x)); % Minimum spike rate (Hz), if a cluster falls below this for too long it gets removed
            ip.addParameter('momentum', [20, 400], @(x)isfloat(x)); % Number of samples to average over (annealed from first to second value)
            ip.addParameter('sigmaMask', 70, @(x)isfloat(x));% was 70 % Spatial constant in um for computing residual variance of spike
            ip.addParameter('ThPre', 8, @(x)isfloat(x)); % Threshold crossings for pre-clustering (in PCA projection space)
            ip.addParameter('CAR',true, @(x)islogical(x)); % Whether to use common average referencing.
            ip.addParameter('max_peels', 50, @(x)isfloat(x)); % Maximum number of spikes to "peel" from the data around each threshold crossing

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
            obj.kilosortOptions = struct();
            for jj = 1 : length(fnames)
                obj.kilosortOptions.(fnames{jj}) = ip.Results.(fnames{jj});
            end
        end

        function xy = getXYSpikeCoordinates(obj, results)
%             wPCA  = results.ops.wPCA;
%             wTEMP  = results.ops.wTEMP;

            tF = results.cProjPC;

            ktid = int32(results.st3(:,2));

            iC = getClosestChannels(results, results.ops.sigmaMask, min(results.ops.Nchan, 32));

            [~,iW] = max(abs(results.dWU(results.ops.nt0min,:,:)),[],2);
            iW = int32(squeeze(iW));

            iC = gather(iC(:,iW));

            uweigh = abs(results.U(:,:,1));
            uweigh = uweigh ./ sum(uweigh,1);
            ycup = sum(uweigh .* results.yc, 1);
            xcup = sum(uweigh .* results.xc, 1);
            
            dmin = median(diff(unique(results.yc)));

            yunq = unique(results.yc);
            mxc = zeros(numel(yunq), 1);
            for j = 1:numel(yunq)
                xc = results.xc(results.yc==yunq(j));
                if numel(xc)>1
                   mxc(j) = median(diff(sort(xc))); 
                end
            end
            dminx = max(5, median(mxc));

            ycenter = (min(results.yc) + dmin-1):(2*dmin):(max(results.yc)+dmin+1);
            xcenter = (min(results.xc) + dminx-1):(2*dminx):(max(results.xc)+dminx+1);
            [xcenter, ycenter] = meshgrid(xcenter, ycenter);
            xcenter = xcenter(:);
            ycenter = ycenter(:);
            
            nst = numel(ktid);

            xy = zeros(nst, 2);

            for j = 1:numel(ycenter)
                y0 = ycenter(j);
                x0 = xcenter(j);    
                xchan = (abs(ycup - y0) < dmin) & (abs(xcup - x0) < dminx);

                itemp = find(xchan);

                tin = ismember(ktid, itemp);

                if sum(tin)<1
                    continue;
                end

                pid = ktid(tin);
                data = tF(tin, :, :);

                ich = unique(iC(:, itemp));
            %     ch_min = ich(1)-1;
            %     ch_max = ich(end);

                nsp = size(data,1);
                dd = zeros(nsp, 3,  numel(ich),  'single');
            %     dd = zeros(nsp, 6,  numel(ich),  'single');
                for k = 1:length(itemp)
                    ix = pid==itemp(k);
                    [~,ia,ib] = intersect(iC(:,itemp(k)), ich);
                    dd(ix, :, ib) = data(ix,:,ia);
                end
                xy(tin, :) =  obj.spikePosition(dd, results.xc(ich), results.yc(ich));
            end
        end
        
        
    end
    
    methods (Static)
        
        
        function xy = spikePosition(dd, yc, xc)
            dT = gpuArray(dd);
            [nspikes, nPC, nchan] = size(dT);
            [~, imax] = max(max(dT.^2, [], 3), [], 2);
            dBest = gpuArray.zeros(nspikes, nchan, 'single');
            for j = 1:nPC
                iX = imax==j;
                dBest(iX, :) = dT(iX, j, :);
            end

            dBest = max(0, dBest);
            dBest = dBest ./ sum(dBest,2);

            ysp = dBest * yc;
            xsp = dBest * xc;

            xy = [xsp, ysp];

            xy = gather(xy);
        end
    end
end
