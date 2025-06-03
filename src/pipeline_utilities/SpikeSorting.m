classdef SpikeSorting < handle
    properties
        gpuInitiated
        saveResults
    end
    
    properties (SetAccess = private)
        dataDirectory
        outputDirectory
        tempDirectory
        fileNumbers % Double file numbers to concatenate and batch process
        fileNames % String file names.
        timeSamples % Number of time samples per concatenated file.
    end
    
    methods
        function obj = SpikeSorting(dataDirectory, outputDirectory, tempDirectory)
            obj.dataDirectory = dataDirectory;
            obj.outputDirectory = outputDirectory;
            obj.tempDirectory = tempDirectory; % This should be a place on the main SSD.

            if ~strcmp(obj.dataDirectory(end),filesep)
                obj.dataDirectory = [obj.dataDirectory,filesep];
            end
            
            if ~strcmp(obj.outputDirectory(end),filesep)
                obj.outputDirectory = [obj.outputDirectory,filesep];
            end
            
            if ~strcmp(obj.tempDirectory(end),filesep)
                obj.tempDirectory = [obj.tempDirectory,filesep];
            end
            
            % Set the GPU device initiation as false.
            obj.gpuInitiated = false;
        end
        
        % Setter methods for properties
        function setDataDirectory(obj, value)
            if nargin > 1
                obj.dataDirectory = value;
            end
        end
        function set.dataDirectory(obj, value)
            if ischar(value)
                obj.dataDirectory = value;
            else
                error('Data directory must be a string/char.');
            end
        end
        
        function setFileNumbers(obj, value)
            if nargin > 1
                obj.fileNumbers = value;
            end
        end
        function set.fileNumbers(obj, value)
            obj.fileNumbers = value;
        end
        
        function setFileNames(obj, value)
            if nargin > 1
                obj.fileNames = value;
            end
        end
        function set.fileNames(obj, value)
            obj.fileNames = value;
        end
        
        
        function setTimeSamples(obj, value)
            if nargin > 1
                obj.timeSamples = value;
            end
        end
        function set.timeSamples(obj, value)
            obj.timeSamples = value;
        end
        
        function setOutputDirectory(obj, value)
            if nargin > 1
                obj.outputDirectory = value;
            end
        end
        
        function set.outputDirectory(obj, value)
            if ischar(value)
                obj.outputDirectory = value;
            else
                error('Output directory must be a string/char.');
            end
        end
        
        % Setter for the temporary directory
        function setTempDirectory(obj, value)
            if nargin > 1
                obj.tempDirectory = value;
            end
        end
        
        function set.tempDirectory(obj, value)
            if ischar(value)
                obj.tempDirectory = value;
            else
                error('Output directory must be a string/char.');
            end
        end
        
        function getFileLengths(obj, fileNumbers)
            fileNumbers = fileNumbers(:)';
            obj.fileNumbers = fileNumbers;
            obj.fileNames = cell(1, length(fileNumbers));
            obj.timeSamples = zeros(1, length(fileNumbers));
            for jj = 1 : length(fileNumbers)
                obj.fileNames{jj} = ['data',obj.fileNumToString(fileNumbers(jj))];
                
                % Get the number of time samples.
                obj.timeSamples(jj) = obj.getNumberOfTimeSamples(obj.fileNames{jj});
            end
        end
        
%         function results = run(obj, dirName, varargin)
%             results = [];
%         end
        
        % Batch process/concatenate files and run through the spike
        % sorting. Keep track of the file numbers and the number of time
        % samples from each file for reconstruction.
        function results = runBatch(obj, fileNumbers, outputName, concatBinaryFiles)
            if ~exist('concatBinaryFiles','var')
                concatBinaryFiles = true;
            end
            % Make sure the file numbers are single vector.
            fileNumbers = fileNumbers(:)';
            obj.fileNumbers = fileNumbers;
            obj.fileNames = cell(1, length(fileNumbers));
            obj.timeSamples = zeros(1, length(fileNumbers));
            for jj = 1 : length(fileNumbers)
                obj.fileNames{jj} = ['data',obj.fileNumToString(fileNumbers(jj))];
                
                % Get the number of time samples.
                obj.timeSamples(jj) = obj.getNumberOfTimeSamples(obj.fileNames{jj});
            end
            
            % Concatenate the raw binary files.
            if concatBinaryFiles
                obj.concatenateBinaryFiles(fileNumbers, outputName);
            end
            
            % Run spike sorting.
            results = obj.run(outputName);
        end
        
        % This method concatenates binary files.
        function concatenateBinaryFiles(obj, fileNumbers, outputName)
            if ~exist('outputName', 'var')
                outputName = 'tmp';
            end
            
            if ispc
                command = 'type ';
            else
                command = 'cat ';
            end
            
            for jj = 1 : length(fileNumbers)
                fileStr = ['data',obj.fileNumToString(fileNumbers(jj))];
                tmp = [obj.outputDirectory,'BinData',filesep,fileStr,'.bin'];
                command = [command, tmp, ' ']; %#ok<AGROW>
            end
            if ~exist([obj.outputDirectory,outputName,filesep],'dir')
                system(['mkdir ', obj.outputDirectory,outputName,filesep]);
            end
            if ~exist([obj.tempDirectory,outputName,filesep],'dir')
                system(['mkdir ', obj.tempDirectory,outputName,filesep]);
            end
            command = [command, '> ',[obj.tempDirectory,outputName,filesep],outputName,'.bin'];
            % Execute the command.
            status = system(command);
            disp(status);
        end
        
        
        function litkeToBinaryBatch(obj, runFilter)
            if ~exist('runFilter','var')
                runFilter = true;
            end
            listing = dir(obj.dataDirectory);
            % Loop through the directories.
            for jj = 1 : length(listing)
                if contains(listing(jj).name,'data') && listing(jj).isdir
                    disp(['Processing data file ',num2str(jj),' of ', num2str(length(listing))]);
                    try 
                        obj.litkeToBinary(listing(jj).name, runFilter);
                    catch
                    end
                end
            end
        end
        
        function timeSamples = getNumberOfTimeSamples(obj, dirName)
            % Open the Litke/Symphony data folder for reading.
            h = edu.ucsc.neurobiology.vision.io.RawDataFile(fullfile(fileparts(obj.dataDirectory),dirName));

            % Get the header.
            header = h.getHeader();
            
            timeSamples = double(header.getNumberOfSamples());
        end

        function litkeToBinary(obj, dirName, runFilter)
            if ~exist('runFilter','var')
                runFilter = true;
            end
            
            % Open the Litke/Symphony data folder for reading.
            h = edu.ucsc.neurobiology.vision.io.RawDataFile(fullfile(fileparts(obj.dataDirectory),dirName));

            % Get the header.
            header = h.getHeader();
            
            blockSize = 200000;
            numBlocks = ceil(header.getNumberOfSamples()/blockSize);
            
            % If the output directory does not exist, create it.
            outdir = fullfile(fileparts(obj.outputDirectory),'BinData'); % fullfile(fileparts(obj.outputDirectory),dirName);
            if ~exist(outdir,'dir')
                disp(['Output directory ',outdir,' does not exist, creating it...'])
                mkdir(outdir);
            end

            % Open a file for writing.
            fileID = fopen(fullfile(fileparts(obj.outputDirectory),'BinData',[dirName,'.bin']),'w');
            
            if runFilter
                if numBlocks == 1
                    data = h.getData(0, header.getNumberOfSamples());
                    data = data(:,2:end)';
                    % Wavelet filter the channels.
                    data = DB4Filter(data, 6);
                    fwrite(fileID,data,'int16','ieee-le');
                else
                    samplesWritten = 0; % Keep track of the number of samples written out.
                    for jj = 1 : numBlocks
%                         disp(['Processing data in block ',num2str(jj),' of ', num2str(numBlocks)]);
                        remSamples = header.getNumberOfSamples() - (jj-1)*blockSize;
                        lastSample = min(blockSize, remSamples);
                        startSample = (jj-1)*blockSize;
                        data = h.getData(startSample, lastSample);
                        data = data(:,2:end)';

                        % Determine the number of samples to write.
                        if jj == 1
                            wSamples = round(0.75 * blockSize);

                            lastIdx = round(0.5 * blockSize);
                            dBuffer = data(:,lastIdx+1:end);
                            samplesWritten = samplesWritten + (wSamples - size(dBuffer,2));
                            % Wavelet filter the channels.
                            data = DB4Filter(data, 6);
                            fwrite(fileID, data(:,1:wSamples), 'int16', 'ieee-le');
                        elseif jj == numBlocks
                            data = DB4Filter([dBuffer,data], 6);
                            clear dBuffer;
                            fwrite(fileID, data(:,samplesWritten+1:end), 'int16', 'ieee-le');
                        else
                            newD = [dBuffer,data];
                            dBuffer = data(:,lastIdx+1:end);
                            data = DB4Filter(newD, 6);
                            clear newD;
                            fwrite(fileID, data(:,samplesWritten+(1:blockSize)), 'int16', 'ieee-le');
                        end
                    end
                end
            else
                for jj = 1 : numBlocks
                    disp(['Processing data in block ',num2str(jj),' of ', num2str(numBlocks)]);
                    remSamples = header.getNumberOfSamples() - (jj-1)*blockSize;
                    lastSample = min(blockSize, remSamples);
                    startSample = (jj-1)*blockSize;
                    % foo = h.getData(startSample, header.getNumberOfSamples());
                    data = h.getData(startSample, lastSample);
                    data = data(:,2:end)';
                    

                    fwrite(fileID,data,'int16','ieee-le');
                end
            end
            fclose(fileID);
        end
        
        function setGPUDevice(obj, device)
            if ~exist('device','var')
                device = gpuDeviceCount;
            end
            
            if length(device) > 1
                % Check if the parpool is open.
                p = gcp('nocreate');
                if isempty(p)
                    parpool('local', numel(device));
                end
                spmd
                    gpuDevice(device(labindex));
                end
            else
                gpuDevice(device);
            end
            
            
            obj.gpuInitiated = true;
        end
        
        
    end
    
    methods (Static)
        
        % Wavelet filter the raw data. Columns are time points, rows are
        % channels.
        function fdata = waveletFilter(data, maxlevel)
            fdata = DB4Filter(data, maxlevel);
        end
        
        function fileStr = fileNumToString(fileNum)
            if fileNum < 10
                fileStr = ['00',num2str(fileNum)];
            elseif fileNum < 100
                fileStr = ['0',num2str(fileNum)];
            else
                fileStr = num2str(fileNum);
            end
        end
        
        function queryGPUDevices()
            for ii = 1:gpuDeviceCount
                g = gpuDevice(ii);
                fprintf(1,'Device %i has ComputeCapability %s \n', ...
                        g.Index,g.ComputeCapability)
            end
        end
    end
end