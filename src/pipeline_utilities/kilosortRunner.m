function [results, obj] = kilosortRunner(experimentName, fileNumbers, outputName, device)

if ispc
	dataDirectory = ['E:\Data\rawdata\',experimentName,'\'];
	outputDirectory = ['E:\Data\kilosort\',experimentName,'\'];
    tempDirectory = ['C:\Users\Mike\Documents\KilosortTemp\',experimentName,'\'];
elseif ismac
    dataDirectory = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Data/20220420C/';
    outputDirectory = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Sort/20220420C/';
    tempDirectory = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Sort/20220420C/';
else
    %dataDirectory = '/home/mike/ftp/files/20220420C/';
    dataDirectory = ['/usr/share/pool/rawdata/MEA_Data/',experimentName,'/'];
    outputDirectory = ['/usr/share/pool/SortedData/',experimentName,'/']; %['/home/mike/ftp/files/ksort_out/',experimentName,'/'];
    tempDirectory = ['/usr/share/pool/SortedData/',experimentName,'/']; %['/home/mike/ftp/files/ksort_out/',experimentName,'/'];
end

if ~exist('device','var')
    device = 1:gpuDeviceCount;
end

if ischar(fileNumbers)
    listing = dir(dataDirectory);
    fileNumbers = [];
    for jj = 1 : length(listing)
        if contains(listing(jj).name, 'data')
        fileNumbers = [fileNumbers,str2double(strrep(listing(jj).name,'data',''))]; %#ok<AGROW>
        end
    end
end

t = tic;
obj = LitkeToKilosort(dataDirectory, outputDirectory, tempDirectory);

%dirName = 'data048';

% Check to see if the raw data files for sorting have been created.
outdir = fullfile(fileparts(obj.outputDirectory),'BinData');
if ~exist(outdir, 'dir')
    obj.litkeToBinaryBatch(false)
end


% Create the binary file.
% obj.litkeToBinary(dirName);
obj.setGPUDevice(device);

results = obj.runBatch(fileNumbers, outputName);

% See how long a single run takes.
toc(t)
