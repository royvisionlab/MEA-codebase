clearvars;


outdir = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/data/';

addpath(genpath('/Users/michaelmanookin/Google Drive/GitRepos/pillowlab/LNPfitting/'));
ovationDir = '/Users/michaelmanookin/Documents/Data/Ovation_MEA/';
binDir = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Data/';
sortDir = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Sort/';

javaaddpath('/Users/michaelmanookin/Documents/GitRepos/Manookin-Lab/MEA/src/vision7_symphony/Vision.jar');

colors = [
    0 0 255
    204 0 0
    0 0 0
    0 128 255
    0 128 0
    50 205 50
    204 0 0
    255 69 0
    ] / 255;

matName = 'MEA_Ovation';

nonlinearityBins = 100;

loader = edu.washington.rieke.Analysis.getEntityLoader(); 
treeFactory = edu.washington.rieke.Analysis.getEpochTreeFactory();

import auimodel.*
import vuidocument.*

% Epoch list
list = loader.loadEpochList([ovationDir,matName,'.mat'],ovationDir);

keywordSplitter = @(list)splitOnKeywords(list);
keywordSplitter_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, keywordSplitter);
dateSplit = @(list)splitOnExperimentDate(list);
dateSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, dateSplit);
blockSplit = @(list)splitOnBlockStart(list);
blockSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, blockSplit);
startSplit = @(list)splitOnEpochStart(list);
startSplit_java = riekesuite.util.SplitValueFunctionAdapter.buildMap(list, startSplit);

% tree = riekesuite.analysis.buildTree(list, {dateSplit_java,'cell.label','protocolSettings(led)','protocolSettings(stimTime)','protocolSettings(epochGroup:label)'});
% tree = riekesuite.analysis.buildTree(list, {dateSplit_java,'cell.label','protocolSettings(epochGroup:label)'});
% tree = riekesuite.analysis.buildTree(list, {dateSplit_java,'cell.label'});
% tree = riekesuite.analysis.buildTree(list, {'protocolID',dateSplit_java,'protocolSettings(epochGroup:label)'});
tree = riekesuite.analysis.buildTree(list, {dateSplit_java,'protocolID','protocolSettings(epochGroup:label)','protocolSettings(dataFileName)',startSplit_java});


gui = epochTreeGUI(tree);

%%
monitorChannel = 0;
node = gui.getSelectedEpochTreeNodes;
node = node{1};
 
triggersExist = false;
triggersLoaded = false;

experiment = struct();

% Protocol level.
for ii = 1 : length(node.children)
    prot = node.children(ii);
    
    % Get the protocol name.
    experiment.protocol(ii).label = prot.splitValue;
    
    % EpochGroup level
    for jj = 1 : length(prot.children)
        eGroup = prot.children(jj);
        
        experiment.protocol(ii).group(jj).label = eGroup.splitValue;

        % EpochBlock level.
        for kk = 1 : length(eGroup.children)
            eBlock = eGroup.children(kk);
            % Epoch level.
            epochList = eBlock.epochList;
            
            % Get the epoch start times.
            startTimes = zeros(1,length(epochList));
            start_s = zeros(1,length(epochList));
            for mm = 1 : length(epochList)
                foo = epochList.elements(mm).protocolSettings.get('epoch:startTime');
                startTimes(mm) = datenum([1970 1 1 0 0 foo.getTime/1000]);
                start_s(mm) = foo.getTime/1000; % Start time in seconds.
            end
            
            % Sort by start time.
            [~, epochOrder] = sort(startTimes);
            start_s = start_s(epochOrder);

            % Get the MEA data file name.
            filename = epochList.elements(1).protocolSettings.get('dataFileName');
            
%             disp(['block:', char(epochList.elements(1).protocolID), '...', eBlock.splitValue, '; epoch:', epochList.elements(1).protocolSettings.get('dataFileName')])
            
%             FlashedSpatialNoise...20220909C\data003.bin
%             FlashedSpatialNoise...20220909C\data000.bin
%             NaturalImageFlashPlusNoise...20220909C\data001.bin
%             NaturalImageFlashPlusNoise...20220909C\data002.bin
            if contains(char(epochList.elements(1).protocolID), 'FlashedSpatialNoise')
                if strcmp(filename, '20220909C\data000.bin')
                    filename = '20220909C\data026.bin';
                elseif strcmp(filename, '20220909C\data003.bin')
                    filename = '20220909C\data029.bin';
                end
            elseif contains(char(epochList.elements(1).protocolID), 'NaturalImageFlashPlusNoise')
                if strcmp(filename, '20220909C\data001.bin')
                    filename = '20220909C\data027.bin';
                elseif strcmp(filename, '20220909C\data002.bin')
                    filename = '20220909C\data028.bin';
                end
            end
            
            % Split the file name on \
            fname = strsplit(filename,'\');
            
            % Experiment from 20220426 was misnamed
            if strcmp(fname{1},'20200426C')
                fname{1} = '20220426C';
            elseif strcmp(fname{1}, '02230329C')
                fname{1} = '20230329C';
            end
            
            if ~contains(fname{1},'test') && ~isempty(fname{1}) && contains(fname{1},'202')
                binName = [fname{1},filesep,strrep(fname{2},'.bin',''),filesep];
                
                if contains(fname{1}, '202')
                    experimentName = fname{1};
                end
                
                if exist([outdir,'triggers/',experimentName,'_triggers.mat'],'file')
                    triggersExist = true;
                    if ~triggersLoaded
                        triggers = load([outdir,'triggers/',experimentName,'_triggers.mat']);
                        triggersLoaded = true;
                    end
                end

                %------------------------------------------------------------------
                datadir = [binDir,binName];
%                 disp([binName, ': ', char(epochList.elements(1).protocolID)]);
                
                
                if triggersExist
                    tmpName = strrep(fname{2},'.bin','');
                    epochStarts = triggers.epochStarts(tmpName);
                    epochEnds = triggers.epochEnds(tmpName);
                    sampleRate = 2e4;
                    arrayId = triggers.arrayId;
                    numSamples = triggers.numSamples(tmpName);
                    if isempty(epochStarts)
                        disp(['Warning! Epoch triggers empty for file ', binName]);
                        epochStarts = 46969 + round((start_s-start_s(1))*20000);
                    end
                else
                    % Read the header
                    
                    h = edu.ucsc.neurobiology.vision.io.RawDataFile(datadir);
                    header = h.getHeader();
                    sampleRate = header.getSamplingFrequency();

                    % Pull the monitor frames.
                    frameMonitor = -h.getData(monitorChannel, 0, header.getNumberOfSamples());
                    
                    % Get the frame times.
                    [epochStarts, epochEnds] = getFrameTimesMEA(frameMonitor);
                    
                    arrayId = header.getArrayID();
                    numSamples = header.getNumberOfSamples();
                end
                % Make sure this is not an empty file.
                experiment.protocol(ii).group(jj).block(kk).dataFile = datadir;
                % Use this for anything before May 2022 before stripping
                % out the frame monitor...
%                 [frameTimes, epochStarts] = getFrameTimesMEA_old(frameMonitor);

                % Make sure all of the epoch starts are accurate.
                if ~isempty(epochStarts)
                    epoch_diff = median(epochEnds - epochStarts(1:length(epochEnds)));
                    tmp = (start_s - start_s(1)) * sampleRate + epochStarts(1);
                    epochStarts = round(tmp);
                    epochEnds = epochStarts + epoch_diff;

                    unrecorded_epochs = (epochStarts > numSamples);       

                    % TTL directory
                    ttlDir = [sortDir,experimentName,'/TTLTriggers/'];
                    fileName = strrep(fname{2},'.bin','');
                    offset = 0;
                    saveTriggerPickle(ttlDir, fileName, epochStarts, numSamples, arrayId, offset);
                end

                % Get the frame monitor.
                if contains('Frame Monitor',cell(epochList.responseStreamNames()))
                    frameMonitor = riekesuite.getResponseMatrix(epochList, 'Frame Monitor');
                    [fTimes, frameRate] = getFrameTimes(frameMonitor, epochList.elements(1).protocolSettings.get('sampleRate'), 1000);
                    fTimes = fTimes(epochOrder);
                    frameRate = frameRate(epochOrder);
                else
                    fTimes = {};
                end
                
                % Get rid of epochs that don't have accompanying data.
                epochStarts(unrecorded_epochs) = [];
                epochEnds(unrecorded_epochs) = []; 
                fTimes(unrecorded_epochs) = []; 

                experiment.protocol(ii).group(jj).block(kk).epochStarts = epochStarts;
                experiment.protocol(ii).group(jj).block(kk).epochEnds = epochEnds;
                experiment.protocol(ii).group(jj).block(kk).frameTimesMs = fTimes; % Save the frame times in milliseconds

                for mm = 1 : length(epochList)
                    % Make sure the epoch played.
                    if unrecorded_epochs(mm) == 0
                        p = epochList.elements(epochOrder(mm)).protocolSettings;
                        keyNames = sort(p.keySet.toArray.cell);

                        % Get the epoch parameters.
                        parameters = containers.Map;
                        for keyIndex = 1 : length(keyNames)
                            foo = p.get(keyNames{keyIndex});

                            % Deal with special cases.
                            if isa(foo, 'java.util.ArrayList')
                                foo = convertArraylist(foo);
                            elseif isa(foo, 'java.util.Date')
                                foo = convertJavaDate(foo);
                            end
                            parameters(keyNames{keyIndex}) = foo;
                        end

                        experiment.protocol(ii).group(jj).block(kk).epoch(mm).parameters = parameters;
                    end
                end
            end
        end
    end
end

% Save the experiment structure to a mat file.
save([outdir,experimentName,'.mat'],'experiment');

% Save to a JSON file.
experimentToJSON(experiment, outdir, experimentName);

disp('All done!!');

% Utilities.
function mArray = convertArraylist(foo)
    iterator = foo.iterator;
    mArray = zeros(1,foo.size);
    count=0;
    while (iterator.hasNext)
        count=count+1;
        mArray(count) = double(iterator.next);
    end
end

function t = convertJavaDate(foo)
    sdf = java.text.SimpleDateFormat('yyyy-MM-dd''T''HH:mm:ss.SSSXXX');
    t = char(sdf.format(foo));
end

% dt = datetime(datenum([1970 1 1 0 0 foo.getTime/1000]),'ConvertFrom','datenum');
% javaDay = dt.Day;
% if javaDay < 10
%     javaDay = ['0', num2str(javaDay)];
% else
%     javaDay = num2str(javaDay);
% end
% javaMonth = dt.Month;
% if javaMonth < 10
%     javaMonth = ['0', num2str(javaMonth)];
% else
%     javaMonth = num2str(javaMonth);
% end
% cellName = [num2str(dt.Year), javaMonth, javaDay, rig(1), 'c', cellNum];