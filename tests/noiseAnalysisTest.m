clearvars;

javaaddpath('../src/vision7_symphony/Vision.jar');

% import java.io.*;
% import java.net.*;
% import java.util.ArrayList;
% import java.util.Collections;
% import java.util.List;
import edu.ucsc.neurobiology.vision.util.*;

sampleRate = 2e4;
binRate = 60;

numXChecks = 80;
numYChecks = 60;
numXStixels = 41;
numYStixels = 31;
preTime = 250;
stimTime = 19500;
tailTime = 250;
frameDwell = 1;
numFrames = 1185;
stepsPerStixel = 2;
chromaticClass = 'RGB';
monitorChannel = 0;

% targetElectrode = 161; threshold = -70;
% targetElectrode = 252; threshold = -70;
targetElectrode = 114; threshold = -70;

% targetBlocks = 34; seeds{1} = [1482380894,1482386986,1482393056,1482399130,1482405196,1482411267];

targetBlocks = 10:23;
seeds{1} = [1481746840, 1481752947, 1481759022];
seeds{2} = [1481769234, 1481775342, 1481781423];
seeds{3} = [1481793244, 1481799347, 1481805422];
seeds{4} = [1481816130, 1481822249, 1481828324, 1481834404, 1481840489, 1481846565];
seeds{5} = [1481857637, 1481863758, 1481869842, 1481875919, 1481881994, 1481888070];
seeds{6} = [1481898980, 1481905094, 1481911171, 1481917261, 1481923332, 1481929413];
seeds{7} = [1481941091, 1481947194, 1481953268, 1481959350, 1481965430, 1481971505];
seeds{8} = [1481982492, 1481988594, 1481994662, 1482000749, 1482006815, 1482012888];
seeds{9} = [1482024422, 1482030519, 1482036601, 1482042668, 1482048744, 1482054815];
seeds{10} =[1482063893, 1482070002, 1482076077, 1482082159, 1482088226, 1482094297];
seeds{11} =[1482106064, 1482112174, 1482118246, 1482124329, 1482130400, 1482136484];
seeds{12} =[1482145401, 1482151528, 1482157607, 1482163683, 1482169762, 1482175844];
seeds{13} =[1482186788, 1482192888, 1482198960, 1482205038, 1482211110, 1482217194];
seeds{14} =[1482226720, 1482232826, 1482238907, 1482245005, 1482251076, 1482257144];

t = tic;

STA = zeros(numYChecks,numXChecks,20);

for ii = 1 : length(targetBlocks)

datadir = ['/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Data/20220406C/data0',num2str(targetBlocks(ii)),'/'];

listing = dir(datadir);


h = edu.ucsc.neurobiology.vision.io.RawDataFile(datadir);
header = h.getHeader();


startSample = 0;
% Get the frame monitor
frameMonitor = -h.getData(monitorChannel, startSample, header.getNumberOfSamples());

% Get the frame times.
[frameTimes, epochStarts] = getFrameTimesMEA(frameMonitor);
disp(['Pulled the frame times in ',num2str(toc(t)),' seconds']);

%% Pull data from the target electrode.

h = edu.ucsc.neurobiology.vision.io.RawDataFile(datadir);
d = h.getData(targetElectrode, startSample, header.getNumberOfSamples());

% Wavelet filter.
w = wavefilter(d(:)', 6);

% Flip over the data for negative thresholds.
if threshold < 0
    spDir = -1;
else
    spDir = 1;
end

% Get the threshold crosses.
spikeTimes = getThresCross(w, threshold, spDir);

spikes = zeros(size(w));
spikes(spikeTimes) = 1;

disp(['Detected spikes in ',num2str(toc(t)),' seconds']);

%% Get the stimulus.


preFrames = round(preTime*1e-3 * 60);
stimFrames = round(stimTime*1e-3 * 60); %numFrames; round(stimTime*1e-3 * 60);


for jj = 1 : length(epochStarts)
    if jj < length(epochStarts)
        whichFrames = find(frameTimes >= epochStarts(jj) & frameTimes < epochStarts(jj+1));
    else
        whichFrames = find(frameTimes >= epochStarts(jj));
    end
    % Get the frames when the stimulus was on.
    fTimes = frameTimes(whichFrames(preFrames+(1:stimFrames+1)));
    
    [Y, B] = getFastNoiseFrames(numXStixels, numYStixels, numXChecks, numYChecks, chromaticClass, stimFrames, stepsPerStixel, seeds{ii}(jj));
    
    % Deal with the frame drops.
%     actualFrames = ceil((fTimes(end)-fTimes(1))/sampleRate * 60);
%     Y = zeros([size(ytmp,1:2),actualFrames]);
    
    response = zeros(1,stimFrames);
    for kk = 1 : length(fTimes)-1
        response(kk) = sum(spikes(fTimes(kk):fTimes(kk+1)-1));
    end
    
    for kk = 1 : numYChecks
        for mm = 1 : numXChecks
            tmp = ifft( fft([zeros(1,10),response,zeros(1,10)]) .* ...
                conj(fft([zeros(1,10),squeeze(Y(kk,mm,:))',zeros(1,10)])) );
            STA(kk,mm,:) = squeeze(STA(kk,mm,:)) + real(tmp(1:20))';
        end
    end
    
end
end

figure(1); clf;
imagesc(mean(STA(:,:,2:3),3));
axis image;

return

[ftimes, frate] = manookinlab.ovation.getFrameTimesFromMonitor(frameMonitor', sampleRate, binRate);


% foo = h.getData(startSample, header.getNumberOfSamples());
foo = h.getData(startSample, min(40000,header.getNumberOfSamples()));

plot(foo(:,162))





