function [fTimes, frameRate] = getFrameTimes(frameMonitor, sampleRate, binRate)
%frameTimes in data points
%frameRate, mean over all flips, is frames/dataPoint

fRate = 59.994;
fTimes = cell(1, size(frameMonitor,1));
frameRate = zeros(1, size(frameMonitor,1));
for k = 1 : size(frameMonitor,1)
    fMonitor = frameMonitor(k,:);
    
    % Bandpass filter.
%     fMonitor = bandPassFilter(fMonitor,12,120,1/sampleRate);
    fMonitor = lowPassFilter(fMonitor,250,1/sampleRate);
    
    %shift & scale s.t. fr monitor lives on [0,1]
    fMonitor = fMonitor - min(fMonitor);
    fMonitor = fMonitor./max(fMonitor);
    
    ups = [1, getThresCross(fMonitor,0.5,1)];
    downs = getThresCross(fMonitor,0.5,-1);
    
    frameTimes = round(sort([ups'; downs']));

    if (sampleRate / mean(diff(frameTimes(3:end)))) > 110
        nframes = ceil(size(frameMonitor,2)/sampleRate*fRate);
        frameTimes = [1 sampleRate/fRate*ones(1,nframes-1)];
        frameTimes = round(cumsum(frameTimes));
%         skip = round(sampleRate / mean(diff(frameTimes(3:end))) / 60);
%         frameTimes = frameTimes(1 : skip : end);
    end
    
    % Get the frame rate.
    frameRate(k) = sampleRate / mean(diff(frameTimes));

    if binRate < sampleRate
        frameTimes = ceil(frameTimes*binRate/sampleRate);
    end
    fTimes{k} = frameTimes;
end
