function [epochStarts, epochEnds] = getFrameTimesMEA(fMonitor, transitionThreshold)

if nargin < 2
    transitionThreshold = 2000;
end

% Low-pass filter
fMonitor = lowPassFilter(fMonitor, 360, 5e-5);

fMonitor = fMonitor(:)';

fMonitor = fMonitor - min(fMonitor);
fMonitor = fMonitor ./ max(fMonitor);
ups = getThresCross(fMonitor, 0.5, 1);
downs = getThresCross(fMonitor, 0.5, -1);

% Check if we're running a file from before splitting out the frame flips
% from the epoch ttls...
frameTimes = round(sort([ups'; downs']));

if length(frameTimes) < 2
    epochStarts = ups;
    epochEnds = downs;
elseif mean(diff(frameTimes)) < 5000
    % Find the Epoch transitions.
    dtimes = [transitionThreshold; diff(frameTimes)];

    % Get the start index (largest transitions).
    startIdx = find(dtimes >= transitionThreshold);
    % Get the epoch end indeces.
    endIdx = [startIdx(2:end)-1; length(frameTimes)];
    epochStarts = frameTimes(startIdx);
    epochEnds = frameTimes(endIdx);
else
    epochStarts = ups;
    epochEnds = downs;
end



