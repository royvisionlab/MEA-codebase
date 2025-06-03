function [frameTimes, epochStarts] = getFrameTimesMEA_old(fMonitor, transitionThreshold)

if nargin < 2
    transitionThreshold = 2000;
end

% Low-pass filter
fMonitor = lowPassFilter(fMonitor, 360, 5e-5);

fMonitor = fMonitor(:)';

fMonitor = fMonitor - min(fMonitor);
fMonitor = fMonitor./max(fMonitor);
ups = getThresCross(fMonitor,0.5,1);
downs = getThresCross(fMonitor,0.5,-1);

frameTimes = round(sort([ups'; downs']));

% Find the Epoch transitions.
dtimes = [transitionThreshold; diff(frameTimes)];

epochStarts = frameTimes(dtimes >= transitionThreshold);
