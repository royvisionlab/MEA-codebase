function [out1,out2,out3] = getFastNoiseFrames(numXStixels, numYStixels, numXChecks, numYChecks, chromaticClass, numFrames, stepsPerStixel, seed)
% 

% Seed the random number generator.
noiseStream = RandStream('mt19937ar', 'Seed', seed);

% Generate the larger grid of stixels.
if strcmpi(chromaticClass, 'BY')
    gridValues = 2*(noiseStream.rand(numYStixels,numXStixels,numFrames*2) > 0.5)-1;
    % Replicate/expand the grid along the spatial dimensions.
    fullGrid = zeros(numYStixels*stepsPerStixel,numXStixels*stepsPerStixel,numFrames*2);
elseif strcmpi(chromaticClass, 'RGB')
    gridValues = 2*(noiseStream.rand(numYStixels,numXStixels,numFrames*3) > 0.5)-1;
    % Replicate/expand the grid along the spatial dimensions.
    fullGrid = zeros(numYStixels*stepsPerStixel,numXStixels*stepsPerStixel,numFrames*3);
else
    gridValues = 2*(noiseStream.rand(numYStixels,numXStixels,numFrames) > 0.5)-1;
    % Replicate/expand the grid along the spatial dimensions.
    fullGrid = zeros(numYStixels*stepsPerStixel,numXStixels*stepsPerStixel,numFrames);
end

for k = 1 : numYStixels*stepsPerStixel
    yindex = ceil(k/stepsPerStixel);
    for m = 1 : numXStixels*stepsPerStixel
        xindex = ceil(m/stepsPerStixel);
        fullGrid(k,m,:) = gridValues(yindex,xindex,:);
    end
end

% Generate the motion trajectory of the larger stixels.
xyStream = RandStream('mt19937ar', 'Seed', seed); % reseed
steps = round((stepsPerStixel-1)*xyStream.rand(2,numFrames));
xSteps = steps(1,:);
ySteps = steps(2,:);

% Get the frame values for the finer grid.
if strcmpi(chromaticClass, 'BY')
    frameValues = zeros(numYChecks,numXChecks,numFrames*2);
    for k = 1 : numFrames*2
        frameValues(:,:,k) = fullGrid((1:numYChecks)+ySteps(ceil(k/2)),(1:numXChecks)+xSteps(ceil(k/2)),k);
    end
    out1 = frameValues(:,:,1:2:end); % Yellow
    if nargout > 1
        out2 = frameValues(:,:,2:2:end); % Blue
        out3 = [];
    else
        out2 = [];
        out3 = [];
    end
elseif strcmpi(chromaticClass, 'RGB')
    frameValues = zeros(numYChecks,numXChecks,numFrames*3);
    for k = 1 : numFrames*3
        frameValues(:,:,k) = fullGrid((1:numYChecks)+ySteps(ceil(k/3)),(1:numXChecks)+xSteps(ceil(k/3)),k);
    end
    out1 = frameValues(:,:,1:3:end);
    if nargout == 1
        out2 = [];
        out3 = [];
    elseif nargout == 2
        out2 = frameValues(:,:,2:3:end);
        out3 = [];
    else
        out2 = frameValues(:,:,2:3:end);
        out3 = frameValues(:,:,3:3:end);
    end
else
    frameValues = zeros(numYChecks,numXChecks,numFrames);
    for k = 1 : numFrames
        frameValues(:,:,k) = fullGrid((1:numYChecks)+ySteps(k),(1:numXChecks)+xSteps(k),k);
    end
    out1 = frameValues;
    if nargout > 1
        out2 = [];
        out3 = [];
    end
end


