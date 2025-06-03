


% noiseStream = RandStream('mt19937ar','Seed',seed);
% 
% % Generate random initial positions for the dots.
% positions = ceil(noiseStream.rand(numDots,2) .* (ones(numDots,1)*screenSize));

stimFrames = 300;
motionPerFrame = 2;
spaceConstant = 5;
numDots = 1000;
screenSize = [40,30]*2;
seed = 1;
correlationFrames = 12;
splitContrasts = true;

tic;
M = getXYDotTrajectories(stimFrames,motionPerFrame,spaceConstant,numDots,screenSize,seed,correlationFrames,splitContrasts);
toc

%%
for j = 1 : stimFrames
    imagesc(M(:,:,j));
    drawnow;
end

