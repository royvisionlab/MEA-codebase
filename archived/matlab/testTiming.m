

h = edu.ucsc.neurobiology.vision.io.RawDataFile('/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Data/20220420C/data034');
header = h.getHeader();
f = h.getData(0,0, header.getNumberOfSamples());



e = zeros(size(f)); e(epochStarts)=2200;

figure(1); clf;
hold on
plot(-f);
plot(e);
hold off;

%%
% Downsample to 1kHz.
f1 = f(1:20:end);

fr = [];
for jj = 1 : length(epochStarts)
    % Get the frame times in msec and convert to samples.
    frameTimes = block(bb).frameTimesMs{jj};
    frameTimes = frameTimes - frameTimes(1);

    % Get the frame timing relative to the clock (in msec).
    reTimes = epochStarts(jj)/sampleRate * 1e3 + frameTimes;

    targetFrames = preFrames+(1 : find(frameTimes <= preTime+stimTime, 1, 'last')+1);
    if length(targetFrames) > (numFrames+1)
        targetFrames = targetFrames(1: (numFrames+1));
    end

    % Get the frames when the stimulus was on.
    fTimes = round(reTimes(targetFrames));
    fr = [fr; fTimes];
end
t = zeros(size(f1));
t(fr) = 2200;

figure(2); clf;
hold on
plot(-f1)
plot(t);
hold off;

%%
c = [];
for gg = 1 : length(groupIdx)
    block = pp.group(groupIdx(gg)).block;

    for bb = 1 : length(block)
        dataFileName = strsplit(block(bb).dataFile,'/');
        dataFileName = dataFileName{end-1};
        obj.readNPYResults(dataFileName);
        clusterIdx = find(obj.spikeIdentity == 145);
        c = [c, clusterIdx];
    end
end
