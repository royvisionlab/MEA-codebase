function [epochStarts, epochEnds] = PullEpochTriggers(experimentName)

monitorChannel = 0;

hostname = getenv('HOSTNAME');
if isempty(hostname)
    [~, hostname] = system('hostname');
end
if contains(hostname,'hyak','IgnoreCase',true) || contains(hostname,'klone','IgnoreCase',true)
    binDir = ['/mmfs1/gscratch/scrubbed/retina/data/raw/',experimentName,'/'];
    outdir = '/mmfs1/gscratch/retina/data/raw/';
    javaaddpath('/gscratch/retina/GitRepos/manookin-lab/MEA/src/vision7_symphony/Vision.jar');
elseif contains(hostname,'obsidian','IgnoreCase',true)
    binDir = ['/data/data/raw/',experimentName,'/'];
    outdir = '/data/';
else
    binDir = ['/usr/share/pool/data/raw/',experimentName,'/'];
    outdir = '/home/mike/ftp/files/';
end


epochStarts = containers.Map;
epochEnds = containers.Map;
numSamples = containers.Map;

% Get the directory listing.
listing = dir(binDir);

for jj = 1 : length(listing)
    if contains(listing(jj).name, 'data')
        disp(['Processing file ',listing(jj).name,'...']);
        % Read the header
        datadir = [binDir, listing(jj).name];
        h = edu.ucsc.neurobiology.vision.io.RawDataFile(datadir);
        header = h.getHeader();
        sampleRate = header.getSamplingFrequency();
        arrayId = header.getArrayID();
        
        % Pull the monitor frames.
        frameMonitor = -h.getData(monitorChannel, 0, header.getNumberOfSamples());
        
        % Get the frame times.
        try
            [starts, ends] = getFrameTimesMEA(frameMonitor);
        catch
            starts = [];
            ends = [];
        end
        epochStarts(listing(jj).name) = starts;
        epochEnds(listing(jj).name) = ends;
        numSamples(listing(jj).name) = header.getNumberOfSamples();
    end
end
    
save([outdir,experimentName,'_triggers.mat'],'epochEnds','epochStarts','sampleRate','arrayId','numSamples');



