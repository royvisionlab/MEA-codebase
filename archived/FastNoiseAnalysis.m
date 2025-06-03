clearvars;

sampleRate = 2e4;

% experimentName = '20220420C'; protocolIdx = 3; groupIdx = 1:2;
experimentName = '20220426C'; protocolIdx = 3; groupIdx = 1:2;


if ispc
	dataDir = 'E:\Data\rawdata\';
	sortDir = 'E:\Data\kilosort\';
elseif ismac
    dataDir = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Data/';
    sortDir = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Sort/';
else
    dataDir = '/run/user/1000/gvfs/smb-share:server=maverick.local,share=sambashare/rawdata/MEA_Data/';
    sortDir = '/home/mike/ftp/files/ksort_out/';
end


staFrames = 30;

% Load the metadata file.
load([dataDir,experimentName,'.mat']);

return
obj = SpikeAnalysis([sortDir,experimentName,filesep]);


pp = experiment.protocol(protocolIdx);

% Get the total block count to keep track of progress
numBlocks = 0;
for gg = 1 : length(groupIdx)
    numBlocks = numBlocks + length(pp.group(groupIdx(gg)).block);
end

t = tic;
blockCount = 0;
% Loop through the blocks.
for gg = 1 : length(groupIdx)
    block = pp.group(groupIdx(gg)).block;

    for bb = 1 : length(block)
        disp(['Processing block ',num2str(blockCount+1),' of ',num2str(numBlocks)]);
        
        % Pull key stimulus parameters from the block.
        preTime = block(bb).epoch(1).parameters('preTime');
        stimTime = block(bb).epoch(1).parameters('stimTime');
        numXStixels = block(bb).epoch(1).parameters('numXStixels'); 
        numYStixels = block(bb).epoch(1).parameters('numYStixels');
        numXChecks = block(bb).epoch(1).parameters('numXChecks');
        numYChecks = block(bb).epoch(1).parameters('numYChecks');
        chromaticClass = block(bb).epoch(1).parameters('chromaticClass');
        numFrames = block(bb).epoch(1).parameters('numFrames');
        stepsPerStixel = block(bb).epoch(1).parameters('stepsPerStixel');
        dataFileName = strsplit(block(bb).dataFile,'/');
        dataFileName = dataFileName{end-1};

        preFrames = round(preTime*1e-3 * 60);
        stimFrames = round(stimTime*1e-3 * 60);

        epochStarts = block(bb).epochStarts';

        % Pull the spikes for the block.
        obj.readNPYResults(dataFileName);
        
        if blockCount == 0 % Initialize
            sY = zeros([length(obj.spikeTimes),60,80,staFrames]);
            sB = zeros([length(obj.spikeTimes),60,80,staFrames]);
        end
        blockCount = blockCount + 1;

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

            [Y, B] = getFastNoiseFrames(numXStixels, numYStixels, numXChecks, numYChecks, chromaticClass, numFrames, stepsPerStixel, block(bb).epoch(jj).parameters('seed'));

            Y = Y(:,:,1:length(fTimes)-1);
            Y(:,:,1:10) = 0;
            if ~isempty(B)
                B = B(:,:,1:length(fTimes)-1);
                B(:,:,1:10) = 0;
            end

            parfor electrode = 1 : length(obj.spikeTimes)
                % Spike times for the cluster in msec.
%                 clusterIdx = find(obj.spikeIdentity == 145);
                spikes = double(obj.spikeTimes{electrode})/sampleRate * 1e3; 

                response = zeros(1,length(fTimes)-1);
                for kk = 1 : length(fTimes)-1
                    response(kk) = sum(spikes >= fTimes(kk) & spikes < fTimes(kk+1)-1);
%                     response(kk) = sum(spikes(fTimes(kk):fTimes(kk+1)-1));
                end
                response(1:10) = 0;

                % Analyze the achromatic frames.
                if sum(response) > 0
                    for kk = 1 : numYChecks
                        for mm = 1 : numXChecks
                            tmp = fft([zeros(1,2),response,zeros(1,10)]) .* ...
                                conj(fft([zeros(1,2),squeeze(Y(kk,mm,:))',zeros(1,10)])) ;
                            tmp = real( ifft( tmp ) );
                            sY(electrode,kk,mm,:) = squeeze(sY(electrode,kk,mm,:)) + tmp(1:staFrames)';
                        end
                    end

                    if ~strcmp(chromaticClass, 'achromatic')
                        for kk = 1 : numYChecks
                            for mm = 1 : numXChecks
                                tmp = fft([zeros(1,2),response,zeros(1,10)]) .* ...
                                    conj(fft([zeros(1,2),squeeze(B(kk,mm,:))',zeros(1,10)])) ;
                                tmp = real( ifft( tmp ) );
                                sB(electrode,kk,mm,:) = squeeze(sB(electrode,kk,mm,:)) + tmp(1:staFrames)';
                            end
                        end
                    end
                end
            end
        end
        toc(t)
    end
end

%%
outdir = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/20220420C/';

f = figure(31); 
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'Position',[pos(1),pos(2),8.5,11]);
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[8.5, 11])

for jj = 1 : size(sY,1)
    sta = squeeze(sY(jj,:,:,:));
    [sfY, tfY] = getStaSummary(sta);
    sfY = sfY/std(sfY,[],'all');
    
    sta = squeeze(sB(jj,:,:,:));
    [sfB, tfB] = getStaSummary(sta);
    sfB = sfB / std(sfB,[],'all');
    
    % Scale the RFs by the same factor.
    sc = max(max(abs(sfY(:))),max(abs(sfB(:))));
    sfY = sfY / sc;
    sfB = sfB / sc;
    
    % Convert to 8-bit uint
    sfY = uint8(255 * (0.5*sfY + 0.5));
    sfB = uint8(255 * (0.5*sfB + 0.5));
    
    imgY = zeros([size(sfY),3],'uint8');
    imgY(:,:,1:2) = repmat(sfY,[1,1,2]);
    
    imgB = zeros([size(sfB),3],'uint8');
    imgB(:,:,3) = sfB;
    
    % Plot the results.
    clf(f);
    subplot(3,1,1);
    image(sfY); axis image;
    set(gca,'XTick',[]); set(gca,'YTick',[]);
    title(['Yellow Rf; Cluster: ',num2str(jj)]);
    
    subplot(3,1,2);
    image(sfB); axis image;
    set(gca,'XTick',[]); set(gca,'YTick',[]);
    title('Blue RF');
    
    t = (0:length(tfY)-1)/60;
    subplot(3,1,3);
    hold on
    plot(t,tfY/norm(tfY),'Color',[1,0.5,0],'LineWidth',2);
    plot(t,tfB/norm(tfB),'Color',[0,0,1],'LineWidth',2);
    hold off; 
    xlabel('time (s)'); ylabel('weight');
    title(['Cluster: ',num2str(jj)]);
    
    if jj == 1
        print(f, [outdir,'20220420C_STA'], '-dpsc');
    else
        print(f, [outdir,'20220420C_STA'], '-dpsc', '-append');
    end
end


function [sfilter,tfilter] = getStaSummary(sta)
% Find the highest variance points.
    V = zeros(size(sta,1:2));
    for o = 1 : size(V,1)
        for p = 1 : size(V,2)
            V(o,p) = var(squeeze(sta(o,p,:)));
        end
    end
    
    vtmp = V;
    vtmp(V < mean(V(:)) + 4*std(V,[],'all')) = 0;
    tf = zeros(size(sta,3),1);
    for o = 1 : size(V,1)
        for p = 1 : size(V,2)
            tf = tf + vtmp(o,p)*squeeze(sta(o,p,:));
        end
    end
    
    [row,col] = find(V > mean(V(:)) + 4*std(V,[],'all'));
    tfilter = squeeze(mean(sta(row,col,:),1:2));
    
    tfilter = tfilter(:); % Make sure it's a column vector.

    tfilter = tfilter / norm(tfilter);
    tfilter = tfilter / max(abs(tfilter));

    % Take the inner product with the STA to get the space filter.
    sfilter = zeros(size(sta,[1,2]));
    for p = 1 : size(sta,1)
        for q = 1 : size(sta,2)
            sfilter(p,q) = tfilter' * squeeze(sta(p,q,:));
        end
    end
end



