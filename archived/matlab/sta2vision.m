

% Reverse the time dimension for Vision.
sta = sta(:,end:-1:1,:,:,:);
sta = permute(sta,[1,3,4,2,5]);

% STA params
params.version = 32; % don't change
params.max = 10000; % don't change
params.stix_size = 66; % stixel size in microns
params.refresh = 16.6550; % frame refresh in ms

if isa(sta, 'cell')
else
    n_cells = length(cluster_id); %size(sta,1); % number of cells
    params.width = size(sta,3); % STA width in STIXELS
    params.height = size(sta,2); % STA height in STIXELS
    params.depth = size(sta,4); % How many frames back per STA
end

neuron_list = cluster_id; %1:n_cells;

% formula for how much space each STA will take (in bytes)
one_stasize = 8 + 4 + params.depth * (4 + 4 + 8 + params.width * params.height * 3 * 2 * 4);

% open file
% fileID = fopen('/Volumes/DataDisk/data/2015-09-23-7/kilosort_data007/data007/data007_test.sta','w');
% fileID = fopen('/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Sort/20220426C/data066/kilosort2/data066.sta','w');
% fileID = fopen('/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Sort/20220420C/data056/kilosort2/data056.sta','w');
fileID = fopen('/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/20220518/noise_by/noise_by.sta','w');

% write header (total 164 bytes)
fwrite(fileID,params.version,'int32','ieee-be');
fwrite(fileID,n_cells,'int32','ieee-be');
fwrite(fileID,params.width,'int32','ieee-be');
fwrite(fileID,params.height,'int32','ieee-be');
fwrite(fileID,params.depth,'int32','ieee-be');

fwrite(fileID,params.stix_size,'real*8','ieee-be');
fwrite(fileID,params.refresh,'real*8','ieee-be');



% identify first location 
location_list(1) = n_cells*4+n_cells*8+164;
% write first neuron ID and location after skipping 164 bytes from the
% beginning of the file (important!)
fseek(fileID, 0, 'bof')
fwrite(fileID,neuron_list(1),'int32',164,'ieee-be');
fwrite(fileID,location_list(1),'int64','ieee-be');

% write out the rest of Neuron IDs and locations
for cnt = 2:n_cells
    fwrite(fileID,neuron_list(cnt),'int32','ieee-be');
    fwrite(fileID,location_list(1)+one_stasize*(cnt-1),'int64','ieee-be');
end

% write each STA
for kk = 1:n_cells
    
    % i store STAs in a cell array, each cell is width x height x 3 x depth
    if isa(sta,'cell')
        mySTA = STA_array{kk};
    else
        mySTA = squeeze(sta(kk,:,:,:,:));
        mySTA = permute(mySTA,[3,1,2,4]);
        % Scale
        mySTA = mySTA - repmat(min(mySTA,[],1:3),[size(mySTA,1:3),1]);
        mySTA = mySTA ./ repmat(max(mySTA,[],1:3),[size(mySTA,1:3),1]);
    end
    

    % for each STA, firsth write refresh and depth
    fwrite(fileID,params.refresh,'real*8','ieee-be');
    fwrite(fileID,params.depth,'int32','ieee-be');
    
    % go through frames one by one
    for j = 1:params.depth
        
        % for each frame, first write width, height, stixel size
        fwrite(fileID,params.width,'int32','ieee-be');
        fwrite(fileID,params.height,'int32','ieee-be');
        fwrite(fileID,params.stix_size,'real*8','ieee-be');
        
        one_frame = squeeze(mySTA(j,:,:,:));
        one_frame = permute(one_frame, [3, 2, 1]);
        one_frame = one_frame(:);
        one_frame = [one_frame, zeros(size(one_frame))]';
        one_frame = one_frame(:);
        % then write Values  for this frame (both real and error ones)
        fwrite(fileID,one_frame,'real*4','ieee-be'); % real frame
        
        % then write Values  for this frame (both real and error ones)
        
%         cc = 1;cr = 1; % row and column indices (width and height)
%         for i = 1:params.width*params.height
%             
%             % here, each STA should have 3 color frames (identical in case
%             % of BW); after each real STA frame, Vision expects an Error
%             % frame that it does not really use so I just input random
%             % numbers there. Stupidly it doubles the size of the file, but
%             % what can you do beyond re-writing Vision?
%             
% 
%             fwrite(fileID,mySTA(j,cc,cr,1),'real*4','ieee-be'); % real frame
%             fwrite(fileID,0.0,'real*4','ieee-be'); % error frame
%             fwrite(fileID,mySTA(j,cc,cr,2),'real*4','ieee-be'); % real frame
%             fwrite(fileID,0.0,'real*4','ieee-be'); % error frame
%             fwrite(fileID,mySTA(j,cc,cr,3),'real*4','ieee-be'); % real frame
%             fwrite(fileID,0.0,'real*4','ieee-be'); % error frame
%             
%             cr = cr+1;
%             % reset width index 
%             if cr==params.width+1
%                 cc = cc+1;
%                 cr = 1;
%             end
%         end
        
    end
    
end

fclose(fileID);