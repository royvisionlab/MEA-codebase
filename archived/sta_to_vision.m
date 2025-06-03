function params = sta_to_vision(filepath, sta, cluster_id, varargin)
% sta_to_vision(filepath, sta, cluster_id, varargin)
% sta: either a cell array with the order [x,y,color, t] 
%    or a matrix: [cell,t,x,y,color]
% cluster_id: cell id for each cluster
% Optional:
% stixel_size: stixel edge length in microns.
% frame_refresh: temporal refresh of sta in msec (e.g. 1000/120 for 120 Hz)


ip = inputParser();
ip.addParameter('version', 32, @(x)isfloat(x));
ip.addParameter('max', 10000, @(x)isfloat(x));
ip.addParameter('stixel_size', 30.0, @(x)isfloat(x)); % Stixel edge length in microns
ip.addParameter('frame_refresh', 1000/120, @(x)isfloat(x)); % Frame refresh

% Parse the inputs.
ip.parse(varargin{:});

% Get the field names from the input parser.
fnames = fieldnames(ip.Results);

% Create the parameters structure.
params = struct();
for jj = 1 : length(fnames)
    params.(fnames{jj}) = ip.Results.(fnames{jj});
end

params.staSize = size(sta);
params.clusterSize = size(cluster_id);

if isa(sta, 'cell')
else
    n_cells = length(cluster_id); %size(sta,1); % number of cells
    params.width = size(sta,4); % STA width in STIXELS
    params.height = size(sta,3); % STA height in STIXELS
    params.depth = size(sta,2); % How many frames back per STA
    
    % Reverse the time dimension for Vision.
    sta = sta(:,end:-1:1,:,:,:);
    sta = permute(sta,[1,3,4,5,2]);
    if min(sta(:)) > -0.001
        sta = 0.5*sta + 0.5; % Scale between 0-1 for Vision...
    end
        
end

neuron_list = cluster_id; %1:n_cells;

% formula for how much space each STA will take (in bytes)
one_stasize = 8 + 4 + params.depth * (4 + 4 + 8 + params.width * params.height * 3 * 2 * 4);

fileID = fopen(filepath,'wb');

% write header (total 164 bytes)
fwrite(fileID,params.version,'int32','ieee-be');
fwrite(fileID,n_cells,'int32','ieee-be');
fwrite(fileID,params.width,'int32','ieee-be');
fwrite(fileID,params.height,'int32','ieee-be');
fwrite(fileID,params.depth,'int32','ieee-be');

fwrite(fileID,params.stixel_size,'real*8','ieee-be');
fwrite(fileID,params.frame_refresh,'real*8','ieee-be');



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
        mySTA = sta{kk};
    else
        mySTA = squeeze(sta(kk,:,:,:,:));
%         mySTA = permute(mySTA,[2,3,4,1]);
%         mySTA = permute(mySTA,[3,1,2,4]);
        % Scale
%         mySTA = mySTA - repmat(min(mySTA,[],1:3),[size(mySTA,1:3),1]);
%         mySTA = mySTA ./ repmat(max(mySTA,[],1:3),[size(mySTA,1:3),1]);
    end

    % for each STA, firsth write refresh and depth
    fwrite(fileID,params.frame_refresh,'real*8','ieee-be');
    fwrite(fileID,params.depth,'int32','ieee-be');
    
    % go through frames one by one
    for j = 1:params.depth
        
        % for each frame, first write width, height, stixel size
        fwrite(fileID,params.width,'int32','ieee-be');
        fwrite(fileID,params.height,'int32','ieee-be');
        fwrite(fileID,params.stixel_size,'real*8','ieee-be');
        
        one_frame = mySTA(:,:,:,j);
        one_frame = permute(one_frame, [3, 2, 1]);
        one_frame = one_frame(:);
        one_frame = [one_frame, zeros(size(one_frame))]';
        one_frame = one_frame(:);
        % then write Values  for this frame (both real and error ones)
        fwrite(fileID,one_frame,'real*4','ieee-be'); % real frame
    end
    
end

fclose(fileID);


