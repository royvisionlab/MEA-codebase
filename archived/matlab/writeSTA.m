% STA file binary format:
% header 164 bytes
% location list (neuron ID, location)
% for each STA, first write refresh and depth, then frames. For each frame,
% first write Width, Height, and Stixel size, then values. Values for each
% frame include Real Values for color channel, Error values for color
% channel.
​
​
​
​
% STA params
sta.version = 32; % don't change
sta.max = 10000; % don't change
sta.width = 80; % STA width in STIXELS
sta.height = 40; % STA height in STIXELS
sta.depth = 50; % How many frames back per STA
sta.stix_size = 46.4000; % stixel size in microns
sta.refresh = 16.6550; % frame refresh in ms
​
n_cells = 1205; % number of cells
​
​
​
neuron_list = 1:n_cells;
​
% formula for how much space each STA will take (in bytes)
one_stasize = 8 + 4 + sta.depth * (4 + 4 + 8 + sta.width * sta.height * 3 * 2 * 4);
​
% open file
fileID = fopen('/Volumes/DataDisk/data/2015-09-23-7/kilosort_data007/data007/data007_test.sta','w');
​
% write header (total 164 bytes)
fwrite(fileID,sta.version,'int32','ieee-be');
fwrite(fileID,n_cells,'int32','ieee-be');
fwrite(fileID,sta.width,'int32','ieee-be');
fwrite(fileID,sta.height,'int32','ieee-be');
fwrite(fileID,sta.depth,'int32','ieee-be');
​
fwrite(fileID,sta.stix_size,'real*8','ieee-be');
fwrite(fileID,sta.refresh,'real*8','ieee-be');
​
​
​
% identify first location 
location_list(1) = n_cells*4+n_cells*8+164;
% write first neuron ID and location after skipping 164 bytes from the
% beginning of the file (important!)
fseek(fileID, 0, 'bof')
fwrite(fileID,neuron_list(1),'int32',164,'ieee-be');
fwrite(fileID,location_list(1),'int64','ieee-be');
​
% write out the rest of Neuron IDs and locations
for cnt = 2:n_cells
    fwrite(fileID,neuron_list(cnt),'int32','ieee-be');
    fwrite(fileID,location_list(1)+one_stasize*(cnt-1),'int64','ieee-be');
end
​
% write each STA
for kk = 1:n_cells
    
    % i store STAs in a cell array, each cell is width x height x 3 x depth
    mySTA = STA_array{kk};
​
    % for each STA, firsth write refresh and depth
    fwrite(fileID,sta.refresh,'real*8','ieee-be');
    fwrite(fileID,sta.depth,'int32','ieee-be');
    
    % go through frames one by one
    for j = 1:sta.depth
        
        % for each frame, first write width, height, stixel size
        fwrite(fileID,sta.width,'int32','ieee-be');
        fwrite(fileID,sta.height,'int32','ieee-be');
        fwrite(fileID,sta.stix_size,'real*8','ieee-be');
        
        
        % then write Values  for this frame (both real and error ones)
        
        cc = 1;cr = 1; % row and column indices (width and height)
        for i = 1:sta.width*sta.height
            
            % here, each STA should have 3 color frames (identical in case
            % of BW); after each real STA frame, Vision expects an Error
            % frame that it does not really use so I just input random
            % numbers there. Stupidly it doubles the size of the file, but
            % what can you do beyond re-writing Vision?
            
​
            fwrite(fileID,mySTA(cc,cr,1,j),'real*4','ieee-be'); % real frame
            fwrite(fileID,rand(1)/1000,'real*4','ieee-be'); % error frame
            fwrite(fileID,mySTA(cc,cr,2,j),'real*4','ieee-be'); % real frame
            fwrite(fileID,rand(1)/1000,'real*4','ieee-be'); % error frame
            fwrite(fileID,mySTA(cc,cr,3,j),'real*4','ieee-be'); % real frame
            fwrite(fileID,rand(1)/1000,'real*4','ieee-be'); % error frame
            
            cr = cr+1;
            % reset width index 
            if cr==sta.width+1
                cc = cc+1;
                cr = 1;
            end
        end
        
    end
    
end
​
fclose(fileID);
​