% params_file = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/20220823C/noise/noise_old.params';
% new_params_file = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/20220823C/noise/noise.params';
% params_file = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/20220531C/20220531C/noise/noise_old.params';
% new_params_file = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/20220531C/20220531C/noise/noise.params';
% params_file = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/20220809C/yass/all/all_old.params';
% new_params_file = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/20220809C/yass/all/all.params';
params_file = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/20220531C/kilosort2/noise/noise_old.params';
new_params_file = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/20220531C/kilosort2/noise/noise.params';

if iscell(cluster_id)
    cluster_id = cell2mat(cluster_id);
end


timeIdx=1:61;

%
old_fid = fopen(params_file);
new_fid = fopen(new_params_file, 'w');

% first we copy the params file... there must be a better way to do this
fseek(old_fid, 0, 'eof');
file_size = ftell(old_fid);
fseek(old_fid, 0, 'bof');
while fseek(new_fid,file_size, 'bof')<0
    fseek(new_fid, 0, 'eof');
    k = fread(old_fid,1,'int8','ieee-be');
    fwrite(new_fid,k,'int8','ieee-be');
end

% read first three intergers in params file 
fseek(old_fid, 0, 'bof');
nParams = fread(old_fid, 1,'int32', 'b');
nNeurons = fread(old_fid, 1,'int32', 'b');
maxNeurons = fread(old_fid, 1,'int32', 'b');

% get params list (first 1500 chars are reserved for parameters description)
kk = [];
s = [];
for i = 1:1500
    k = fread(old_fid, 1,'int8', 'b');
    kk = [kk k];
    s = [s char(k)];
end

total_end= 0;
for i = 1:1500-2
    if kk(i)==0 && kk(i+1)==0 && kk(i+2)==0
        total_end = i+2+find(kk(i+3:end)==0,1);
    end
end
s = s(1:total_end);

% transform the parameters line into a list of "parameter" - "data type"
ends = [regexp(s,'\w\W') length(s)];
starts = regexp(s,'\W\w')+1;
params_list = cell(1,length(starts));
for i = 1:length(starts)
    params_list{i} = s(starts(i):ends(i));
end
params_list = reshape(params_list,[2 length(params_list)/2])';

% total number of entries 
totalN = nParams*nNeurons;


for i = 1:50
    fseek(old_fid, total_end+i, 'bof');
    k = fread(old_fid, 1, 'int32', 'b');
    if k>2000000 && k<3000000
        break
    end
end
fseek(old_fid, total_end+i, 'bof');
% real locations are shifted by 2 bytes from recorded locations, i guess just to fuck with us 
locs = fread(old_fid, totalN, 'int32', 'b')+2;

% number of cells (should be the same in the old file and new file)
n_cells = size(timecourse_matrix,1);

%this loop with read Vision parameters, replace some of them with new stuff
%and write them to the new file
%
for cellID = 0:n_cells-1
    for i = 1:nParams
        fseek(old_fid, locs(i+cellID*nParams), 'bof');
        fseek(new_fid, locs(i+cellID*nParams), 'bof');
        if  strcmp(params_list{i,2}, 'String') % string
            ss = fread(old_fid, 100,'int8', 'b');
            for j = 1:100
                fwrite(new_fid,ss(j),'int8','ieee-be');
            end
            ss = char(ss);
            ss = ss';
            a = regexp(ss,'All');
            b = regexp(ss, '\W');
            c = regexp(ss, '\/');
            b = setdiff(b,c);
            cellClass = ss(a:b(find(b>a,1))-1);
            
            
        elseif strcmp(params_list{i,2}, 'DoubleArray') % array
            for j = 1:3
                k = fread(old_fid, 1,'int16', 'b');
                fwrite(new_fid,k,'int16','ieee-be');
                shift = j-1;
                arr = k;
            end
            
            k = fread(old_fid, 1,'int16', 'b');
            fwrite(new_fid,k,'int16','ieee-be');
            
            fseek(old_fid, locs(i+cellID*nParams)+2+4, 'bof');
            fseek(new_fid, locs(i+cellID*nParams)+2+4, 'bof');
            
            if arr>0 && arr<=10000
                my_array = fread(old_fid, arr,'real*8', 'b');
                if strcmp(params_list{i,1}, 'GreenTimeCourse')
                    my_array = squeeze(timecourse_matrix(cellIdx,timeIdx,2));
                    my_array = fliplr(my_array(:)');
                end
                if strcmp(params_list{i,1}, 'RedTimeCourse')
                    my_array = squeeze(timecourse_matrix(cellIdx,timeIdx,1));
                    my_array = fliplr(my_array(:)');
                end
                if strcmp(params_list{i,1}, 'BlueTimeCourse')
                    my_array = squeeze(timecourse_matrix(cellIdx,timeIdx,3));
                    my_array = fliplr(my_array(:)');
                end
                if strcmp(params_list{i,1}, 'Auto')
%                     my_array = double(acf(cellIdx, 1:length(my_array))); %acf(cellIdx, :);
                end
                
                for j = 1:length(my_array)
                    fwrite(new_fid,my_array(j),'real*8','ieee-be');
                end
            end
            
        else % just a double
            k = fread(old_fid, 1,'real*8', 'b');
            if strcmp(params_list{i,1}, 'ID')
                visionID = k;
                my_unit_id = k;
                cellIdx = find(cluster_id == my_unit_id);
                verts = squeeze(hull_vertices(cellIdx,:,:));
                if all(verts(:)==0)
                    verts = 10*ones(size(verts));
                end
                verts(verts(:,1)==0,:) = [];
                verts = verts / 5; % need to fix this in the python code...
                if isempty(verts)
                    verts = zeros(2,2);
                end
            end
            if strcmp(params_list{i,1}, 'x0')
                k = max(1,gauss_params(cellIdx,3));
            elseif strcmp(params_list{i,1}, 'y0')
                k = size(sta,3)-max(1,gauss_params(cellIdx,4));
                k = max(1,k);
            elseif strcmp(params_list{i,1}, 'SigmaX')
                k = max(1,abs(gauss_params(cellIdx,5)));
            elseif strcmp(params_list{i,1}, 'SigmaY')
                k = max(1,abs(gauss_params(cellIdx,6)));
            elseif strcmp(params_list{i,1}, 'Theta')
                k = gauss_params(cellIdx,2);
            elseif strcmp(params_list{i,1}, 'nSpikes')
%                 k = sum(isi(cellIdx,:));
            elseif strcmp(params_list{i,1}, 'acfBinning')
%                 k = 1.0;
            end
            fwrite(new_fid,k,'real*8','ieee-be');
            
        end
    end
end


fclose(new_fid);
fclose(old_fid);