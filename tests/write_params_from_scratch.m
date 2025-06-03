% PARAMs

disp('Params')
dummy_params_file = ['/Volumes/Scratch/Users/alexth/dummies/data999.old_params'];
new_params_file = [new_path,'data999.params'];

old_fid = fopen(dummy_params_file, 'r');
new_fid = fopen(new_params_file, 'w');
[fid,msg] = fopen(new_params_file,'w');

nParams = 61;
nNeurons = ncells;
maxNeurons = 10000;

fwrite(new_fid,nParams,'int32','ieee-be');
fwrite(new_fid,nNeurons,'int32','ieee-be');
fwrite(new_fid,maxNeurons,'int32','ieee-be');

fseek(old_fid, 0, 'bof');
fread(old_fid, 3,'int32', 'b');
kk = [];
s = [];
for i = 1:nParams*21-12
    k = fread(old_fid, 1,'int8', 'b');
    kk = [kk k];
    s = [s char(k)];
end

fwrite(new_fid, kk,'int8','ieee-be');


ends = [regexp(s,'\w\W') length(s)];
starts = regexp(s,'\W\w')+1;
params_list = {};
for i = 1:length(starts)
    params_list{i} = s(starts(i):ends(i));
end
params_list = reshape(params_list,[2 length(params_list)/2])';

totalN = nParams*nNeurons;

fseek(old_fid, nParams*21-1, 'bof');
locs = fread(old_fid, totalN, 'int32', 'b')+2;

fseek(old_fid, locs(2+0*nParams), 'bof');
cellNameString = fread(old_fid,7,'int8','ieee-be');

my_sizes = diff(locs(1:62));
my_sizes(9:10) = 6;


clear new_locs
new_locs(1) = 2441282;
for i = 1:61
    new_locs(i+1) = new_locs(i)+my_sizes(i);
end
new_locs = new_locs';
for cellID = 1:n_cells-1
    for i = 1:61
        new_locs(cellID*61+i+1) = new_locs(cellID*61+i)+my_sizes(i);
    end
end
new_locs(end) = [];

fseek(new_fid, nParams*21-1, 'bof');
fwrite(new_fid,new_locs-2,'int32','ieee-be');

end_of_fid = new_locs(end)+8;
a = end_of_fid-ftell(new_fid);
fwrite(new_fid,zeros(1,a),'int8','ieee-be');

% string: 327
% array: varies 319 - not empty, 260 - empty
% double: 200

for cellID = 0:n_cells-1
%     cellID
    for i = 1:nParams
        fseek(new_fid, new_locs(i+cellID*nParams)-2, 'bof');
        ftell(new_fid);
        
        
        if  strcmp(params_list{i,2}, 'String') % string
            fwrite(new_fid,327,'int16','ieee-be');
            fwrite(new_fid,cellNameString','int8','ieee-be');
            
        elseif strcmp(params_list{i,2}, 'DoubleArray') % array
            
            arr = 0;
            switch params_list{i,1}
                case 'RedTimeCourse'
                    fwrite(new_fid,319,'int16','ieee-be');
                    fwrite(new_fid,404,'int16','ieee-be');
                    fwrite(new_fid,0,'int16','ieee-be');
                    fwrite(new_fid,50,'int16','ieee-be');
                    arr = 1;
                case 'GreenTimeCourse'
                    fwrite(new_fid,319,'int16','ieee-be');
                    fwrite(new_fid,404,'int16','ieee-be');
                    fwrite(new_fid,0,'int16','ieee-be');
                    fwrite(new_fid,50,'int16','ieee-be');
                    arr = 1;
                case 'BlueTimeCourse'
                    fwrite(new_fid,319,'int16','ieee-be');
                    fwrite(new_fid,404,'int16','ieee-be');
                    fwrite(new_fid,0,'int16','ieee-be');
                    fwrite(new_fid,50,'int16','ieee-be');
                    arr = 1;
                case 'Auto'
                    fwrite(new_fid,319,'int16','ieee-be');
                    fwrite(new_fid,1604,'int16','ieee-be');
                    fwrite(new_fid,0,'int16','ieee-be');
                    fwrite(new_fid,200,'int16','ieee-be');
                    arr = 1;
                otherwise
                    
                    fwrite(new_fid,260,'int16','ieee-be');
                    fwrite(new_fid,0,'int16','ieee-be');
                    fwrite(new_fid,0,'int16','ieee-be');
            end
            
            if arr==1
                if strcmp(params_list{i,1}, 'GreenTimeCourse')
                    if max(abs(data.time_course{my_unit_id,end}(:,2)))>max(abs(data.time_course{my_unit_id,end}(:,5)))
                        my_array = data.time_course{my_unit_id,end}(:,2);
                    else
                        my_array = data.time_course{my_unit_id,end}(:,5);
                    end
                    
                end
                if strcmp(params_list{i,1}, 'RedTimeCourse')
                    if max(abs(data.time_course{my_unit_id,end}(:,2)))>max(abs(data.time_course{my_unit_id,end}(:,5)))
                        my_array = data.time_course{my_unit_id,end}(:,5);
                    else
                        my_array = data.time_course{my_unit_id,end}(:,2);
                    end
                end
                if strcmp(params_list{i,1}, 'BlueTimeCourse')
                    if max(abs(data.time_course{my_unit_id,end}(:,3)))>max(abs(data.time_course{my_unit_id,end}(:,6)))
                        my_array = data.time_course{my_unit_id,end}(:,3);
                    else
                        my_array = data.time_course{my_unit_id,end}(:,6);
                    end
                end
                if strcmp(params_list{i,1}, 'Auto')
                    my_array = double(data.acf{my_unit_id,end});
                end
                
                fwrite(new_fid,my_array,'real*8','ieee-be');
                
            end
            
        else % just a double
            
            if strcmp(params_list{i,1}, 'ID')
                k = cellID+1;
                visionID = k;
                my_unit_id = k;
                sr = data.sta{my_unit_id,end};
                cont = new_data.contours{my_unit_id};
                if size(cont,1)==1
                    cont = cont{1};
                end
                
            elseif strcmp(params_list{i,1}, 'x0')
                tmp = cont;%data.contours{my_unit_id,1};
                k = mean(tmp(:,1));
                
            elseif strcmp(params_list{i,1}, 'y0')
                tmp = cont;%data.contours{my_unit_id,1};
                k = size(data.sta{1,1},1)-mean(tmp(:,2));
            elseif strcmp(params_list{i,1}, 'SigmaX')
                tmp = cont;%data.contours{my_unit_id,1};
                k = (max(tmp(:,1))-min(tmp(:,1)))/2;
            elseif strcmp(params_list{i,1}, 'SigmaY')
                tmp = cont;%data.contours{my_unit_id,1};
                k = (max(tmp(:,2))-min(tmp(:,2)))/2;
            elseif strcmp(params_list{i,1}, 'Theta')
                k = 0;
            elseif strcmp(params_list{i,1}, 'nSpikes')
                k = double(sum(data.spikes(my_unit_id,:)));
            elseif strcmp(params_list{i,1}, 'tOffset')
                k = 0;
            elseif strcmp(params_list{i,1}, 'blueness') || strcmp(params_list{i,1}, 'amp1') ...
                    || strcmp(params_list{i,1}, 'amp2') || strcmp(params_list{i,1}, 'amp3')
                k = 0;
            elseif strcmp(params_list{i,1}, 'acfBinning')
                k = 0.5;
            elseif strcmp(params_list{i,1}, 'acfMean')
                k = 50.872;
            elseif strcmp(params_list{i,1}, 'acfRMS')
                k = 43.1999;
            else
                k = 0;
            end
            
            fwrite(new_fid,200,'int16','ieee-be');
            fwrite(new_fid,k,'real*8','ieee-be');
            
        end
    end
end

fclose(new_fid);
fclose(old_fid);