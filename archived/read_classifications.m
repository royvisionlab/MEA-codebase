

filepath = '/Users/michaelmanookin/Documents/Data/rawdata/MEA_Data/Analysis/20220531C/noise/noise.classification.txt';

clusterIds = [];
cell_class = {};

fileId = fopen(filepath,'rt');

tline = fgetl(fileId);
while ischar(tline)
    parts = strsplit(tline,'  ');
    clusterIds = [clusterIds, str2double(parts{1})]; %#ok<AGROW>
    cell_class = [cell_class; parts{2}]; %#ok<AGROW>
    
    tline = fgetl(fileId);
end

fclose(fileId);

targetIds = clusterIds(contains(cell_class,'bt','IgnoreCase',true));
targetIds = clusterIds(contains(cell_class,'SmoothMonoOn','IgnoreCase',true));



