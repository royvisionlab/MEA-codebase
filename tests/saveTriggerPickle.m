function saveTriggerPickle(ttlDir, fileName, triggerTimes, numSamples, arrayId, offset)

if ~exist(ttlDir, 'dir')
    mkdir( ttlDir );
end


pickle_file_path = [ttlDir,fileName,'.p'];

% Format the trigger times.
% triggerString = 'np.array([';
% for jj = 1 : length(triggerTimes)
%     triggerString = [triggerString,num2str(triggerTimes(jj))]; %#ok<AGROW>
%     if jj < length(triggerTimes)
%         triggerString = [triggerString,',']; %#ok<AGROW>
%     end
% end
% triggerString = [triggerString,'])']; 

triggerString = 'ttl = np.array([';
for jj = 1 : length(triggerTimes)
    triggerString = [triggerString,num2str(triggerTimes(jj))]; %#ok<AGROW>
    if jj < length(triggerTimes)
        triggerString = [triggerString,',']; %#ok<AGROW>
    end
end
triggerString = [triggerString,'])']; 


% if ~exist(filename,'file')
%     error('%s is not a file',filename);
% end
% outname = [tempname() '.mat'];
% pyscript = ['import pickle;import sys;import scipy.io;file=open("' filename '", "rb");dat=pickle.load(file);file.close();scipy.io.savemat("' outname '.dat")'];
% system(['python -c "' pyscript '"']);
% a = load(outname);

% pyscript = [ 
%     'import pickle; import numpy as np; trigger_data_dict = {''trigger_times'': ', triggerString ,','...
%     '''array_id'': 504, ''n_samples'': ', num2str(numSamples), ',', ...
%     '''neuron_spike_time_offset'': 0 }; ',...
%     'with open(',pickle_file_path,', ''wb'') as pfile:    pickle.dump(trigger_data_dict, pfile)'
%     ];

% pyscript = [ 
%     'import pickle; import numpy as np; trigger_data_dict = {''trigger_times'': ', num2str(10) ,','...
%     '''array_id'': 504, ''n_samples'': ', num2str(numSamples), ',', ...
%     '''neuron_spike_time_offset'': 0 }; ',...
%     'pickle.dump(trigger_data_dict, open(''',pickle_file_path,''', ''wb''))'
%     ];

% pyscript = [ 
%     'import pickle; import numpy.array as np; trigger_data_dict = {''trigger_times'': ', triggerString ,','...
%     '''array_id'': 504, ''n_samples'': ', num2str(numSamples), ',', ...
%     '''neuron_spike_time_offset'': 0 }; ',...
%     'pickle.dump(trigger_data_dict, open(''',pickle_file_path,''', ''wb''))'
%     ];

% pyscript = [ 
%     'import pickle; import numpy as np;',triggerString,'; trigger_data_dict = {''trigger_times'': ttl,'...
%     '''array_id'': 504, ''n_samples'': ', num2str(numSamples), ',', ...
%     '''neuron_spike_time_offset'': 0 }; ',...
%     'pickle.dump(trigger_data_dict, open(''',pickle_file_path,''', ''wb''))'
%     ];

pyscript = [ 
    'import numpy as np;import export_triggers;',triggerString,'; export_triggers.save_trigger_pickle_file(''',pickle_file_path,''',',num2str(arrayId),',',num2str(numSamples),',ttl,',num2str(offset),')'
    ];

% system(['python -c "' pyscript '"']);
system(['/Users/michaelmanookin/opt/anaconda3/bin/python -c "' pyscript '"']);


% trigger_data_dict = {
%     'trigger_times': trigger_times,
%     'array_id': array_id,
%     'n_samples': n_samples,
%     'neuron_spike_time_offset': neuron_spike_time_offset
% }

% with open(pickle_file_path, 'wb') as pfile:
%     pickle.dump(trigger_data_dict, pfile)



