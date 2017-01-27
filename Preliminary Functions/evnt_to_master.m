function [LogFile_annot] = evnt_to_master(rootdir, subj, stimOnset, logfile_address);
%%This file loads the event matrix for each of the blocks of stimulus
%%recording, and then appends this information into an array of the
%%dimensions of the other fields in annot. The resulting array is 0 except
%%for the rows in which a stimulus is given.


logfile_info = [];      
for i = 1:4    % MAGIC CONSTANT for each of the four sessions
    filename = [rootdir '/data/' subj '/behavioral_data/output_run' num2str(i) '.mat'];
    LogFile = load(filename);
    block_data = eval(['LogFile.logfile.' logfile_address]);
    logfile_info = [logfile_info, block_data];
end
if isnumeric(logfile_info)
    LogFile_annot = zeros(size(stimOnset));
else
    LogFile_annot = repmat({''},size(stimOnset));
end
LogFile_annot(strcmpi(stimOnset, 'S')) = logfile_info;  % Add stimulus information to annotations for stimulus events

end