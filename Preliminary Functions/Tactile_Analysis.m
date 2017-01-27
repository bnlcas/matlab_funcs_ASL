subj = 'CH';
hemi = 'lh';
task = 'ASL';
blocks = [6];
timeLim = [-2 2];
zscore_timeLim = [-1 -0.5];
rootdir = '/Users/changlab/Documents';

extract_ERPs_flag = 1;
plot_ERPs_flag = 0;

plot_conditions = {'lexical','transitional'};
%% load annotations and set up ERPs

% set up annot
fprintf('Creating annotation structure\n');




%% load audio
    audio_data = audioread([rootdir '/data/raw_data/' subj '/CH_B6/Analog/analog2.wav']);
    
%% load master annotation
if ~exist('ContactDataTable','var');
    fprintf('Loading master table\n');
%    ContactDataTable = readtable(['/Users/changlab/Documents/data/CH/task_info/Audio Timing - Sheet1.csv']);
    ContactDataTable = csvimport(['/Users/changlab/Documents/data/CH/task_info/Audio Timing - Sheet1.csv']);
    ContactDataTable = cell2table(ContactDataTable(2:end,:),'VariableNames',{'StimLoc','OnsetTime', 'StimOn','Notes'});
end

% set up annot
fprintf('Creating annotation structure\n');
annot.contact_loc = table2array(ContactDataTable(:,1));  % .....................  Place of Contact

annot.start_ms = table2array(ContactDataTable(:,2)); %Given in Seconds
annot.stimOnset= table2array(ContactDataTable(:,3));



%% EXTRACT ERPs

    [ERPs_tact]= Tactile_extract_ERPs(rootdir, subj, timeLim, zscore_timeLim, annot);
    ERPs_tact.subj = subj;
    ERPs_tact.hemi = hemi;
    ERPs_tact.task = task;
    ERPs_tact.annot = annot;
    

    
%clearvars -except ERPs
 