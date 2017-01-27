subj = 'CH';
hemi = 'lh';
task = 'ASL';
blocks = [5];
timeLim = [-2 6];
zscore_timeLim = [-1 -0.5];
rootdir = '/Users/changlab/Documents';

extract_ERPs_flag = 1;
plot_ERPs_flag = 0;

plot_conditions = {'lexical','transitional'};
%% load annotations and set up ERPs

% set up annot
fprintf('Creating annotation structure\n');




%% load audio
    audio_data = audioread([rootdir '/data/raw_data/' subj '/CH_B5/Analog/analog2.wav']);
    
%% load master annotation
if ~exist('MasterTable','var');
    fprintf('Loading master table\n');
    MasterTable = readtable(['/Users/changlab/Documents/data/CH/task_info/Audio Timing - Sheet1.csv']);
end

% set up annot
fprintf('Creating annotation structure\n');
annot.notes = table2array(MasterTable(:,9));  % .....................  Place of Contact

annot.start_time = table2array(MasterTable(:,5)); %Given in Seconds
annot.Movement = table2array(MasterTable(:,7));
annot.handshape = table2array(MasterTable(:,8));
annot.full_handshape = table2array(MasterTable(:,9));


%% EXTRACT ERPs

    [ERPs]= Movement_extract_ERPs(rootdir, subj, timeLim, zscore_timeLim, annot);
    ERPs.subj = subj;
    ERPs.hemi = hemi;
    ERPs.task = task;
    ERPs.annot = annot;
    

    

 