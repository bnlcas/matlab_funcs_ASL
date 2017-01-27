subj = 'CH';
hemi = 'lh';
task = 'ASL';
blocks = [1 2 3 4];
timeLim = [-2 2];
%%*********************************************
zscore_timeLim = [-1 -0.5];
%zscore_timeLim = [-2 -1.5];
%%*********************************************
rootdir = '/Users/changlab/Documents';

extract_ERPs_flag = 1;
plot_ERPs_flag = 0;

plot_conditions = {'lexical','transitional'};
%% load annotations and set up ERPs

% load master annotation
if ~exist('masterTable','var');
    fprintf('Loading master table\n');
    masterTable = readtable([rootdir '/data/' subj '/task_info/ECOG_contrasts_v2.xlsx']);
end

% set up annot
fprintf('Creating annotation structure\n');
annot.trl = table2array(masterTable(:,2));  % .....................  Trial Number
annot.log = table2array(masterTable(:,3));
annot.start_samp = table2array(masterTable(:,4));
annot.start_ms = ASLtoECOG_Time_Convert(annot,annot.start_samp);
annot.end_samp = table2array(masterTable(:,6));
annot.end_ms = ASLtoECOG_Time_Convert(annot,annot.end_samp);
annot.dur_samp = table2array(masterTable(:,10));
annot.dur_ms = ASLtoECOG_Time_Convert(annot,annot.dur_samp);
annot.stimOnset = table2array(masterTable(:,11));
annot.movOnset = table2array(masterTable(:,12));
annot.lexTrans = strrep(table2array(masterTable(:,13)),' ','');
annot.filledLexTrans = strrep(table2array(masterTable(:,14)),' ','');
annot.vidName = table2array(masterTable(:,15));
annot.stimType = table2array(masterTable(:,16));
annot.respType = strrep(table2array(masterTable(:,17)),' ','');
annot.gloss = table2array(masterTable(:,18));
annot.signType = strrep(table2array(masterTable(:,19)),' ','');
annot.intMov = strrep(table2array(masterTable(:,20)),' ','');
annot.handshape = strrep(table2array(masterTable(:,21)),' ','');
annot.loc = strrep(table2array(masterTable(:,22)),' ','');
annot.movPath = strrep(table2array(masterTable(:,23)),' ','');
annot.movDir = strrep(table2array(masterTable(:,24)),' ','');
annot.contact = strrep(table2array(masterTable(:,25)),' ','');

%Add event data (word frequency, age of acquistion, condition_number) to
%annot structure
annot.frequency = evnt_to_master(rootdir, subj, annot.stimOnset, 'freq.nat');
annot.aoa = evnt_to_master(rootdir, subj, annot.stimOnset, 'aoa.nat');
annot.cond_num = evnt_to_master(rootdir, subj, annot.stimOnset, 'cond_num');
annot.word_name = evnt_to_master(rootdir, subj, annot.stimOnset, 'word_label_short');

% Add English Freqency
word_table = load('/Users/changlab/Documents/changrepo/matlab/analysis/ASL/Linguistic Analysis/Word_Stats_table.mat');
word_table = word_table.word_table; % Wonkedy...
word_list = table2cell(word_table(:,1));
freq_eng = zeros(size(annot.frequency));
for i = 1:length(word_list)
    freq_eng(strcmpi(annot.word_name, word_list(i))) = word_table.Freq_Eng(i);
end
freq_eng(freq_eng~=0) = log(freq_eng(freq_eng~=0)); % Take the Log of Frequency
annot.frequency_eng = freq_eng;
    
% Blank NaNs
annot.frequency(isnan(annot.frequency)) = 0;
annot.aoa(isnan(annot.aoa)) = 0;

% Fill in Frequency & AoA for a lexical sign
i = 1;
while i < length(annot.frequency)
    if annot.frequency(i) ~= 0
        i = i+1;
        while ~strcmpi(annot.stimOnset(i),'S')
            annot.frequency(i) = annot.frequency(i-1);
            annot.aoa(i) = annot.aoa(i-1);
            annot.frequency_eng(i) = annot.frequency_eng(i-1);
            i = i+1;
        end
    else
        i = i+1;
    end
end




%% EXTRACT ERPs

if extract_ERPs_flag
    [ERPs]= ASL_extract_ERPs(rootdir, subj, timeLim, zscore_timeLim, annot);
    ERPs.subj = subj;
    ERPs.hemi = hemi;
    ERPs.task = task;
    ERPs.annot = annot;
    
    %% Add a filled in Annotation such that empty spaces inherit the Annotation of the predecessors:
    %ERPs.filled_annot = annot_fill(annot); % Fills empty rows of annotations with preceding values
    %ERPs.grouped_annot =  group_annot_alt(annot); % longest duration of each class in each linguistic gesture
    %ERPS.temp_grouped_annot = annot_group_time_prox(annot); % Groups Tags that occur within 50 ms start time
    %ERPs.filled_change_annot = annot_fill_tochange(annot); % Groups Tags by filling in except in changing
    ERPs.grouped_fill_annot = annot_group_fill(annot);
    
    save([rootdir '/data/' ERPs.subj '/analysis/' ERPs.subj '_ERPs.mat'],'ERPs','-v7.3');
    
    clearvars -except ERPs
    
    Data_Tag = is_good_trial(ERPs) & strcmpi(ERPs.annot.filledLexTrans,'lexical'); % Restricting Boolean true for good linguistic ERPs

else
    load([rootdir '/data/' ERPs.subj '/analysis/' ERPs.subj '_ERPs.mat']);
end

%%
