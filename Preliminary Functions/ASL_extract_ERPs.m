%% Load ECoG data to create ERPs
% Written by Matt Leonard
% Changed by Max on 9/8/2014
% Load ECoG data, look for time of stimuli, crops out the Ecog Data at the time of stimuli,
% Looks for trials with artifacts (calls Find_bad_trials function), crop them out, looks for bad channels, report them but does not crop them out, Z-scores the data by electrode. 
% step1: get all trials in a structure
% step 2: oraganize structure by stimulus

% Update
% - Max 11/5/2014: bug with erasing BadTrials, changes indexing and
% therefore just replaced Ecog by NaN for BadTrials. 

%Update ii
% - Ben Lucas 6/11/2015: Could not find section replacing badTrails in Ecog
% with NaN. Since this will lead to issues with some functions, I have
% created a new cell (ecog_trim) with these sections removed altogether.

function [ERPs]= ASL_extract_ERPs(datarootdir, subj, timeLim, zscore_timeLim, annot)


%% Usage
% datarootdir='/Users/Max/Documents/MATLAB/EC_data/Phoneme_restoration';
% subj='EC63';


%% Load raw data

% calls structure where all events (600) were identified by the
% DetectEventsQuick function. This file is stored in the analysis pathway

%***THIS SECTION HAS BEEN MODIFIED FOR TAKING START TIMES FROM A TABLE NOT
%INCLUDED IN THE RAW DATA*****************************************
% allEventTimes=load([datarootdir '/' subj '/evnt.mat']); %phonrest_evnt.mat']);
% allEventTimes = allEventTimes.evnt;
for i = 1:length(annot.trl)
    allEventTimes(i).name = annot.vidName{i};
    allEventTimes(i).StartTime = annot.start_ms(i);
    allEventTimes(i).block = annot.log{i};
end
allEventStartTimes = annot.start_ms; %[allEventTimes.StartTime];
blocks=annot.log; %{allEventTimes.block};
current_block=blocks{1};
block_length=sum(strcmp(blocks,current_block));number_of_blocks=ceil(length(blocks)/block_length);
fprintf('Starting analysis for %d blocks \n',number_of_blocks);

% reads htk file containining Hilbert transforms in 8 frequency bands as well as sampling frequency (resampled to 400 during preprocessing)
% output is mean High-gamma for each 256 electrodes (columns) over time (rows)
% Load first block (current block)
[ecogDS.data,ecogDS.sampFreq] = readhtks([datarootdir '/data/raw_data/' subj '/' subj '_' current_block '/HilbAA_70to150_8band']);
        ecogDS.data = resample(ecogDS.data,1,4);% resampling 400 Hz sampled Hilbert transforms to 100 Hz for future data handling
        ecogDS.data = ecogDS.data';
        ecogDS.sampFreq = ecogDS.sampFreq/4;

% Define timeframe for analysis
if isempty(timeLim)
    timeLim=[-.500 1.50];% take 500 ms before and 1500 ms after event.
end


%Can use a baseline block for normalization (Matt)
% if ~isempty(bsln_block)
%     [bsln.data,bsln.sampFreq] = readhtks([rootdir '/raw_data/' subj '/' subj '_b' num2str(bsln_block) '/' hilb_band]);
%     bsln.data = resample(bsln.data,1,4);
%     bsln.data = bsln.data';
%     bsln.sampFreq = bsln.sampFreq/4;
% end


%% Preallocation
% 3D matrix with 256 raws, each raw being ecog data from one channel  for
% interval around stimulus (-500 ms to +500 ms sampled at 100 Hz will gives
% 100 data points, layers are different trials
data.ecog = zeros(size(ecogDS.data,1),length(timeLim(1)*ecogDS.sampFreq:timeLim(2)*ecogDS.sampFreq),length(allEventStartTimes));

% set number of blocks to at least 1
k=1;l=0;
BadTrials=cell(1,2);
BadChans = cell(number_of_blocks,2);
%% Loops 
% loops all 600 trials (i) and all 256 channels (j) to  get 1 second Ecog
% around stimuli

            
Is_Stim = strcmpi(annot.stimOnset, 'S'); %% Used to EXTRACT PreStim Normalized ERPs
%Is_Stim = true(size(annot.stimOnset));

for i = 1:length(allEventStartTimes)% typically 600 trials for phoneme restoration task
    %StimStartTime=round(allEventStartTimes(i));
    
    StimStartTime = round(allEventStartTimes(i)*ecogDS.sampFreq);
    datTimes = [round(StimStartTime + timeLim(1)*ecogDS.sampFreq),   round(StimStartTime + timeLim(2)*ecogDS.sampFreq)];% window around stim onset
    
    
    %% ***** Used IN PRESTIM NORM *****
    if Is_Stim(i)
        StartTime_Norm = round(allEventStartTimes(i)*ecogDS.sampFreq);
        datTimes_Norm = [round(StimStartTime + timeLim(1)*ecogDS.sampFreq),   round(StimStartTime + timeLim(2)*ecogDS.sampFreq)];% window around stim onset
    end
    %% *********
    
    
    
    %%
        % looks at transition from one block to the next
        if or(strcmp(allEventTimes(1,i).block,current_block)==0,i==length(allEventTimes))% if the block has changed for this trial or we are at the end of allEventTimes
        
        % Save information  collected during prior block
        badChans = textread([datarootdir '/data/raw_data/' subj '/' subj '_' current_block '/Artifacts/badChannels.txt']);% Calls list of channels with artifacts identified during prelim Analysis (usually flat) 
        BadChans{k,1}=current_block;BadChans{k,2}=badChans;
        
        % BadTrials
        badTrials = Find_bad_trials(datarootdir,subj,current_block,[allEventTimes(i-block_length:i-1).StartTime],timeLim);% calls Find bad trials function
        fprintf('Found %d bad trial(s) in current block \n',length(badTrials));
            if ~isempty(badTrials) %If there are bad trials, put them in an array with their corresponding stimulus
                badTrials=badTrials+i-(block_length+1);% to account for blocks of past trials
                for m=1:length(badTrials)
                BadTrials{m+l,1}= allEventTimes(1,badTrials(m)).name;
                BadTrials{m+l,2}= badTrials(m);
                BadTrials{m+l,3}= current_block;
                end
                l=l+length(badTrials);
            end
     
        k=k+1;
 
            if i<length(allEventTimes) % if this is not the last trial
                current_block=allEventTimes(1,i).block;% change current block to block after transition
                [ecogDS.data,ecogDS.sampFreq] = readhtks([datarootdir '/data/raw_data/' subj '/' subj '_' current_block '/HilbAA_70to150_8band']);
                % load EcoG for next block
                ecogDS.data = resample(ecogDS.data,1,4);% resampling 400 Hz sampled Hilbert transforms to 100 Hz for future data handling
                ecogDS.data = ecogDS.data';ecogDS.sampFreq = ecogDS.sampFreq/4;
            end
        end
        
        
    %% Normalize and Collect data
        Raw_data=ecogDS.data(:,datTimes(1):datTimes(2));
%         Normalizing_data=Raw_data(:,1:round(abs(timeLim(1)/2)*ecogDS.sampFreq)); %takes average of the timewindow taken prior to the stimulus onset
    %    Normalizing_data=Raw_data(:,round(abs(zscore_timeLim(1))*ecogDS.sampFreq):round(abs(zscore_timeLim(1))*ecogDS.sampFreq)+round(abs(zscore_timeLim(2))*ecogDS.sampFreq));
    %% ***** Used IN PRESTIM NORM *****
    if Is_Stim(i)
         Normalizing_data=Raw_data(:,round(abs(zscore_timeLim(1))*ecogDS.sampFreq):round(abs(zscore_timeLim(1))*ecogDS.sampFreq)+round(abs(zscore_timeLim(2))*ecogDS.sampFreq));
    end
    %% ********* 

            
        Normalized_data=gdivide(gsubtract(Raw_data,mean(Normalizing_data,2)),std(Normalizing_data,[],2)); %gsubtract and gdivide treat each row (channel) individually
       
        % Background Subtract:
%         Baseline=Normalized_data(:,round(abs(zscore_timeLim(1))*ecogDS.sampFreq):round(abs(zscore_timeLim(1))*ecogDS.sampFreq)+round(abs(zscore_timeLim(2))*ecogDS.sampFreq));
%         Normalized_data=gsubtract(Normalized_data,mean(Baseline,2)); %gsubtract and gdivide treat each row (channel) individually

        data.ecog(:,:,i) = Normalized_data;
        % dim=2 set mean and std over rows and not columns, [] set n-1 as
        % normalizing factor for z-score
        %(row?mean(row))./std(row) [] set flag to 0 and uses n-1 for std deviation
        data.stims{i}=allEventTimes(i).name;

end

%% TRIM BAD TRIALS
% bad_trials = cell2mat(BadTrials(:,2));
% trim=data.ecog;
% for i = 1:length(bad_trials)
%     trim(:,:,(bad_trials(i) - (i-1))) = [];    % deletes frames with indexes in bad_trials
% end                                                 % adjusts by number of frames deleted (i)
% data.ecog_trim = trim;   


%% COLLECT INFORMATION INTO OUTPUT STRUCT
%         data.Blocks = blocks;
        data.BadChans = BadChans;
        data.BadTrials = BadTrials;
        data.time_axis = timeLim(1)*1000:10:timeLim(2)*1000;    % This assumes timeLim to be in seconds
        data.stim_onset= abs(timeLim(1));
        data.fs=100;
        
        % add Subsets for Word Frequnecy, Condition Number and Age of
        % Acquistion from the LogFile_data
%         data.cond_num = logfile_data(1,:);
%         data.freq = logfile_data(2,:);
%         data.aoa = logfile_data(3,:);
        
        
        ERPs=data;

end

