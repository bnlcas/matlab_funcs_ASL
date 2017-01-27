%%
% subject number
subj = 'CH';
% block number
block = [1 2 3 4];
% whether to epoch to production onset
prod_flag = 0;

% rootdir
rootdir = '/Users/changlab/Documents/data';
% raw data directory
rawdatadir = [rootdir '/raw_data'];
% % logfile directory
logfile_dir = [rootdir '/' subj '/behavioral_data'];

% frequency band
hilb_band = 'HilbAA_70to150_8band'; % HilbAA_70to150_8band, HilbReal_4to200_40band
band_name = 'hg'; % theta, hg, 40band

% time window
tlim = [-1 3];

% force flag
force_flag = 0;

%% Prepare events file
% block

% if ~exist(outEventfile) || force_flag

%     [evntTime,allConditions,evnt] = run_ECogFindEvents_phonotactic(soundfile_dir,rawdatadir,evnt_block,anin_htkfile,outEventfile,outCondfile);
evnt = [];
logfile_data = [];
for block_num = 1:length(block)
    fprintf('Block [%d] of [%d]\n',block_num,length(block));
    evnt_block{1} = [subj '_B' num2str(block(block_num))];
    outEventfile = [rawdatadir '/' subj '/' evnt_block{1} '/Analog/allEventTimes.mat'];
    outCondfile = [rawdatadir '/' subj '/' evnt_block{1} '/Analog/allConditions.mat'];
    anin_file = [rawdatadir '/' subj '/' evnt_block{1} '/Analog/analog1.wav'];
    load([logfile_dir '/output_run' num2str(block(block_num)) '.mat']);
    
    logfile_data = [logfile_data, [logfile.cond_num; logfile.freq.nat; logfile.aoa.nat]];

    evntTime = [];
    if ~prod_flag
        [anin,anin_sf] = audioread(anin_file);
        anin = abs(anin - mean(anin));
        
        thresh = 0.02;
        dead_time = anin_sf;
        i=1;
        j = 1;
        thx_idx = [];
        while i<length(anin)
            if anin(i)>thresh
                thx_idx(j)=i;
                i=i+dead_time;
                j = j + 1;
            end
            i = i+1;
        end
        
        
        %thx = find(diff(anin) > std(diff(anin)));
        %thx_idx = find(diff(thx) > std(diff(thx)));
        evntTime = thx_idx; %[thx(1) ; thx(thx_idx+1)];
        fprintf('Found %d events....\n',length(evntTime));
        
        allConditions = logfile.word_label;
        
        for i = 1:length(allConditions)
            allEventTimes{i,1} = evntTime(i) ./ anin_sf;
            allEventTimes{i,2} = logfile.word_label{i};
            allEventTimes{i,3} = logfile.cond_name{i};
            allEventTimes{i,4} = logfile.cond_num(i);
            tmpevnt(i).word_label = logfile.word_label{i};
            tmpevnt(i).StartTime = evntTime(i) ./ anin_sf;
            tmpevnt(i).block = ['B' num2str(block(block_num))];
        end
        

    else
        fprintf('Found %d events....\n',length(logfile.mvt_onset));
        allConditions = logfile.word_label;
        %Could this for loop be eliminate?
        %allEventTimes{:,1}=logfile.mvt_onset(i);...
        for i = 1:length(allConditions)
            allEventTimes{i,1} = logfile.mvt_onset(i);
            allEventTimes{i,2} = logfile.word_label{i};
            allEventTimes{i,3} = logfile.cond_name{i};
            allEventTimes{i,4} = logfile.cond_num(i);
            % allEventTimes{i,5} = logfile.cond_num(i); 
            % allEventTimes{i,6} = logfile.freq.nat(i);
            % allEventTimes{i,7} = logfile.aoa.nat(i);
        end
        
        allConditions = logfile.word_label;
        
    end
    
    fprintf('%s: Re-saving allEventTimes file....\n',mfilename);
    save(outEventfile,'allEventTimes');
    save(outCondfile,'allConditions');
    % else
    %     fprintf('\n%s: allEventTimes file %s already exists, loading....\n',mfilename,outEventfile);
    %     load(outEventfile);
    %     load(outCondfile);
    % end
    
    evnt = [evnt tmpevnt];
end
%% MAKE ERPs STRUCTURE

Trans_OnDur_Times = Get_TransLexStartTimes(MaterTable, 'transitional');
Lex_OnDur_Times = Get_TransLexStartTimes(MaterTable, 'lexical');

[ERPs_long] = Extract_ERPs_step1(rawdatadir, subj, tlim, logfile_data);
%[ERPs_Trans]= Extract_ERPs_Masterfile_step1(rawdatadir, subj, tlim, logfile_data, Trans_OnDur_Times);
%[ERPs_Lex]= Extract_ERPs_Masterfile_step1(rawdatadir, subj, tlim, logfile_data, Lex_OnDur_Times);

% %% TRIM BAD TRIALS
% bad_trials = cell2mat(ERPs.BadTrials(:,2));
% ERPs_Trim=ERPs.ecog;
% for i = 1:length(bad_trials)
%     ERPs_Trim(:,:,(bad_trials(i) - (i-1))) = [];    % deletes frames with indexes in bad_trials
% end                                                 % adjusts by number of frames deleted (i)
%     

%% PLOT ERPs
plotMeanERPs(ERPs)
