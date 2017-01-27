function [] = Duration_Timing_BoxPlots(ERPs)
%% Creates box plot for the stimulus on set timings for Graphic #5...

annot = ERPs.annot;

%% Get Stim on Indecies
Stim_On = strcmpi(annot.stimOnset,'S');

start_times = annot.start_samp/10; % start times in (ms)
end_times = annot.end_samp/10; % end times in (ms)

stim_start = start_times(Stim_On);
stim_end = end_times(Stim_On);

% Get Movement Stat Times
Mov_On = strcmpi(annot.movOnset,'M');
[Movement_Tag, Filled_Movement_Tag] = Get_First_Subsequent_Tag(Stim_On, Mov_On);
Stim_start_mov = start_times(Filled_Movement_Tag);
Mov_start = start_times(Movement_Tag);

%% Get Lexical Start Times
Lex_On = strcmpi(annot.lexTrans,'lexical');
[Lex_Tag, Filled_Lex_Tag] = Get_First_Subsequent_Tag(Stim_On, Lex_On);
Stim_start_lex = start_times(Filled_Lex_Tag);
Lex_start = start_times(Lex_Tag);

%% Get Lexical to Movement_Start Times 
Mov_On = Movement_Tag;
[Movement_Lex_Tag, Filled_Movement_Lex_Tag] = Get_First_Subsequent_Tag(Mov_On, Lex_On);


%% Plot Forward ScatterPlot

figure; boxplot([stim_end - stim_start, Mov_start-Stim_start_mov, Lex_start-Stim_start_lex], 'orientation','horizontal', 'outliersize',0)
ax = gca;
ax.YTick = 1:3;
ax.YTickLabels = {'stimulus duration', 'movement reaction time', 'linguistic reaction time'}
xlabel('Time from Stimulus Onset (ms)')
xlim([0 2200])
title('Stimulus Locked Movement Sequence')

%% Plot Backward ScatterPlot
figure; boxplot([start_times(Filled_Movement_Lex_Tag)-start_times(Movement_Lex_Tag), end_times(Filled_Lex_Tag) - Lex_start, start_times(Filled_Lex_Tag) - Lex_start], 'orientation','horizontal', 'outliersize',0)
ax = gca
ax.YTick = 1:3;
ax.YTickLabels = {'Pre-Linguistic Movement duration', 'End of Stimulus', 'Start of Stimulus'}
xlabel('Time from Linguistic Onset (ms)')
xlim([-2200 0])
title('Lexical Movement Locked Event Sequence')