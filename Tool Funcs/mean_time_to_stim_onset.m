function mean_stim_lag = mean_time_to_stim_onset(ERPs, inds)
%% Finds the mean time to the preceding stimulus onset for the indecies in question (in MS)

stim_inds = find(strcmpi(ERPs.annot.stimOnset,'S'));

start_times = ERPs.annot.start_samp/10;
stim_lag = zeros(size(inds));

for i = 1:length(inds)
    closest_stim_ind = max(stim_inds(stim_inds < inds(i)));
    stim_lag(i) = start_times(closest_stim_ind) - start_times(inds(i));
end

mean_stim_lag = mean(stim_lag);