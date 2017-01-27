function [eq_inds] = Equalize_Duration_Distribution(ERPs, Data_Tag, labels1, labels2)
%% This function takes the ECoG table, as well as two data Booleans and Returns 
% Returns of 2 x n vector whose colums represent indecies for 
% ERPs such that the distribution of durations of those two ERPs are
% equivalent, given some binning width

%% set equalizer parameter (handshapes or locations or...)
durations = ERPs.annot.dur_samp;
bin_width = 25; %50 ms wide bin for durations being called the same
dur_bound = 500; %500 ms upper limit on ERPs (used in the case of fs)
num_bins = ceil(dur_bound/bin_width); % 

%% restrict parameter and labels to Data_Tag erps by setting other values to false
%eq_param = eq_param(Data_Tag);
labels1(~Data_Tag) = false;
labels2(~Data_Tag) = false;

durations1 = durations(labels1);
durations2 = durations(labels2);

ind1 = [];
ind2 = [];

%% loop through bins an find equalized listing of each labeling
for i = 1:num_bins
    min_dur = (i-1)*bin_width;
    max_dur = i*bin_width - 1;
    
    %Get list of whether events are in 
    in_bin = (durations < max_dur) & (durations > min_dur);
    
    bin_inds1 = find(in_bin & labels1);
    bin_inds2 = find(in_bin & labels2);
    
    eq_num = min(length(bin_inds1), length(bin_inds2));
    
    % Randomize and Select eqivalent number;
    shuffle = randperm(length(bin_inds1));
    bin_inds1 = bin_inds1(shuffle(1:eq_num));
    shuffle = randperm(length(bin_inds2));
    bin_inds2 = bin_inds2(shuffle(1:eq_num));

    ind1 = [ind1; bin_inds1];
    ind2 = [ind2; bin_inds2];
end
eq_inds= [ind1, ind2];


