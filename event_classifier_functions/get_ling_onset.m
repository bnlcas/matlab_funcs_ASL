function [is_onset] = get_ling_onset(ERPs)
%% Returns the first ERPs that is timestamped at the beginning of Linguistic Movement
% Note that this will also coincide with the unique erps list of the
% grouped_filled_annot scheme, since that function also ascribes uniqueness
% to the first ERPs among ERPs with equiavlent start times
lexical_onset = strcmpi(ERPs.annot.lexTrans,'lexical');
onset_start_times = ERPs.annot.start_samp(lexical_onset);

%% Find ERPs with the same onset time as the lexical onset
is_onset = false(size(lexical_onset));
for i = 1:length(onset_start_times)
    inds = find((onset_start_times(i) == ERPs.annot.start_samp),1);
    is_onset(inds) = true;
end

a = 1;

end