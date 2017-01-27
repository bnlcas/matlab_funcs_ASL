function [Cleaned] = Remove_Bad_Trials_Ecog(ERPs, ecog_data)
%%This function takes a 256 x timepts x trials matrix of ecog_data from the ERPs
%%and elinimates the trials which are considered bad.

Cleaned = ecog_data;

bad_trials = cell2mat(ERPs.BadTrials(:,2));
removed = 0;
for i = 1:length(bad_trials)
    Cleaned(:,:,bad_trials(i)-removed) = [];
    removed = removed + 1;
end