function [is_good] = is_good_trial(ERPs)
%% This function returns a Trails length boolean array that is true
% when the Trial is not a bad channel

%% Get a List of Bad Trials
bad_trials = cell2mat(ERPs.BadTrials(:,2));
is_good = true(size(ERPs.ecog,3),1);
is_good(bad_trials) = false;

end